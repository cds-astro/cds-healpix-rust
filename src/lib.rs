//! Implementation of the HEALPix framework.  
//! See papers:  
//! * Gorsky2005: "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
//!                Distributed on the Sphere", Górski, K. M. et al., 2005; 2005ApJ...622..759G.
//! * Calabretta2004: "Mapping on the HEALPix grid", Calabretta, M. R., 2004; 2004astro.ph.12607C
//! * Calabretta2007: "Mapping on the HEALPix grid", Calabretta, M. R. et Roukema, B. F., 2007; 2007MNRAS.381..865C
//! * Reinecke2015: "Efficient data structures for masks on 2D grids",  Reinecke, M. et Hivon, E., 2015; 2015A&A...580A.132R

// #![cfg_attr(test, feature(test))]
// #![cfg_attr(test)]


// #[cfg(test)]
// extern crate test;

use std::sync::Once;
use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};

/// Constant = sqrt(6).
/// 
/// ```rust
/// use cdshealpix::{SQRT6};
/// assert_eq!(6_f64.sqrt(), SQRT6);
/// ```
pub const SQRT6: f64 = 2.449_489_742_783_178_f64;
const ONE_OVER_SQRT6: f64 = 0.408_248_290_463_863_f64;
const HALF: f64 = 0.5_f64;

/// Upper limit on sqrt(3(1-|z|)) to consider that we are not near from the poles
const EPS_POLE: f64 = 1e-13_f64;

/// Constant = 2 * pi.
/// 
/// ```rust
/// use cdshealpix::{TWICE_PI};
/// use std::f64::consts::PI;
/// assert_eq!(2f64 * PI, TWICE_PI);
/// ```
pub const TWICE_PI: f64 = 2.0 * PI;

/// Constant = 4 / pi.
/// 
/// ```rust
/// use cdshealpix::{FOUR_OVER_PI};
/// use std::f64::consts::PI;
/// assert_eq!(4f64 / PI, FOUR_OVER_PI);
/// ```
pub const FOUR_OVER_PI: f64 = 4_f64 / PI;

/// Constant = pi / 4.
///
/// ```rust
/// use cdshealpix::{PI_OVER_FOUR};
/// use std::f64::consts::PI;
/// assert_eq!(PI / 4f64, PI_OVER_FOUR);
/// ```
pub const PI_OVER_FOUR: f64 = 0.25_f64 * PI;

/// Constant = 29, i.e. the largest possible depth we can store on a signed positive long
/// (4 bits for base cells + 2 bits per depth + 2 remaining bits (1 use in the unique notation).
/// 
/// ```rust
/// use cdshealpix::{DEPTH_MAX};
/// assert_eq!(29, DEPTH_MAX);
/// ```
pub const DEPTH_MAX: u8 = 29;

/// Constant = nside(29), i.e. the largest possible nside available when we store HEALPix hash
/// on a u64.
///
/// ```rust
/// use cdshealpix::{DEPTH_MAX, NSIDE_MAX, nside};
/// assert_eq!(nside(DEPTH_MAX), NSIDE_MAX);
/// ```
pub const NSIDE_MAX: u32 = 536870912;

/// Limit on the latitude (in radians) between the equatorial region and the polar caps.
/// Equals asin(2/3) = 0.7297276562269663 radians ~= 41,81 degrees.
/// Written $\theta_X$ in Calabretta2007.
/// 
/// ```rust
/// use cdshealpix::{TRANSITION_LATITUDE};
/// assert_eq!(f64::asin(2f64 / 3f64), TRANSITION_LATITUDE);
/// ```
pub const TRANSITION_LATITUDE: f64 = 0.729_727_656_226_966_3_f64; // asin(2/3)
/// Limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 2/3, see Eq. (1) in Gorsky2005.
pub const TRANSITION_Z: f64 = 2_f64 / 3_f64;
/// Inverse of the limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
/// Equals 1/(2/3) = 1.5, see Eq. (1) in Gorsky2005.
pub const ONE_OVER_TRANSITION_Z: f64 = 1.5_f64;

/// Mask to keep only the f64 sign
pub const F64_SIGN_BIT_MASK: u64 = 0x8000000000000000;
/// Equals !F64_SIGN_BIT_MASK (the inverse of the f64 sign mask)
pub const F64_BUT_SIGN_BIT_MASK: u64 = 0x7FFFFFFFFFFFFFFF;

/// For each HEALPix depth, stores the smallest distance from an edge of a cell to the opposite
/// edge of the same cell. If the radius of a cone is smaller than this distance, we know that
/// it will overlap maximum 9 pixels (the pixel containing the center of the cone plus
/// the 8 neighbours).  
/// In practice, this distance if the distance between the point of coordinate
/// (0, TRANSITION_LATITUDE) and it nearest point on the Northeast edge of the
/// cell of base hash 0 and coordinates in the base hash (x=0, y=nside-1).  
/// IMPORTANT REMARK:  
/// - this value is larger than the smallest center to vertex distance
/// - this value x2 is larger than the smallest diagonal (NS or EW)
/// - this value x2 is larger than the smallest edge
/// - BUT there is no case in which the value is larger than the four center-to-vertex distance
/// - BUT there is no case in which the value x2 is larger than both diagonals
/// => this radius is smaller than the smaller circumcircle radius (=> no cone having the smaller
/// -edge-to-opposite-edge-radius radius can contains the 4 vertices of a cell (but 3 is ok)
/// vertices
static SMALLER_EDGE2OPEDGE_DIST: [f64; 30] =  [
    0.8410686705685088,    // depth = 0
    0.37723631722170053,   // depth = 1
    0.18256386461918295,   // depth = 2
    0.09000432499034523,   // depth = 3
    0.04470553761855741,   // depth = 4
    0.02228115704023076,   // depth = 5
    0.011122977211214961,  // depth = 6
    0.005557125022105058,  // depth = 7
    0.0027774761500209185, // depth = 8
    0.0013884670480328143, // depth = 9
    6.941658374603201E-4,  // depth = 10
    3.4706600585087755E-4, // depth = 11
    1.7352877579970442E-4, // depth = 12
    8.676333125510362E-5,  // depth = 13
    4.338140148342286E-5,  // depth = 14
    2.1690634707822447E-5, // depth = 15
    1.084530084565172E-5,  // depth = 16
    5.422646295795749E-6,  // depth = 17
    2.711322116099695E-6,  // depth = 18
    1.3556608000873442E-6, // depth = 19
    6.778303355805395E-7,  // depth = 20
    3.389151516386149E-7,  // depth = 21
    1.69457571754776E-7,   // depth = 22
    8.472878485272006E-8,  // depth = 23
    4.236439215502565E-8,  // depth = 24
    2.1182195982014308E-8, // depth = 25
    1.0591097960375205E-8, // depth = 26
    5.295548939447981E-9,  // depth = 27
    2.647774429917369E-9,  // depth = 28
    1.3238871881399636E-9  // depth = 29
];

/// Latitude, in the equatorial region, for which the distance from the cell center to its four
/// vertices is almost equal on the sky (i.e. the shape of the cell on the sky is close to a square).
/// The larger the depth, the better the approximation (based on differential calculus).
/// > dX = dY = 1 / nside (center to vertex distance)
/// > X = 4/pi * lon     => dX = 4/pi dlon
/// > Y = 3/2 * sin(lat) => dY = 3/2 * cos(lat) dlat
/// > dlon * cos(lat) = dlat (same distance on the sky)
/// > => cos^2(lat) = 2/3 * 4/pi
/// > => lat = arccos(sqrt(2/3 * 4/pi)) ~= 22.88 deg ~= 0.39934 rad
/// 
/// ```rust
/// use cdshealpix::{TRANSITION_Z, FOUR_OVER_PI, LAT_OF_SQUARE_CELL};
/// assert!(f64::abs(f64::acos(f64::sqrt(TRANSITION_Z * FOUR_OVER_PI)) - LAT_OF_SQUARE_CELL) < 1e-15_f64);
/// ```
pub static LAT_OF_SQUARE_CELL: f64 = 0.399_340_199_478_977_75_f64;
/// Simply the consine of LAT_OF_SQUARE_CELL
static COS_LAT_OF_SQUARE_CELL: f64 = 0.921_317_731_923_561_3_f64;

/// Array storing pre-computed values for each of the 30 possible depth (from 0 to 29)
/// Info: I would have prefered to compute those quantities at compilation time, and thus have
/// a `static CSTS_C2V: [ConstantsC2V; 30]`. Unfortunately:
/// - macro do not seems work with static arrays
/// - const fn is not stable and can only use const fn (so no min/max/sin/ ...) :o/
static mut CSTS_C2V: [Option<ConstantsC2V>; 30] = [
  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
];
// Found here: https://stackoverflow.com/questions/28656387/initialize-a-large-fixed-size-array-with-non-copy-types
// I wanted to use it to set a static array. So far it is not possible.
// I hope in the future to be able to compute such an array at compilation time.
/*macro_rules! make_array {
  ($n: expr, $constructor: expr) => {
    {
      let mut items: [_; $n] = mem::uninitialized();
      for (i, place) in items.iter_mut().enumerate() {
          ptr::write(place, $constructor(i as u8));
      }
      items
    }
  }
}
static CSTS_C2V: [ConstantsC2V; 30] = unsafe { make_array!(30, |depth| new_cst_c2v(depth)) };
*/

/// See the get_or_create function, each object is used for the lazy instantiation of the 
/// layer of the corresponding depht.
/// Info: Unfortunatly Default::default(); do no work with static arrays :o/
static CSTS_C2V_INIT: [Once; 30] = [
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new()
];
/// Lazy factory method: instantiate a new Layer at the first call for a given depth; after the 
/// first call, returns an already instantiated Layer.
/// # Info
/// This method resort to a double-checked lock, ensuring thread-safety.
fn get_or_create(depth: u8) -> &'static ConstantsC2V {
  unsafe {
    // Inspired from the Option get_or_insert_with method, modified to ensure thread safety with
    // https://doc.rust-lang.org/std/sync/struct.Once.html
    // This implements a double-checked lock
    match CSTS_C2V[depth as usize] {
      Some(ref v) =>  return v,
      None => {
        CSTS_C2V_INIT[depth as usize].call_once(|| {
          CSTS_C2V[depth as usize] = Some(ConstantsC2V::new(depth));
        });
      },
    }
    match CSTS_C2V[depth as usize] {
      Some(ref v) => v,
      _ => unreachable!(),
    }
  }
}

struct ConstantsC2V {
  slope_npc: f64,
  intercept_npc: f64,
  slope_eqr: f64,
  intercept_eqr: f64,
  coeff_x2_eqr: f64,
  coeff_cst_eqr: f64,
}

impl ConstantsC2V {
  fn new(depth: u8) -> ConstantsC2V {
    let nside = nside_unsafe(depth);
    let dist_cw = 1.0_f64 / (nside as f64); // Center to West (or East) vertex distance on the transition latitude
    let one_min_dist_cw = 1.0_f64 - dist_cw;
    // NPC, see comment of function largest_c2v_dist_in_npc()
    let lat_north = f64::asin(1_f64 - (pow2(one_min_dist_cw) / 3_f64));
    let d_min_npc = lat_north - TRANSITION_LATITUDE;
    let d_max_npc: f64 = sphe_dist(
      squared_half_segment(
        FRAC_PI_4 * dist_cw, d_min_npc,
        lat_north.cos(),  TRANSITION_LATITUDE.cos()));
    debug_assert!(d_min_npc < d_max_npc);
    // - linear approx
    let slope_npc: f64 = (d_max_npc - d_min_npc) / (FRAC_PI_4 * one_min_dist_cw); // a = (yB - yA) / (xB - xA); with xA = 0
    let intercept_npc: f64 = d_min_npc;                                           // b = yA - a * xA          ; with xA = 0
    // EQR TOP, see comment of function largest_c2v_dist_in_eqr_top()
    let d_min_eqr_top = PI_OVER_FOUR * dist_cw * COS_LAT_OF_SQUARE_CELL; // pi/4 * 1/nside * cos(LAT_OF_SQUARE_CELL)
    let d_max_eqr_top = TRANSITION_LATITUDE - f64::asin(one_min_dist_cw * TRANSITION_Z);
    //debug_assert!((d_min_npc - d_max_eqr_top).abs() < 1e-5, "depth: {}, 1: {}; 2: {}", depth, d_min_npc.to_degrees() * 3600.0, d_max_eqr_top.to_degrees() * 3600.0);
    debug_assert!(d_min_eqr_top < d_max_eqr_top);
    let slope_eqr: f64 = (d_max_eqr_top - d_min_eqr_top) / (TRANSITION_LATITUDE - LAT_OF_SQUARE_CELL); // a = (yB - yA) / (xB - xA)
    let intercept_eqr: f64 = d_min_eqr_top - slope_eqr * LAT_OF_SQUARE_CELL;                           // b = yA - a * xA
    // EQR BOTTOM,  see comment of function largest_c2v_dist_in_eqr_bottom()
    let d_max = PI_OVER_FOUR * dist_cw;
    let coeff_cst_eqr: f64 = d_max;
    let coeff_x2_eqr: f64 = (d_min_eqr_top - d_max) / pow2(LAT_OF_SQUARE_CELL);
    // Struct creation
    ConstantsC2V {
      slope_npc,
      intercept_npc,
      slope_eqr,
      intercept_eqr,
      coeff_x2_eqr,
      coeff_cst_eqr
    }
  }
}

#[inline]
// #[allow(dead_code)]
pub fn haversine_dist(p1_lon: f64, p1_lat: f64, p2_lon: f64, p2_lat: f64) -> f64 {
  let shs = squared_half_segment(
    p2_lon - p1_lon, p2_lat - p1_lat, 
    p1_lat.cos(), p2_lat.cos());
  sphe_dist(shs)
}

/// Returns the angular distance corresponding to the given squared half great-circle arc segment
#[inline]
fn sphe_dist(squared_half_segment: f64) -> f64 {
  squared_half_segment.sqrt().asin().twice()
}

/// Returns `(s/2)^2` with `s` the segment (i.e. the Euclidean distance) between 
/// the two given points  `P1` and `P2` on the unit-sphere.
/// We recall that `s = 2 sin(ad/2)` with `ad` the angular distance between the two points.
/// # Input
/// - `dlon` the longitude difference, i.e. (P2.lon - P1.lon), in radians
/// - `dlat` the latitude difference, i.e. (P2.lat - P1.lat), in radians
/// - `cos_lat1` cosine of the latitude of the first point
/// - `cos_lat2` cosine of the latitude of the second point
#[inline]
fn squared_half_segment(dlon: f64, dlat: f64, cos_lat1: f64, cos_lat2: f64) -> f64 {
  dlat.half().sin().pow2() + cos_lat1 * cos_lat2 * dlon.half().sin().pow2()
}

#[inline]
fn to_squared_half_segment(spherical_distance: f64) -> f64 {
  spherical_distance.half().sin().pow2()
}

#[inline]
fn pow2(x: f64) -> f64 {
  x * x
}

/// Simple trait used to implements `pow2`, `twice` and `half` on f64.
pub trait Customf64 {
  fn pow2(self) -> f64;
  fn twice(self) -> f64;
  fn half(self) -> f64;
  fn div_eucl(self, rhs: f64) -> f64;
}

impl Customf64 for f64 {
  /// Returns x^2
  #[inline]
  fn pow2(self) -> f64 {
    self * self // or powi ?
  }
  /// Returns 2 * x
  #[inline]
  fn twice(self) -> f64 {
    2.0 * self // self + self (I hope the compiler know the simple shift bit to be used for x2)
  }
  /// Returns x / 2
  #[inline]
  fn half(self) -> f64 {
    0.5 * self
  }

  /// [Duplicated code](https://doc.rust-lang.org/std/primitive.f64.html#method.div_euc), because
  /// it is unstable so far.
  #[inline]
  fn div_eucl(self, rhs: f64) -> f64 {
    let q = (self / rhs).trunc();
    if self % rhs < 0.0 {
      return if rhs > 0.0 { q - 1.0 } else { q + 1.0 }
    }
    q
  }
}
// All types that implement `f64` get methods defined in `Customf64` for free.
// impl<F: f64> Customf64 for F {}*/


/// Returns an upper limit on the distance between a cell center around the given position 
/// and its furthest vertex.
/// # Params
/// - `depth` the depth of the cell
/// - `lon` the longitude of the point on the unit-sphere, in radians 
/// - `lat` the latitude of the point on the unit-sphere, in radians
/// 
/// # Result
/// The following plot shows, for the depth 8, the real largest distances (in red) and the result
/// of this method (in blue).  
/// WARNING: the units of `(lon, lat)` on the plot are *degrees*, while the distance is in *mas*  
/// Credit: plot made using [TOPCAT](http://www.star.bris.ac.uk/~mbt/topcat/)
/// ![CenterToVertexDist](https://raw.githubusercontent.com/cds-astro/cds-healpix-rust/master/resources/4doc/d_center_vertex.png)
/// 
pub fn largest_center_to_vertex_distance(depth: u8, lon: f64, lat: f64) -> f64 {
  // Specific case for depth 0
  if depth == 0 {
    return FRAC_PI_2 - TRANSITION_LATITUDE;
  }
  // Regular case
  let lat_abs = lat.abs();
  if lat_abs >= TRANSITION_LATITUDE {
    largest_c2v_dist_in_npc(lon, get_or_create(depth))
  } else if lat_abs >= LAT_OF_SQUARE_CELL {
    largest_c2v_dist_in_eqr_top(lat_abs, get_or_create(depth))
  } else {
    largest_c2v_dist_in_eqr_bottom(lat_abs, get_or_create(depth))
  }
}

/// Returns an upper limit on the distance between a cell center and it furthest vertex, for
/// all the cells in the region covered by a cone of given center and radius.  
/// It is an extension of [largest_center_to_vertex_distance](#fn.largest_center_to_vertex_distance)
/// # Params
/// - `depth` the depth of the cell
/// - `lon` the longitude of the point on the unit-sphere, in radians 
/// - `lat` the latitude of the point on the unit-sphere, in radians
/// - `radius` the radius of the cone, in radians
///
pub fn largest_center_to_vertex_distance_with_radius(depth: u8, lon: f64, lat: f64, radius: f64) -> f64 {
  // Specific case for depth 0
  if depth == 0 {
    return FRAC_PI_2 - TRANSITION_LATITUDE;
  }
  // Regular case
  let lat_abs = lat.abs();
  let lat_max = lat_abs + radius;
  let lat_min = lat_abs - radius;
  if lat_max >= TRANSITION_LATITUDE {
    largest_c2v_dist_in_npc_with_radius(lon, radius, get_or_create(depth))
  } else if lat_min >= LAT_OF_SQUARE_CELL {
    largest_c2v_dist_in_eqr_top_with_radius(lat_abs, radius, get_or_create(depth))
  } else if lat_max <= LAT_OF_SQUARE_CELL {
    largest_c2v_dist_in_eqr_bottom_with_radius(lat_abs, radius, get_or_create(depth))
  } else {
    let csts = get_or_create(depth);
    f64::max(
      largest_c2v_dist_in_eqr_top_with_radius(lat_abs, radius, csts),
      largest_c2v_dist_in_eqr_bottom_with_radius(lat_abs, radius, csts)
    )
  }
}

/// Same as [largest_center_to_vertex_distance_with_radius](#fn.largest_center_to_vertex_distance_with_radius)
/// but making the computation for several depths at the same time.
pub fn largest_center_to_vertex_distances_with_radius(mut from_depth: u8, to_depth: u8, lon: f64, lat: f64, radius: f64) -> Box<[f64]> {
  let mut vec: Vec<f64> = Vec::with_capacity((to_depth - from_depth) as usize);
  // Specific case for depth 0
  if from_depth == 0 {
    vec.push(FRAC_PI_2 - TRANSITION_LATITUDE);
    from_depth = 1_u8;
  }
  // Regular case
  let lat_abs = lat.abs();
  let lat_max = lat_abs + radius;
  let lat_min = lat_abs - radius;
  if lat_max >= TRANSITION_LATITUDE {
    let mut lon = (FRAC_PI_4 - (lon % FRAC_PI_2)).abs();
    lon = f64::min(lon + radius, FRAC_PI_4);
    for depth in from_depth..to_depth {
      let csts = get_or_create(depth);
      vec.push(linear_approx(lon, csts.slope_npc, csts.intercept_npc));
    }
  } else if lat_min >= LAT_OF_SQUARE_CELL {
    for depth in from_depth..to_depth {
      vec.push(largest_c2v_dist_in_eqr_top(lat_max, get_or_create(depth)));
    }
  } else if lat_max <= LAT_OF_SQUARE_CELL {
    let val_min = f64::max(lat_min, 0_f64);
    for depth in from_depth..to_depth {
      vec.push(largest_c2v_dist_in_eqr_bottom(val_min, get_or_create(depth))); 
    }
  } else {
    let val_max = f64::min(lat_max, TRANSITION_LATITUDE);
    let val_min = f64::max(lat_min, 0_f64);
    for depth in from_depth..to_depth {
      let csts = get_or_create(depth);
      vec.push(
        f64::max(
          largest_c2v_dist_in_eqr_top(val_max, csts),
          largest_c2v_dist_in_eqr_bottom(val_min, csts)
        )
      );
    }
  }
  vec.into_boxed_slice()
}

/// Returns an upper limit on distances between the center of a cell and its furthest vertex.  
/// We assumes that the cell center is located in the North polar cap region (OR IS ON THE
/// TRANSITION LATITUDE).  
/// We use a linear upper limit based on the longitude.
/// - At the transition latitude, we note CN the distance at a base cell border between:
///   - Cell center (C): (x_c = 1/nside, y_c = 1) => (lon_c = pi/4 * 1/nside, lat_c = TRANSITION_LATITUDE)
///   - North vertex (N): (x_n = 0, y_c = 2 - (1 + 1/nside)) => (lon_n = 0, lat_n = asin(1 - ((1 - 1/nside)^2 / 3)) 
///   - using the Haversine formula and SC defined below:
/// > dMax = CN = 2 * asin(sqrt( sin^2(SC/2)  + sin^2(pi/8) * cos(pi/4 * 1/nside)))
/// - At the transition latitude, we not SC the distance at a base cell center between:
///   - South vertex (S): (lon_s = pi/4, lat_s = TRANSITION_LATITUDE)
///   - Cell center: (lon_c = pi/4,  lat_c = asin(1 - ((1 - 1/nside)^2 / 3))
/// > dMin = SC =  asin(1 - ((1 - 1/nside)^2 / 3)) - TRANSITION_LATITUDE
/// - finally, linear approx:
/// > d = ((lon % pi/2) - pi/4) * (dMax - dMin)/(pi/4 * (1 - 1/nside)) + dMin
#[inline]
fn largest_c2v_dist_in_npc(lon: f64, csts: &ConstantsC2V) -> f64 {
  let lon = (FRAC_PI_4 - (lon % FRAC_PI_2)).abs();
  debug_assert!(0_f64 <= lon && lon <= FRAC_PI_4);
  linear_approx(lon, csts.slope_npc, csts.intercept_npc)
}
/// Same as the above method, but taking into account an additional radius
#[inline]
fn largest_c2v_dist_in_npc_with_radius(lon: f64, radius: f64, csts: &ConstantsC2V) -> f64 {
  debug_assert!(0_f64 < radius);
  let mut lon = (FRAC_PI_4 - (lon % FRAC_PI_2)).abs();
  debug_assert!(0_f64 <= lon && lon <= FRAC_PI_4);
  lon = f64::min(lon + radius, FRAC_PI_4);
  linear_approx(lon, csts.slope_npc, csts.intercept_npc)
}

/// Returns an upper limit on distance between the center of a cell and its furthest vertex.  
/// We assumes that the cell center is located in the equatorial region,
/// above the latitude at which cells are squares.
/// We use a linear upper limit based on the latitude.
/// - At the latitude in which cells are squares (lat = LAT_OF_SQUARE_CELL), we note CE the 
///   Center-to-East distance: CE (=CW) =  pi/4 * 1/nside * cos(LAT_OF_SQUARE_CELL)
/// > dMin = pi/4 * 1/nside * cos(LAT_OF_SQUARE_CELL)
/// - At the latitude under the transition latitude, we note CN the distance with N on the 
///   transition latitude and C at latitude: y_c = 1 - 1/nside => lat_c = arcsin(y_c * 2/3)
/// > dMax = TRANSITION_LATITUDE - arcsin(y_c * 2/3)
/// - finally, linear approx:
/// > d = (lat - LAT_OF_SQUARE_CELL) * (dMax - dMin)/(TRANSITION_LATITUDE - LAT_OF_SQUARE_CELL) + dMin
#[inline]
fn largest_c2v_dist_in_eqr_top(lat_abs: f64, csts: &ConstantsC2V) -> f64 {
  debug_assert!(LAT_OF_SQUARE_CELL <= lat_abs && lat_abs < TRANSITION_LATITUDE);
  linear_approx(lat_abs, csts.slope_eqr, csts.intercept_eqr)
}
/// Same as the above method, but taking into account an additional radius
#[inline]
fn largest_c2v_dist_in_eqr_top_with_radius(lat_abs: f64, radius: f64, csts: &ConstantsC2V) -> f64 {
  debug_assert!(0_f64 < radius);
  debug_assert!(LAT_OF_SQUARE_CELL <= lat_abs && lat_abs < TRANSITION_LATITUDE);
  largest_c2v_dist_in_eqr_top(f64::min(lat_abs + radius, TRANSITION_LATITUDE), csts)
}

/// Returns an upper limit on distance between the center of a cell and its furthest vertex.  
/// We assumes that the cell center is located in the equatorial region,
/// bellow the latitude at which cells are squares.
/// We use a parabola approximation upper limit based on the latitude.
/// - At lat = 0, d_max = pi/4 * 1/nside
/// - At lat = LAT_OF_SQUARE_CELL, d_min = pi/4 * 1/nside * cos(LAT_OF_SQUARE_CELL)
/// - Parabola approx: a * lat^2 + b = dist
///   - At lat = 0, b = d_max
///   - At lat = LAT_OF_SQUARE_CELL, a * LAT_OF_SQUARE_CELL^2 + d_max = d_min
///   - => a = (d_min - d_max) / LAT_OF_SQUARE_CELL^2
///   - => a = d_max (cos(LAT_OF_SQUARE_CELL) - 1) / LAT_OF_SQUARE_CELL^2
///   - => dist =  d_max * (1 + (cos(LAT_OF_SQUARE_CELL) - 1) / LAT_OF_SQUARE_CELL^2 * lat^2)
///   - => dist =  d_max * (1 - (1 - cos(LAT_OF_SQUARE_CELL)) / LAT_OF_SQUARE_CELL^2 * lat^2)
#[inline]
fn largest_c2v_dist_in_eqr_bottom(lat_abs: f64, csts: &ConstantsC2V) -> f64 {
  debug_assert!(0_f64 <= lat_abs && lat_abs <= LAT_OF_SQUARE_CELL);
  csts.coeff_x2_eqr * pow2(lat_abs) + csts.coeff_cst_eqr
}
/// Same as the above method, but taking into account an additional radius.
#[inline]
fn largest_c2v_dist_in_eqr_bottom_with_radius(lat_abs: f64, radius: f64, csts: &ConstantsC2V) -> f64 {
  debug_assert!(0_f64 < radius);
  debug_assert!(0_f64 <= lat_abs && lat_abs <= LAT_OF_SQUARE_CELL);
  largest_c2v_dist_in_eqr_bottom(f64::max(lat_abs - radius, 0_f64), csts)
}

#[inline]
fn linear_approx(x: f64, slope: f64, intercept: f64) -> f64 {
  slope * x + intercept
}

/// Returns, for the given depth, the number of cells along both axis of a base-resolution cell.
///
/// # Input
/// - `depth` must be in `[0, 29]`
///
/// # Output
/// - `nside` = 2^`depth`
///
/// # Panics
/// If `depth` is not valid (see [is_depth](fn.is_depth.html)), this method panics.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{nside};
///
/// assert_eq!(1, nside(0));
/// assert_eq!(2, nside(1));
/// assert_eq!(4, nside(2));
/// assert_eq!(8, nside(3));
/// assert_eq!(16, nside(4));
/// assert_eq!(32, nside(5));
/// assert_eq!(64, nside(6));
/// assert_eq!(128, nside(7));
/// assert_eq!(256, nside(8));
/// assert_eq!(512, nside(9));
/// assert_eq!(1024, nside(10));
/// assert_eq!(2048, nside(11));
/// assert_eq!(4096, nside(12));
/// assert_eq!(8192, nside(13));
/// assert_eq!(16384, nside(14));
/// assert_eq!(32768, nside(15));
/// assert_eq!(65536, nside(16));
/// assert_eq!(131072, nside(17));
/// assert_eq!(262144, nside(18));
/// assert_eq!(524288, nside(19));
/// assert_eq!(1048576, nside(20));
/// assert_eq!(2097152, nside(21));
/// assert_eq!(4194304, nside(22));
/// assert_eq!(8388608, nside(23));
/// assert_eq!(16777216, nside(24));
/// assert_eq!(33554432, nside(25));
/// assert_eq!(67108864, nside(26));
/// assert_eq!(134217728, nside(27));
/// assert_eq!(268435456, nside(28));
/// assert_eq!(536870912, nside(29));
/// // Using a for loop...
/// for depth in 0..29 {
///     assert_eq!(2u32.pow(depth), nside(depth as u8));
/// }
/// ```
#[inline]
pub fn nside(depth: u8) -> u32 {
    check_depth(depth);
    nside_unsafe(depth)
}

/// Same as [nside](fn.nside.html) except that this version does not check the argument, and thus
/// does not panics if the argument is illegal.
#[inline]
pub const fn nside_unsafe(depth: u8) -> u32 {
    1_u32 << depth
}


/// Returns, for the given difference of depth, the number of cells small cells the large cell
/// contains. If the small cell level is 0, the result is the sqaured nside.
///
/// # Input
/// - `delta_depth` must be in `[0, 29]`
///
/// # Output
/// - `nside^2` = 2^2*`delta_depth`
///
/// # Panics
/// If `delta_depth` is not in `[0, 29]`.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{nside_square};
/// 
/// for delta_depth in 0..29_u8 {
///     assert_eq!(2u64.pow(2 * delta_depth as u32), nside_square(delta_depth));
/// }
/// ```
#[inline]
pub fn nside_square(delta_depth: u8) -> u64 {
  check_depth(delta_depth);
  nside_square_unsafe(delta_depth)
}

/// Same as [nside_square](fn.nside_square.html) except that this version does not check the argument, 
/// and thus does not panics if the argument is illegal.
#[inline]
pub const fn nside_square_unsafe(delta_depth: u8) -> u64 {
  1_u64 << (delta_depth << 1)
}

#[inline]
fn check_depth(depth: u8) {
    assert!(is_depth(depth), "Expected depth in [0, 29]");
}

/// Returns `true` if the given argument is a valid depth, i.e. if it is <= [DEPTH_MAX](constant.DEPTH_MAX.html). 
#[inline]
pub const fn is_depth(depth: u8) -> bool {
    depth <= DEPTH_MAX
}

/// Returns, for the given `nside`, the number of subdivision of a base-resolution cell (i.e. the depth).
/// For the NESTED scheme only.
/// 
/// # Input
/// - `nside` must be a power of 2 in `[0, 2^29]`
///
/// # Output
/// - `depth` = `log2(nside)`
///
/// # Panics
/// If `nside` is not valid (see [is_nside](fn.is_nside.html)), this method panics.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{nside, depth};
///
/// for d in 0..29 {
///     assert_eq!(d, depth(nside(d as u8)));
/// }
/// ```
#[inline]
pub fn depth(nside: u32) -> u8 {
    check_nside(nside);
    depth_unsafe(nside)
}

/// Same as [depth](fn.depth.html) except that this version does not check the argument, and thus
/// does not panics if the argument is illegal.
#[inline]
pub const fn depth_unsafe(nside: u32) -> u8 {
    nside.trailing_zeros() as u8
}

#[inline]
fn check_nside(nside: u32) {
    assert!(is_nside(nside), "Nside must be a power of 2 in [1-2^29]");
}

/// Returns `true` if the given argument is a valid `nside` for the NESTED scheme, i.e. 
/// if it is a power of 2, is != 0 and is <= [NSIDE_MAX](constant.NSIDE_MAX.html). 
#[inline]
pub fn is_nside(nside: u32) -> bool {
    is_pow_of_2(nside) && nside > 0 && nside <= NSIDE_MAX
}

/// Determines if an integer is a power of two, including 0.
/// Taken from the "Bit Twiddling Hacks" web page of Sean Eron Anderson.
#[inline]
fn is_pow_of_2(x: u32) -> bool {
    (x & (x - 1)) == 0
}

/// Returns the number of distinct hash value (the number of cells or pixel the unit sphere is
/// devided in) at the given `depth`.
/// 
/// # Input
/// - `depth` must be in `[0, 29]`
///
/// # Output
/// - `n_hash` = `12 * nside^2`
///
/// # Panics
/// If `depth` is not valid (see [is_depth](fn.is_depth.html)), this method panics.
///
/// # Examples
/// 
/// ```rust
/// use cdshealpix::{n_hash};
/// 
/// assert_eq!(12u64, n_hash(0u8));
/// assert_eq!(48u64, n_hash(1u8));
/// assert_eq!(192u64, n_hash(2u8));
/// assert_eq!(768u64, n_hash(3u8));
/// assert_eq!(3072u64, n_hash(4u8));
/// assert_eq!(12288u64, n_hash(5u8));
/// assert_eq!(49152u64, n_hash(6u8));
/// assert_eq!(196608u64, n_hash(7u8));
/// assert_eq!(786432u64, n_hash(8u8));
/// assert_eq!(3145728u64, n_hash(9u8));
/// assert_eq!(12582912u64, n_hash(10u8));
/// assert_eq!(50331648u64, n_hash(11u8));
/// assert_eq!(201326592u64, n_hash(12u8));
/// assert_eq!(805306368u64, n_hash(13u8));
/// assert_eq!(3221225472u64, n_hash(14u8));
/// assert_eq!(12884901888u64, n_hash(15u8));
/// assert_eq!(51539607552u64, n_hash(16u8));
/// assert_eq!(206158430208u64, n_hash(17u8));
/// assert_eq!(824633720832u64, n_hash(18u8));
/// assert_eq!(3298534883328u64, n_hash(19u8));
/// assert_eq!(13194139533312u64, n_hash(20u8));
/// assert_eq!(52776558133248u64, n_hash(21u8));
/// assert_eq!(211106232532992u64, n_hash(22u8));
/// assert_eq!(844424930131968u64, n_hash(23u8));
/// assert_eq!(3377699720527872u64, n_hash(24u8));
/// assert_eq!(13510798882111488u64, n_hash(25u8));
/// assert_eq!(54043195528445952u64, n_hash(26u8));
/// assert_eq!(216172782113783808u64, n_hash(27u8));
/// assert_eq!(864691128455135232u64, n_hash(28u8));
/// assert_eq!(3458764513820540928u64, n_hash(29u8));
/// ```
/// 
#[inline]
pub fn n_hash(depth: u8) -> u64 {
  check_depth(depth);
  n_hash_unsafe(depth)
}

/// Same as [n_hash](fn.n_hash.html) except that this version does not panic if the given `depth` is
/// out of range.
#[inline]
pub const fn n_hash_unsafe(depth: u8) -> u64 { 12u64 << (depth << 1u8) }

/// Returns `true` if the function [best_starting_depth](fn.best_starting_depth.html) is valid
/// for the given argument `d_max_rad`. So if `d_max_rad < ~48 deg`. `d_max_rad` is given in radians.
/// 
/// ```rust
/// use cdshealpix::{has_best_starting_depth};
/// use std::f64::consts::PI;
/// 
/// assert!(!has_best_starting_depth(PI / 3f64));
/// assert!(has_best_starting_depth(PI / 4f64));
/// ```
#[inline]
pub fn has_best_starting_depth(d_max_rad: f64) -> bool {
    d_max_rad < SMALLER_EDGE2OPEDGE_DIST[0]
}

/// Returns the the smallest depth (in `[0, 29]`) at which a shape having the given largest distance
/// from its center to a border overlaps a maximum of 9 cells (the cell containing the center of
/// the shape plus the 8 neighbouring cells).  
/// Info: internally, unrolled binary search loop on 30 pre-computed values (one by depth).

/// @return -1 if the given distance is very large (> ~48deg), else returns the smallest depth
/// (in [0, 29]) at which a shape having the given largest distance from its center to a border
/// overlaps a maximum of 9 cells (the cell containing the center of the shape plus the 8
/// neighbouring cells).
/// 
/// # Input
/// - `d_max_rad` largest possible distance, in radians, between the center and the border of a shape
///
/// # Output
/// - `depth` = the smallest depth (in `[0, 29]`) at which a shape having the given largest distance 
/// from its center to a border overlaps a maximum of 9 cells (the cell containing the center of the
/// shape plus the 8 neighbouring cells).
///
/// # Panics
/// If the given distance is very large (> ~48deg), this function is not valid since the 12 base
/// cells could be overlaped by the shape 
/// (see [has_best_starting_depth](fn.has_best_starting_depth.html)). Thus it panics.
///
/// # Examples
///
/// ```rust
/// use cdshealpix::{best_starting_depth};
/// use std::f64::consts::PI;
///
/// assert_eq!(0, best_starting_depth(PI / 4f64)); // 45 deg
/// assert_eq!(5, best_starting_depth(0.0174533)); //  1 deg
/// assert_eq!(7, best_starting_depth(0.0043632)); // 15 arcmin
/// assert_eq!(9, best_starting_depth(0.0013));    // 4.469 arcmin
/// assert_eq!(15, best_starting_depth(1.454E-5)); // 3 arcsec
/// assert_eq!(20, best_starting_depth(6.5E-7));   // 0.134 arcsec
/// assert_eq!(22, best_starting_depth(9.537E-8)); // 20 mas
/// ```
#[inline]
pub fn best_starting_depth(d_max_rad: f64) -> u8 { // Could have used an Option
    assert!(d_max_rad < SMALLER_EDGE2OPEDGE_DIST[0],
            "Too large value, use first function has_best_starting_depth");
    // Unrolled binary search loop
    if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[29] {
        29
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[15] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[22] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[25] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[27] {
                    if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[28] {
                        28
                    } else {
                        27
                    } 
                } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[26] {
                    26
                } else {
                    25
                } 
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[24] {
                24
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[23] {
                23
            } else {
                22
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[18] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[20] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[21] {
                    21
                } else {
                    20
                }
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[19] {
                19
            } else {
                18
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[17] {
            17
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[16] {
            16
        } else {
            15
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[7] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[11] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[13] {
                if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[14] {
                    14
                } else {
                    13
                }
            } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[12] {
                12
            } else {
                11
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[9] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[10] {
                10
            } else {
                9
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[8] {
            8
        } else {
            7
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[3] {
        if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[5] {
            if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[6] {
                6
            } else {
                5
            }
        } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[4] {
            4
        } else {
            3
        }
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[2] {
        2
    } else if d_max_rad < SMALLER_EDGE2OPEDGE_DIST[1] {
        1
    } else {
        0
    }
}

/// Performs the HEALPix projection: `(x, y) = proj(lon, lat)`.  
/// The chosen scale is such that: base cell vertices and center coordinates are integers;
/// the distance from a cell center to its vertices equals one.  
/// This projection is multi-purpose in the sense that if `lon` is in `[-pi, pi]`, then
/// `x` is in `[-4, 4]` and if `lon` is in `[0, 2pi]`, then `x` is in `[0, 8]`.  
/// It means that a same position on the sphere can lead to different positions in the projected 
/// Euclidean plane.
/// 
/// Simplified projection formulae are:
///  - Equatorial region
/// ```math
/// \boxed{
///   \left\{
///     \begin{array}{lcl}
///       X & = & \alpha \times \frac{4}{\pi} \\
///       Y & = & \sin(\delta) \times \frac{3}{2}
///     \end{array}
///   \right.
/// }
/// \Rightarrow
/// \left\{
///   \begin{array}{lcl}
///     \alpha \in [0, 2\pi] & \leadsto &  X \in [0, 8] \\
///     \sin\delta \in [-\frac{2}{3}, \frac{2}{3}] & \leadsto & Y \in [-1, 1]
///   \end{array}
/// \right.
/// ```
///  - Polar caps:
/// ```math
/// \boxed{
///   \left\{
///     \begin{array}{lcl}
///       t & = & \sqrt{3(1-\sin\delta)} \\
///       X & = & (\alpha\frac{4}{\pi} - 1)t+1 \\
///       Y & = & 2 - t
///     \end{array}
///   \right.
/// }
/// \Rightarrow
/// \left\{
///   \begin{array}{l}
///     \alpha \in [0, \frac{\pi}{2}] \\
///     \sin\delta \in ]\frac{2}{3}, 1]
///   \end{array}
/// \right.
/// \leadsto
/// \begin{array}{l}
///    t \in [0, 1[  \\
///    X \in ]0, 2[ \\
///    Y \in ]1, 2]
/// \end{array}
/// ```
/// 
/// It is the responsibility of the caller to homogenize the result according to its needs.
/// ![Proj](https://raw.githubusercontent.com/cds-astro/cds-healpix-rust/master/resources/4doc/hpx_proj.png)
///
/// # Inputs
/// - `lon` longitude in radians, support positive and negative reasonably large values
///   (naive approach, no Cody-Waite nor Payne Hanek range reduction).
/// - `lat` latitude in radians, must be in `[-pi/2, pi/2]`
///
/// # Output
/// - `(x, y)` the projected planar Euclidean coordinates of the point
///    of given coordinates `(lon, lat)` on the unit sphere
///     - `lon` &le; `0` => `x in [-8, 0]`
///     - `lon` &ge; `0` => `x in [0, 8]`
///     - `y in [-2, 2]`
///
/// # Panics
/// If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
///
/// # Examples
/// To obtain the WCS projection (see Calabretta2007), you can write:
/// ```rust
/// use cdshealpix::proj;
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
///
/// let lon = 25.1f64;
/// let lat = 46.7f64;
///
/// let (mut x, mut y) = proj(lon.to_radians(), lat.to_radians());
/// if x > 4f64 {
///   x -= 8f64;
/// }
/// x *= FRAC_PI_4;
/// y *= FRAC_PI_4;
///
/// assert!(-PI <= x && x <= PI);
/// assert!(-FRAC_PI_2 <= y && y <= FRAC_PI_2);
/// ```
///
/// Other test example:
/// ```rust
/// use cdshealpix::{TRANSITION_LATITUDE, proj};
/// use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
/// 
/// let (x, y) = proj(0.0, 0.0);
/// assert_eq!(0f64, x);
/// assert_eq!(0f64, y);
///
/// assert_eq!((0.0, 1.0), proj(0.0 * FRAC_PI_2, TRANSITION_LATITUDE));
///
/// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
///     f64::sqrt((p2.0 - p1.0) * (p2.0 - p1.0) + (p2.1 - p1.1) * (p2.1 - p1.1))
/// }
/// assert!(dist((0.0, 0.0), proj(0.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((0.0, 1.0), proj(0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((1.0, 2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((2.0, 1.0), proj(1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((3.0, 2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((4.0, 1.0), proj(2.0 * FRAC_PI_2 , TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((5.0, 2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((6.0, 1.0), proj(3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((7.0, 2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((0.0, 1.0), proj(4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((0.0, 0.0), proj(4.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((0.0, -1.0), proj(0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((1.0, -2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((2.0, -1.0), proj(1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((3.0, -2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((4.0, -1.0), proj(2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((5.0, -2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((6.0, -1.0), proj(3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((7.0, -2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((0.0, -1.0), proj(4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
///
/// assert!(dist((-0.0, 0.0), proj(-0.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((-0.0, 1.0), proj(-0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-1.0, 2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-2.0, 1.0), proj(-1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-3.0, 2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-4.0, 1.0), proj(-2.0 * FRAC_PI_2 , TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-5.0, 2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-6.0, 1.0), proj(-3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-7.0, 2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
/// assert!(dist((-0.0, 1.0), proj(-4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-0.0, 0.0), proj(-4.0 * FRAC_PI_2, 0.0)) < 1e-15);
/// assert!(dist((-0.0, -1.0), proj(-0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-1.0, -2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-2.0, -1.0), proj(-1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-3.0, -2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-4.0, -1.0), proj(-2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-5.0, -2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-6.0, -1.0), proj(-3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// assert!(dist((-7.0, -2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
/// assert!(dist((-0.0, -1.0), proj(-4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
/// ```
#[inline]
pub fn proj(lon: f64, lat: f64) -> (f64, f64) {
    check_lat(lat);
    let lon = abs_sign_decompose(lon);
    let lat = abs_sign_decompose(lat);
    let x = pm1_offset_decompose(lon.abs * FOUR_OVER_PI);
    let mut xy = (x.pm1, lat.abs);
    if is_in_equatorial_region(lat.abs) {
        proj_cea(&mut xy);
    } else {
        proj_collignon(&mut xy);
    }
    apply_offset_and_signs(&mut xy, x.offset, lon.sign, lat.sign);
    xy
}

/// Returns the hash of the base cell given the coordinates of a points in the Euclidean projection
/// plane.
/// The purpose so far is just to test and compare both speed and precision with the 45deg 
/// rotation solution.
/// 
/// # Input
/// - `(x, y)` the coordinates in the Euclidean projection plane, i.e. $x \in [0, 8[$ 
/// and $y \in [-2, 2]$
/// 
/// # Ouput
/// - `d0h` the hash value of the base cell (i.e. the depth 0 / nside 1 cell)
/// 
/// # Example
/// Simple example based on the center of each base cell.
/// ```rust
/// use cdshealpix::base_cell_from_proj_coo;
/// 
/// assert_eq!(base_cell_from_proj_coo(1.0,  1.0),  0);
/// assert_eq!(base_cell_from_proj_coo(3.0,  1.0),  1);
/// assert_eq!(base_cell_from_proj_coo(5.0,  1.0),  2);
/// assert_eq!(base_cell_from_proj_coo(7.0,  1.0),  3);
/// assert_eq!(base_cell_from_proj_coo(0.0,  0.0),  4);
/// assert_eq!(base_cell_from_proj_coo(2.0,  0.0),  5);
/// assert_eq!(base_cell_from_proj_coo(4.0,  0.0),  6);
/// assert_eq!(base_cell_from_proj_coo(6.0,  0.0),  7);
/// assert_eq!(base_cell_from_proj_coo(1.0, -1.0),  8);
/// assert_eq!(base_cell_from_proj_coo(3.0, -1.0),  9);
/// assert_eq!(base_cell_from_proj_coo(5.0, -1.0), 10);
/// assert_eq!(base_cell_from_proj_coo(7.0, -1.0), 11);
/// ```
pub fn base_cell_from_proj_coo(x: f64, y: f64) -> u8 {
  let mut x = 0.5 * ensures_x_is_positive(x);
  let mut y = 0.5 * (y + 3.0);
  let mut i = x as u8;     debug_assert!(i <  4); // if can be == 4, then (x as u8) & 3
  let mut j = (y as u8) << 1;  debug_assert!(j == 0 || j == 2 || j == 4);
  x -= i as f64;               debug_assert!(0.0 <= x && x < 1.0);
  y -= (j >> 1) as f64;        debug_assert!(0.0 <= y && y < 1.0);
  let in_northwest = (x <= y) as u8;       // 1/0
  let in_southeast = (x >= 1.0 - y) as u8; // 0\1
  i += in_southeast >> in_northwest; // <=> in_southeast & (1 - in_northwest) => 0 or 1
  j += in_northwest + in_southeast;
  if j == 6 { j = 4; } // Very rare case (North pole), so few risks of branch miss-prediction
  debug_assert!(j == 2 || j == 3 || j == 4);
  ((4 - j) << 2) + i
}

/// Unproject the given HEALPix projected points.  
/// This unprojection is multi-purpose in the sense that:
///  - if input `x` in `[-8, 0[`, then output `lon` in `[-2pi, 0]`
///  - if input `x` in `[ 0, 8]`, then output `lon` in `[0, 2pi]`
///  - output `lat` always in `[-pi/2, pi/2]`
/// 
/// # Inputs
///  - `x` the projected coordinate along the x-axis, supports positive and negative reasonably 
///        large values with a naive approach (no Cody-Waite nor Payne Hanek range reduction).
///  - `y` the projected coordinate along te x-axis, must be in `[-2, 2]`
/// 
/// # Output
/// -  `(lon, lat)` in radians, the position on the unit sphere whose projected coordinates are 
///    the input coordinates `(x, y)`.  
///   - if `x <= 0`, then `lon` in `[-2pi, 0]`;
///   - else if `x >= 0`, the  `lon` in `[0, 2pi]`
///   - `lat` always in `[-pi/2, pi/2]`.
///
/// # Panics
/// If `y` **not in** `[-2, 2]`, this method panics.
///
/// # Examples
/// To obtain the WCS un-projection (see Calabretta2007), you can write:
/// ```rust
/// use cdshealpix::{FOUR_OVER_PI, unproj};
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
///
/// let x = 2.1f64;
/// let y = 0.36f64;
///
/// let (mut lon, mut lat) = unproj(x * FOUR_OVER_PI, y * FOUR_OVER_PI);
/// if lon < 0f64 {
///     lon += 2f64 * PI;
/// }
///
/// assert!(0f64 <= lon && lon <= 2f64 * PI);
/// assert!(-FRAC_PI_2 <= lat && lat <= FRAC_PI_2);
/// ```
///
/// Other test example:
/// ```rust
/// use cdshealpix::{TRANSITION_LATITUDE, proj, unproj};
/// use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
///
/// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
///     let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
///     let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
///     2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
/// }
/// 
/// let points: [(f64, f64); 40] = [
///     (0.0 * FRAC_PI_2, 0.0),
///     (1.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (2.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (3.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (4.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (4.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (0.0 * FRAC_PI_2, 0.0),
///     (1.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (2.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (3.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (4.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (4.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-0.0 * FRAC_PI_2, 0.0),
///     (-1.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-2.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-3.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-4.0 * FRAC_PI_2, TRANSITION_LATITUDE),
///     (-0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-4.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2),
///     (-0.0 * FRAC_PI_2, 0.0),
///     (-1.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-2.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-3.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-4.0 * FRAC_PI_2, -TRANSITION_LATITUDE),
///     (-0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2),
///     (-4.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)
/// ];
/// 
/// for (lon, lat) in points.iter() {
///     let (x, y): (f64, f64) = proj(*lon, *lat);
///     assert!(dist((*lon, *lat), unproj(x, y)) < 1e-15);
/// }
/// ```
#[inline]
pub fn unproj(x: f64, y: f64) -> (f64, f64) {
    check_y(y);
    let x = abs_sign_decompose(x);
    let y = abs_sign_decompose(y);
    let lon = pm1_offset_decompose(x.abs);
    let mut lonlat= (lon.pm1, y.abs);
    if is_in_projected_equatorial_region(y.abs) {
        deproj_cea(&mut lonlat);
    } else {
        deproj_collignon(&mut lonlat);
    }
    apply_offset_and_signs(&mut lonlat, lon.offset, x.sign, y.sign);
    lonlat.0 *= FRAC_PI_4;
    lonlat
}

/// In the case we want the projection to return values `x in [0, 8]`, we just have to apply
/// this method to the returned `x` value.
#[inline]
pub(crate) fn ensures_x_is_positive(x: f64) -> f64 {
  if x < 0.0 { x + 8.0 } else { x }
}

/// Verify that the latitude is in [-PI/2, PI/2], panics if not.
#[inline]
fn check_lat(lat: f64) {
    assert!(-FRAC_PI_2 <= lat && lat <= FRAC_PI_2);
}

/// Verify that the latitude is in [-PI/2, PI/2], panics if not.
#[inline]
fn check_lat_res(lat: f64) -> Result<(), String> {
  if -FRAC_PI_2 <= lat && lat <= FRAC_PI_2 {
    Ok(())
  } else {
    Err(format!("Wrong latitude. Expected value in [-pi/2, pi/2]. Actual: {}", lat))
  }
}

/// Verify that the projected y coordinate is in [-2, 2], panics if not.
#[inline]
fn check_y(y: f64) { assert!(-2f64 <= y && y <= 2f64); }

/// Returns `true` if the point of given (absolute value of) latitude is in the equatorial region,
/// and `false` if it is located in one of the two polar caps
#[inline]
pub fn is_in_equatorial_region(abs_lat: f64) -> bool {
    abs_lat <= TRANSITION_LATITUDE
}

/// Returns `true` if the point of given (absolute value of) y coordinate in the projected plane
/// is in the equatorial region, and `false` if it is located in one of the two polar caps
#[inline]
pub fn is_in_projected_equatorial_region(abs_y: f64) -> bool { abs_y <= 1.0 }

// Returns the absolute value of the given double together with its bit of sign
struct AbsAndSign {
    abs: f64,
    sign: u64,
}
#[inline]
pub(crate) fn abs_sign_decompose(x: f64) -> AbsAndSign {
    let bits = f64::to_bits(x);
    AbsAndSign {
        abs: f64::from_bits(bits & F64_BUT_SIGN_BIT_MASK),
        sign: bits & F64_SIGN_BIT_MASK,
    }
}

// Decompose the given positive real value in
// --* an integer offset in [1, 3, 5, 7] (*PI/4) and
// --* a real value in [-1.0, 1.0] (*PI/4)
pub(crate) struct OffsetAndPM1 {
    offset: u8, // = 1, 3, 5 or 7
    pm1: f64,   // in [-1.0, 1.0]
}
#[inline]
pub(crate) fn pm1_offset_decompose(x: f64) -> OffsetAndPM1 {
    let floor: u8 = x as u8;
    let odd_floor: u8 = floor | 1u8;
    OffsetAndPM1 {
        offset: odd_floor & 7u8, // value modulo 8
        pm1: x - (odd_floor as f64),
    }
}

// Cylindrical Equal Area projection
#[inline]
pub(crate) fn proj_cea(xy: &mut (f64, f64)) {
    let (_, ref mut y) = *xy;
    *y = f64::sin(*y) * ONE_OVER_TRANSITION_Z;
}
#[inline]
fn deproj_cea(lonlat: &mut (f64, f64)) {
    let (_, ref mut lat) = *lonlat;
    // Using asin is OK here since |lat*TRANSITION_Z| < 2/3, so not near from 1.
    *lat = f64::asin((*lat) * TRANSITION_Z);
}

// Collignon projection
#[inline]
pub(crate) fn proj_collignon(xy: &mut (f64, f64)) {
    let (ref mut x, ref mut y) = *xy;
    *y = SQRT6 * f64::cos(HALF * *y + FRAC_PI_4);
    *x *= *y;
    *y = 2.0 - *y;
}
#[inline]
fn deproj_collignon(lonlat: &mut (f64, f64)) {
    let (ref mut lon, ref mut lat) = *lonlat;
    *lat = 2.0 - *lat;
    if is_not_near_from_pole(*lat) { // Rare, so few risks of branch miss-prediction
        *lon /= *lat;
        deal_with_numerical_approx_in_edges(lon);
    } // in case of pole, lon = lat = 0 (we avoid NaN due to division by lat=0)
    *lat *= ONE_OVER_SQRT6;
    // Using acos is OK here since lat < 1/sqrt(6), so not near from 1.
    *lat = 2.0 * f64::acos(*lat) - FRAC_PI_2;
}

#[inline]
fn is_not_near_from_pole(sqrt_of_three_time_one_minus_sin_of: f64) -> bool {
    // In case of pole: x = y = 0
    sqrt_of_three_time_one_minus_sin_of > EPS_POLE
}

#[inline]
fn deal_with_numerical_approx_in_edges(lon: &mut f64) {
    if *lon > 1.0 {
        *lon = 1.0;
    } else if *lon < -1.0 {
        *lon = -1.0;
    }
}

// Shift x by the given offset and apply lon and lat signs to x and y respectively
#[inline]
pub(crate) fn apply_offset_and_signs(ab: &mut (f64, f64), off: u8, a_sign: u64, b_sign: u64) {
    let (ref mut a, ref mut b) = *ab;
    *a += off as f64;
    *a = f64::from_bits(f64::to_bits(*a) | a_sign);
    *b = f64::from_bits(f64::to_bits(*b) | b_sign);
}

// Import module compass point
pub mod compass_point;
pub mod external_edge;
use crate::compass_point::{MainWind};
use crate::compass_point::MainWind::*;

/// Compute the base cell value which is the neighbour of the given base cell, in the given direction.  
/// There is no neighbour:
/// - in the North and South directions for the equatorial region cells (i.e. cells 4, 5, 6 and 7)
/// - in the East and West directions for:
///   - the north polar cap cells (i.e. cells 0, 1, 2 and 3)
///   - the south polar cap cells (i.e. cells 8, 9, 10 and 11)  
pub fn neighbour(base_cell: u8, direction: MainWind) -> Option<u8> {
  if direction == MainWind::C {
    Some(base_cell)
  } else {
    let d0h_mod_4 = base_cell & 3_u8;  // <=> base_cell modulo 4
    match base_cell >> 2 { // <=> basce_cell / 4
      0 => npc_neighbour(d0h_mod_4, direction),
      1 => eqr_neighbour(d0h_mod_4, direction),
      2 => spc_neighbour(d0h_mod_4, direction),
      _ => panic!("Base cell must be in [0, 12["),
    }
  }
}

fn npc_neighbour(d0h_mod_4: u8, direction: MainWind) -> Option<u8> {
  match direction {
     S => base_cell_opt(iden(d0h_mod_4), 2),
    SE => base_cell_opt(next(d0h_mod_4), 1),
    SW => base_cell_opt(iden(d0h_mod_4), 1),
    NE => base_cell_opt(next(d0h_mod_4), 0),
    NW => base_cell_opt(prev(d0h_mod_4), 0),
     N => base_cell_opt(oppo(d0h_mod_4), 0),
    _ => None,
  }  
}

fn eqr_neighbour(d0h_mod_4: u8, direction: MainWind) -> Option<u8> {
  match direction {
    SE => base_cell_opt(iden(d0h_mod_4), 2),
     E => base_cell_opt(next(d0h_mod_4), 1),
    SW => base_cell_opt(prev(d0h_mod_4), 2),
    NE => base_cell_opt(iden(d0h_mod_4), 0),
     W => base_cell_opt(prev(d0h_mod_4), 1),
    NW => base_cell_opt(prev(d0h_mod_4), 0),
    _ => None,
  }
}

fn spc_neighbour(d0h_mod_4: u8, direction: MainWind) -> Option<u8> {
  match direction {
     S => base_cell_opt(oppo(d0h_mod_4), 2),
    SE => base_cell_opt(next(d0h_mod_4), 2),
    SW => base_cell_opt(prev(d0h_mod_4), 2),
    NE => base_cell_opt(next(d0h_mod_4), 1),
    NW => base_cell_opt(iden(d0h_mod_4), 1),
     N => base_cell_opt(iden(d0h_mod_4), 0),
    _ => None,
  }
}

/// Returns the direction of a cell on the inner edge of the given base cell from its neighbour 
/// located at the given direction in a different base cell.
/// # Inputs
/// - `base_cell` the base cell containing the sub-cell we are looking for the direction from its
///   neighbour in the given `neighbour_direction`
/// - `inner_direction` the direction of the sub-cell in the edge of the given base cell
/// - `neighbour_direction` direction of the neighbour of the sub-cell from which we are looking 
///    at the direction of the sub-cell
///   
pub fn edge_cell_direction_from_neighbour(base_cell: u8, inner_direction: &MainWind, neighbour_direction: &MainWind) -> MainWind {
  match base_cell >> 2 { // <=> basce_cell / 4
    0 => npc_egde_direction_from_neighbour(inner_direction, neighbour_direction),
    1 => eqr_edge_direction_from_neighbour(inner_direction, neighbour_direction),
    2 => spc_edge_direction_from_neighbour(inner_direction, neighbour_direction),
    _ => panic!("Base cell must be in [0, 12["),
  }
}

fn npc_egde_direction_from_neighbour(inner_direction: &MainWind, neighbour_direction: &MainWind) -> MainWind {
  match neighbour_direction {
    C => panic!("No neighbour in direction {:?}", &neighbour_direction),
    E => match inner_direction {
      N | NE => N,
      E => panic!("No neighbour in direction {:?}", &neighbour_direction),
      S | SE => neighbour_direction.opposite(),
      _ => unreachable!(),
    },
    W => match inner_direction {
      N | NW => N,
      W => panic!("No neighbour in direction {:?}", &neighbour_direction),
      S | SW => neighbour_direction.opposite(),
      _ => unreachable!(),
    },
    NE => {
      assert!(*inner_direction == N || *inner_direction == E || *inner_direction == NE);
      NW
    },
    NW => {
      assert!(*inner_direction == N || *inner_direction == W || *inner_direction == NW);
      NE
    },
    N  => match inner_direction {
      N => N,
      E | NE => W,
      W | NW => E,
      _ => unreachable!(),
    },
    _ => neighbour_direction.opposite(),
  }
}

fn eqr_edge_direction_from_neighbour(_inner_direction: &MainWind, neighbour_direction: &MainWind) -> MainWind {
  neighbour_direction.opposite()
}

fn spc_edge_direction_from_neighbour(inner_direction: &MainWind, neighbour_direction: &MainWind) -> MainWind {
  match neighbour_direction {
    C => panic!("No neighbour in direction {:?}", &neighbour_direction),
    E => match inner_direction {
      S | SE => S,
      E => panic!("No neighbour in direction {:?}", &neighbour_direction),
      N | NE => neighbour_direction.opposite(),
      _ => unreachable!(),
    },
    W => match inner_direction {
      S | SW => S,
      W => panic!("No neighbour in direction {:?}", &neighbour_direction),
      N | NW => neighbour_direction.opposite(),
      _ => unreachable!(),
    },
    SE => {
      assert!(*inner_direction == S || *inner_direction == E || *inner_direction == SE);
      SW
    },
    SW => {
      assert!(*inner_direction == S || *inner_direction == W || *inner_direction == SW);
      SE
    },
    S  => match inner_direction {
      S => S,
      E | SE => W,
      W | SW => E,
      _ => unreachable!(),
    },
    _ => neighbour_direction.opposite(),
  }
}

/// Returns the direction of the given base cell from its neighbour base cell located 
/// in the given direction.
/// # Panics
/// If the base cell has no neighbour in the given direction (i.e. N/S for equatorial cells
/// and E/W for polar caps cells)
pub fn direction_from_neighbour(base_cell: u8, neighbour_direction: &MainWind) -> MainWind {
  match base_cell >> 2 { // <=> basce_cell / 4
    0 => npc_direction_from_neighbour(neighbour_direction),
    1 => eqr_direction_from_neighbour(neighbour_direction),
    2 => spc_direction_from_neighbour(neighbour_direction),
    _ => panic!("Base cell must be in [0, 12["),
  }
}

fn npc_direction_from_neighbour(neighbour_direction: &MainWind) -> MainWind {
  match neighbour_direction {
    E | W | C => panic!("No neighbour in direction {:?}", &neighbour_direction),
    NE => NW,
    NW => NE,
    N  => N,
    _ => neighbour_direction.opposite(),
  }
}

fn eqr_direction_from_neighbour(neighbour_direction: &MainWind) -> MainWind {
  match neighbour_direction {
    S | N | C => panic!("No neighbour in direction {:?}", &neighbour_direction),
    _ => neighbour_direction.opposite(),
  }
}

fn spc_direction_from_neighbour(neighbour_direction: &MainWind) -> MainWind {
  match neighbour_direction {
    E | W | C => panic!("No neighbour in direction {:?}", &neighbour_direction),
    S  => S,
    SE => SW,
    SW => SE,
    _ => neighbour_direction.opposite(),
  }
}

/// Returns (mod4 - 1) in [0, 2], and 3 if mod4 == 0 (i.e. the previous value in [0, 3] range)
#[inline]
fn prev(mod4: u8) -> u8 {
  debug_assert!(mod4 < 4);
  (((mod4 as i8) - 1) & 3) as u8
}

/// Returns (mod4 + 1) in [1, 3], and 0 if mod4 == 3 (i.e. the next value in [0, 3] range)
#[inline]
fn next(mod4: u8) -> u8 {
  debug_assert!(mod4 < 4);
  (mod4 + 1) & 3
}

/// Returns (mod4 + 2) in [2, 3], and 0 if mod4 == 2 and 1 if mod4 == 3 (i.e. the opposite value in [0, 3] range)
#[inline]
fn oppo(mod4: u8) -> u8 {
  debug_assert!(mod4 < 4);
  (mod4 + 2) & 3
}

/// Returns the input value: useless, just used to improve code legibility
#[inline]
fn iden(mod4: u8) -> u8 {
  debug_assert!(mod4 < 4);
  mod4
}

#[inline]
fn base_cell_opt(i: u8, j: u8) -> Option<u8> {
  Some(base_cell(i, j))
}

/// Compute the base cell from its (i, j) coordinates:
/// - i: index along the longitude axis ( = base_cell modulo 4)
/// - j: index along the latitude axis ( = base_cell / 4)
///   - = 0 for the cells covering the north polar cap
///   - = 1 for the cells with are only in the equatorial region
///   - = 2 for the cells covering the south polar cap
#[inline]
fn base_cell(i: u8, j: u8) -> u8 {
  debug_assert!(i < 4 && j < 3);
  (j << 2) + i
}


/// Module containing NESTED scheme methods
pub mod nested;

/// Module containing RING scheme methods
pub mod ring;

/// No need to make those public!
mod xy_geom;
pub mod sph_geom;
mod special_points_finder;

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn testok_nside() {
    assert_eq!(1, nside(0));
    assert_eq!(2, nside(1));
    assert_eq!(4, nside(2));
    assert_eq!(8, nside(3));
    assert_eq!(16, nside(4));
    assert_eq!(32, nside(5));
    assert_eq!(64, nside(6));
    assert_eq!(128, nside(7));
    assert_eq!(256, nside(8));
    assert_eq!(512, nside(9));
    assert_eq!(1024, nside(10));
    assert_eq!(2048, nside(11));
    assert_eq!(4096, nside(12));
    assert_eq!(8192, nside(13));
    assert_eq!(16384, nside(14));
    assert_eq!(32768, nside(15));
    assert_eq!(65536, nside(16));
    assert_eq!(131072, nside(17));
    assert_eq!(262144, nside(18));
    assert_eq!(524288, nside(19));
    assert_eq!(1048576, nside(20));
    assert_eq!(2097152, nside(21));
    assert_eq!(4194304, nside(22));
    assert_eq!(8388608, nside(23));
    assert_eq!(16777216, nside(24));
    assert_eq!(33554432, nside(25));
    assert_eq!(67108864, nside(26));
    assert_eq!(134217728, nside(27));
    assert_eq!(268435456, nside(28));
    assert_eq!(536870912, nside(29));
  }

  #[test]
  #[should_panic]
  fn testpanic_nside() {
    nside(30);
  }

  #[test]
  fn testok_proj() {
    fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
      f64::sqrt((p2.0 - p1.0) * (p2.0 - p1.0) + (p2.1 - p1.1) * (p2.1 - p1.1))
    }
    // println!("{:?}", proj(0.0, TRANSITION_LATITUDE));
    // println!("{}", dist((0.0, 1.0), proj(0.0, TRANSITION_LATITUDE)));
    assert!(dist((0.0, 0.0), proj(0.0 * FRAC_PI_2, 0.0)) < 1e-15);
    assert!(dist((0.0, 1.0), proj(0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((1.0, 2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((2.0, 1.0), proj(1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((3.0, 2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((4.0, 1.0), proj(2.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((5.0, 2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((6.0, 1.0), proj(3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((7.0, 2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((0.0, 1.0), proj(4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((0.0, 0.0), proj(4.0 * FRAC_PI_2, 0.0)) < 1e-15);
    assert!(dist((0.0, -1.0), proj(0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((1.0, -2.0), proj(0.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((2.0, -1.0), proj(1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((3.0, -2.0), proj(1.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((4.0, -1.0), proj(2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((5.0, -2.0), proj(2.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((6.0, -1.0), proj(3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((7.0, -2.0), proj(3.0 * FRAC_PI_2 + FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((0.0, -1.0), proj(4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);

    assert!(dist((-0.0, 0.0), proj(-0.0 * FRAC_PI_2, 0.0)) < 1e-15);
    assert!(dist((-0.0, 1.0), proj(-0.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-1.0, 2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((-2.0, 1.0), proj(-1.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-3.0, 2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((-4.0, 1.0), proj(-2.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-5.0, 2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((-6.0, 1.0), proj(-3.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-7.0, 2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, FRAC_PI_2)) < 1e-15);
    assert!(dist((-0.0, 1.0), proj(-4.0 * FRAC_PI_2, TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-0.0, 0.0), proj(-4.0 * FRAC_PI_2, 0.0)) < 1e-15);
    assert!(dist((-0.0, -1.0), proj(-0.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-1.0, -2.0), proj(-0.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((-2.0, -1.0), proj(-1.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-3.0, -2.0), proj(-1.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((-4.0, -1.0), proj(-2.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-5.0, -2.0), proj(-2.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((-6.0, -1.0), proj(-3.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
    assert!(dist((-7.0, -2.0), proj(-3.0 * FRAC_PI_2 - FRAC_PI_4, -FRAC_PI_2)) < 1e-15);
    assert!(dist((-0.0, -1.0), proj(-4.0 * FRAC_PI_2, -TRANSITION_LATITUDE)) < 1e-15);
  }

  #[test]
  #[should_panic]
  fn testpanic_proj_1() {
    proj(0.0, -1.58);
  }

  #[test]
  #[should_panic]
  fn testpanic_proj_2() {
    proj(3.14159, -1.58);
  }

  #[test]
  fn testok_shs() {
    let ang_dist = 0.24;
    let shs = to_squared_half_segment(ang_dist);
    let ang_dist_2 = sphe_dist(shs);
    assert!((ang_dist - ang_dist_2).abs() < 1e-4);
  }

  #[test]
  fn testok_largest_center_to_vertex_distance() {
    let depth = 8;
    let d1_mas = largest_center_to_vertex_distance(depth, 0.0, 0.0)
      .to_degrees() * 3600.0;
    let d2_mas = largest_center_to_vertex_distance(depth, 0.0, LAT_OF_SQUARE_CELL)
      .to_degrees() * 3600.0;
    let d3_mas = largest_center_to_vertex_distance(depth, 0.0, TRANSITION_LATITUDE)
      .to_degrees() * 3600.0;
    let d4_mas = largest_center_to_vertex_distance(depth, PI / 4.0, TRANSITION_LATITUDE)
      .to_degrees() * 3600.0;
    let d5_mas = largest_center_to_vertex_distance(depth, 0.0, 0.5 * LAT_OF_SQUARE_CELL)
      .to_degrees() * 3600.0;
    let d6_mas = largest_center_to_vertex_distance(depth, 0.0, 0.5 * TRANSITION_LATITUDE)
      .to_degrees() * 3600.0;
    println!("d1: {}", d1_mas); // center parabol
    println!("d2: {}", d2_mas); // smallest value
    println!("d3: {}", d3_mas); // largest in PC (polar caps)
    println!("d4: {}", d4_mas); // smallest in PC = larges in EQR
    println!("d5: {} ~620?", d5_mas); // milieu parabol
    println!("d6: {} ~660?", d6_mas); // milieu droite eqr

    assert!((d1_mas - 632.8125000000000).abs() < 1e-12);
    assert!((d2_mas - 583.0213772328785).abs() < 1e-12);
    assert!((d3_mas - 861.2025252838432).abs() < 1e-12);
    assert!((d4_mas - 720.3786531802374).abs() < 1e-12);
    assert!((d5_mas - 620.3647193082197).abs() < 1e-12);
    assert!((d6_mas - 591.2475292479544).abs() < 1e-12);
  }
  
}
