use std::sync::Once;
use std::ops::Shr;

use super::compass_point::*;
use super::external_edge::*;
use super::*;

/// Array storing pre-computed values for each of the 30 possible depth (from 0 to 29)
/// Info: Unfortunately Default::default(); do no work with static arrays :o/
static mut LAYERS: [Option<Layer>; 30] =  [
  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
];
/// See the get_or_create function, each object is used for the lazy instantiation of the 
/// layer of the corresponding depth.
/// Info: Unfortunately Default::default(); do no work with static arrays :o/
static LAYERS_INIT: [Once; 30] = [
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(), Once::new(),
  Once::new(), Once::new()
];

// const l0: Layer = 

/// Lazy factory method: instantiate a new Layer at the first call for a given depth; after the 
/// first call, returns an already instantiated Layer.
/// # Info
/// This method resort to a double-checked lock, ensuring thread-safety.
pub fn get_or_create(depth: u8) -> &'static Layer {
  unsafe {
    // Inspired from the Option get_or_insert_with method, modified to ensure thread safety with
    // https://doc.rust-lang.org/std/sync/struct.Once.html
    // This implements a double-checked lock
    match LAYERS[depth as usize] {
      Some(ref v) =>  return v,
      None => {
        LAYERS_INIT[depth as usize].call_once(|| {
          LAYERS[depth as usize] = Some(Layer::new(depth));
        });
      },
    }
    match LAYERS[depth as usize] {
      Some(ref v) => v,
      _ => unreachable!(),
    }
  }
}

pub const fn to_range(hash: u64, delta_depth: u8) -> std::ops::Range<u64> {
  let twice_delta_depth = delta_depth << 1;
  (hash << twice_delta_depth)..((hash + 1) << twice_delta_depth)
}


/// Transforms the given NESTED hash value into its uniq representation, i.e. the depth
/// is encoded together with the hash value such that each possible (deph, hash) pair 
/// gives a unique number.
/// In practice, the unique representation uses a sentinel bit (set to one) to code the depth.
/// The sentinel bit (set to 1) is the least significant bit (LSB) among the unused bits,
/// i.e. the most significant bit (MSB) located just after the hash MSB.
/// Said differently, the sentinel bit is the `(1 + 4 + 2*depth)^th` MSB
/// The encoding, in the case of the nested scheme, is thus `0...0sbbbb112233...`, with:
/// * `0...0`: unused bits</li>
/// * `s` : sentinel bit</li>
/// * `bbbb`: the 4 bits coding the base cell
/// * `11`: the 2 bits coding depth 1
/// * `22`: the 2 bits coding depth 2
/// * `33`: the 2 bits coding depth 3
/// * ...
/// 
/// # Example
/// 
/// ```rust
/// use cdshealpix::nested::{get_or_create, Layer};
/// let l0 = get_or_create(0);
/// assert_eq!(l0.to_uniq(0), 16);
/// ```
#[inline]
pub fn to_uniq(depth: u8, hash: u64) -> u64 {
  check_depth(depth);
  to_uniq_unsafe(depth, hash)
}
#[inline]
const fn to_uniq_unsafe(depth: u8, hash: u64) -> u64 {
  (16_u64 << (depth << 1)) | hash
}

/// Same as [to_uniq](fn.to_uniq.html), but
/// following the [IVOA](http://ivoa.net/documents/MOC/) convention.
/// It does not rely on a sentinel bit and use one less bit.
#[inline]
pub fn to_uniq_ivoa(depth: u8, hash: u64) -> u64 {
  check_depth(depth);
  to_uniq_ivoa_unsafe(depth, hash)
}
#[inline]
const fn to_uniq_ivoa_unsafe(depth: u8, hash: u64) -> u64 {
  (4_u64 << (depth << 1)) + hash
}

/// Returns the depth and the hash number from the uniq representation.
/// Inverse operation of [to_uniq](fn.to_uniq.html).
pub fn from_uniq(uniq_hash: u64) -> (u8, u64) {
  let depth = (60 - uniq_hash.leading_zeros()) >> 1;
  let hash = uniq_hash & !(16_u64 << (depth << 1));
  // = uniq_hash - (16_u64 << (depth << 1));
  // = uniq_hash & !highest_one_bit(uniq_hash)); also works, but do not benefit from depth
  (depth as u8, hash)
}

/// Returns the depth and the hash number from the uniq representation.
/// Inverse operation of [to_uniq](fn.to_uniq_ivoa.html).
pub fn from_uniq_ivoa(uniq_hash: u64) -> (u8, u64) {
  let depth = (61 - uniq_hash.leading_zeros()) >> 1;
  let hash = uniq_hash - (4_u64 << (depth << 1));
  (depth as u8, hash)
}

const fn highest_one_bit(mut i: u64) -> u64 {
  i |= i >>  1;
  i |= i >>  2;
  i |= i >>  4;
  i |= i >>  8;
  i |= i >> 16;
  i |= i >> 32;
  i - (i >> 1)
}


/// Conveniency function returning the number of hash value in the [Layer] of the given *depth*.
#[inline]
pub fn n_hash(depth: u8) -> u64 {
  get_or_create(depth).n_hash
}

/// Conveniency function simply calling the [hash](struct.Layer.html#method.hash) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn hash(depth: u8, lon: f64, lat: f64) -> u64 {
  get_or_create(depth).hash(lon, lat)
}

/// Conveniency function simply calling the [hash](struct.Layer.html#method.hash_with_dxdy) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn hash_with_dxdy(depth: u8, lon: f64, lat: f64) -> (u64, f64, f64) {
  get_or_create(depth).hash_with_dxdy(lon, lat)
}

/// Conveniency function simply calling the [center](struct.Layer.html#method.sph_coo) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn sph_coo(depth: u8, hash: u64, dx: f64, dy: f64) -> (f64, f64) { 
  get_or_create(depth).sph_coo(hash, dx, dy)
}

/// Conveniency function simply calling the [center](struct.Layer.html#method.center) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn center(depth: u8, hash: u64) -> (f64, f64) {
  get_or_create(depth).center(hash)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.vertices) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn vertices(depth: u8, hash: u64) -> [(f64, f64); 4] {
  get_or_create(depth).vertices(hash)
}


/// Conveniency function simply calling the [vertices](struct.Layer.html#method.path_along_cell_side) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn path_along_cell_side(depth: u8, hash: u64, from_vertex:  &Cardinal, to_vertex: &Cardinal, include_to_vertex: bool, n_segments: u32) -> Box<[(f64, f64)]> {
  get_or_create(depth).path_along_cell_side(hash, from_vertex, to_vertex, include_to_vertex, n_segments)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.path_along_cell_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn path_along_cell_edge(depth: u8, hash: u64, starting_vertex: &Cardinal, clockwise_direction: bool, n_segments_by_side: u32) -> Box<[(f64, f64)]> {
  get_or_create(depth).path_along_cell_edge(hash, starting_vertex, clockwise_direction, n_segments_by_side)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.grid) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn grid(depth: u8, hash: u64, n_segments_by_side: u16) -> Box<[(f64, f64)]> {
  get_or_create(depth).grid(hash, n_segments_by_side)
}

/// Conveniency function simply calling the [neighbours](struct.Layer.html#method.neighbours) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn neighbours(depth: u8, hash: u64, include_center: bool) -> MainWindMap<u64> {
  get_or_create(depth).neighbours(hash, include_center)
}

/// Conveniency function simply calling the [internal_edge](struct.Layer.html#method.internal_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn internal_edge(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  assert!(depth + delta_depth < DEPTH_MAX);
  Layer::internal_edge(hash, delta_depth)
}

/// Conveniency function simply calling the [internal_edge_sorted](struct.Layer.html#method.internal_edge_sorted) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn internal_edge_sorted(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  assert!(depth + delta_depth < DEPTH_MAX);
  Layer::internal_edge_sorted(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge_struct](struct.Layer.html#method.external_edge_struct) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge_struct(depth: u8, hash: u64, delta_depth: u8) -> ExternalEdge {
  get_or_create(depth).external_edge_struct(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge](struct.Layer.html#method.external_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  get_or_create(depth).external_edge(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge_sorted](struct.Layer.html#method.external_edge_sorted) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge_sorted(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  get_or_create(depth).external_edge_sorted(hash, delta_depth)
}

/// Conveniency function simply calling the [bilinear_interpolation](struct.Layer.html#method.bilinear_interpolation) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn bilinear_interpolation(depth: u8, lon: f64, lat: f64) -> [(u64, f64); 4] {
  get_or_create(depth).bilinear_interpolation(lon, lat)
}

/// Conveniency function simply calling the [cone_coverage_approx](struct.Layer.html#method.cone_coverage_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_approx(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get_or_create(depth).cone_coverage_approx(cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_coverage_flat](struct.Layer.html#method.cone_coverage_approx) method
/// of the [Layer] of the given *depth* and retrieving a flat array.
#[inline]
pub fn cone_coverage_approx_flat(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> Box<[u64]> {
  get_or_create(depth).cone_coverage_approx(cone_lon, cone_lat, cone_radius).to_flat_array()
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.cone_coverage_approx_custom) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_approx_custom(depth: u8, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get_or_create(depth).cone_coverage_approx_custom(delta_depth, cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [elliptical_cone_coverage](struct.Layer.html#method.elliptical_cone_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn elliptical_cone_coverage(depth: u8, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
  get_or_create(depth).elliptical_cone_coverage( lon, lat, a, b, pa)
}

/// Conveniency function simply calling the [elliptical_cone_coverage_custom](struct.Layer.html#method.elliptical_cone_coverage_custom) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn elliptical_cone_coverage_custom(depth: u8, delta_depth: u8, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
  get_or_create(depth).elliptical_cone_coverage_custom(delta_depth, lon, lat, a, b, pa)
}

/// Conveniency function simply calling the [polygon_coverage_approx](struct.Layer.html#method.polygon_coverage_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn polygon_coverage(depth: u8, vertices: &[(f64, f64)], exact_solution: bool) -> BMOC {
  get_or_create(depth).polygon_coverage(vertices, exact_solution)
}


pub mod bmoc;
pub mod zordercurve;
mod gpu;

use self::zordercurve::{ZOrderCurve, get_zoc};
use self::bmoc::*;
use super::ring::{triangular_number_x4};
use super::sph_geom::coo3d::*;
use super::sph_geom::{Polygon};
use super::sph_geom::cone::{Cone};
use super::sph_geom::elliptical_cone::EllipticalCone;
use super::{proj, direction_from_neighbour, edge_cell_direction_from_neighbour};
use super::compass_point::{Cardinal, CardinalSet, CardinalMap};
use super::compass_point::MainWind::{S, SE, E, SW, C, NE, W, NW, N};
use super::special_points_finder::{arc_special_points};


/// Defines an HEALPix layer in the NESTED scheme.
/// A layer is simply an utility structure containing all constants and methods related
/// to a given depth.
pub struct Layer { // Why not creating all of them at compilation using a macro?!
  depth: u8,
  nside: u32,
  nside_minus_1: u32,
  n_hash: u64,
  twice_depth: u8,
  d0h_mask: u64,
  x_mask: u64,
  y_mask: u64,
  xy_mask: u64,
  nside_remainder_mask: u64, // = nside - 1
  time_half_nside: i64,
  one_over_nside: f64,
  z_order_curve: &'static dyn ZOrderCurve,
}

impl Layer {
  
  fn new(depth: u8) -> Layer {
    let twice_depth: u8 = depth << 1u8;
    let nside: u32 = 1u32 << depth;
    let mut x_mask = 0u64;
    let mut xy_mask = 0u64;
    let mut time_half_nside = -1_i64 << 52;
    if depth > 0 {
      x_mask = 0x5555555555555555u64 >> (64 - twice_depth); // ...0101
      xy_mask = (1u64 << twice_depth) - 1u64;
      time_half_nside = ((depth - 1) as i64) << 52;
    }
    Layer {
      depth,
      nside,
      nside_minus_1: nside - 1,
      n_hash: super::n_hash_unsafe(depth),
      twice_depth,
      d0h_mask: 15_u64 << twice_depth,
      x_mask,
      y_mask: x_mask << 1,
      xy_mask,
      nside_remainder_mask: xy_mask >> depth, // = nside - 1
      time_half_nside,
      one_over_nside: 1f64 / super::nside_unsafe(depth) as f64,
      z_order_curve: get_zoc(depth)
    }
  }

  /// Returns the depth of the Layer (i.e. the HEALPix *order*)
  #[inline]
  pub fn depth(&self) -> u8 { self.depth }

  /// Returns the number of hash value of the Layer, i.e. the number of cells.
  #[inline]
  pub fn n_hash(&self) -> u64 {
    self.n_hash
  }
  
  /// Returns the cell number (hash value) associated with the given position on the unit sphere
  /// # Inputs
  /// - `lon`: longitude in radians, support reasonably large positive and negative values
  ///          producing accurate results with a naive range reduction like modulo 2*pi
  ///          (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
  /// - `lat`: latitude in radians, must be in `[-pi/2, pi/2]`
  /// # Output
  /// - the cell number (hash value) associated with the given position on the unit sphere,
  ///   in `[0, 12*nside^2[`
  /// # Panics
  ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
  /// # Examples
  /// ```rust
  /// use cdshealpix::{nside};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 12_u8;
  /// let nside = nside(depth) as u64;
  /// let nested12: &Layer = get_or_create(depth);
  /// assert_eq!(nside * nside - 1, nested12.hash(12.5_f64.to_radians(), 89.99999_f64.to_radians()));
  /// ```
  pub fn hash(&self, lon: f64, lat: f64) -> u64 {
    self.hash_v2(lon, lat)
  }

  /* The clean way to do, but changes the API...
  pub fn hash(&self, lon: f64, lat: f64) -> Result<u64, Error> {
    check lon, lat
    hash_uncheked(lon, lat)
  }
  pub fn hash_uncheked(&self, lon: f64, lat: f64) -> u64 {
    ...
  }*/
  
  pub fn hash_v1(&self, lon: f64, lat: f64) -> u64 {
    let mut xy = proj(lon, lat);
    xy.0 = ensures_x_is_positive(xy.0);
    self.shift_rotate_scale(&mut xy);
    let mut ij = discretize(xy);
    let ij_d0c = self.base_cell_coos(&ij);
    let d0h_bits = self.depth0_bits(ij_d0c.0, ij_d0c.1, ij, xy/*, lon, lat*/);
    self.to_coos_in_base_cell(&mut ij);
    self.build_hash(d0h_bits, ij.0 as u32, ij.1 as u32)
  }

  
  
  pub fn hash_v2(&self, lon: f64, lat: f64) -> u64 {
    check_lat(lat);
    let (d0h, l_in_d0c, h_in_d0c) = Layer::d0h_lh_in_d0c(lon, lat);
    // Coords inside the base cell
    //  - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
    let i = f64::from_bits((self.time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64) as u32;
    let j = f64::from_bits((self.time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64) as u32;
    //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
    let i = if i == self.nside { self.nside_minus_1 } else { i };
    let j = if j == self.nside { self.nside_minus_1 } else { j };
    self.build_hash_from_parts(d0h, i, j)
  }
  
  pub fn hash_dxdy_v2(&self, lon: f64, lat: f64) -> (u64, f64, f64) {
    let (d0h, l_in_d0c, h_in_d0c) = Layer::d0h_lh_in_d0c(lon, lat);
    // Coords inside the base cell time nside/2
    let x = f64::from_bits((self.time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64); 
    debug_assert!(x <= (0.0 - 1e-15) || x < (self.nside as f64 + 1e-15), format!("x: {}, x_proj: {}; y_proj: {}", &x, &h_in_d0c, &l_in_d0c));
    let y = f64::from_bits((self.time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64);
    debug_assert!(y <= (0.0 - 1e-15) || y < (self.nside as f64 + 1e-15), format!("y: {}, x_proj: {}; y_proj: {}", &y, &h_in_d0c, &l_in_d0c));
    // - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
    let i = x as u32;
    let j = y as u32;
    //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
    let i = if i == self.nside { self.nside_minus_1 } else { i };
    let j = if j == self.nside { self.nside_minus_1 } else { j };
    (
      self.build_hash_from_parts(d0h, i, j),
      (x - (i as f64)),
      (y - (j as f64)),
    )
  }

  #[inline]
  fn d0h_lh_in_d0c(lon: f64, lat: f64) -> (u8, f64, f64) {
    let (x_pm1, q) = Layer::xpm1_and_q(lon);
    if lat > TRANSITION_LATITUDE {
      // North polar cap, Collignon projection.
      // - set the origin to (PI/4, 0)
      let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 + PI_OVER_FOUR).cos();
      let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
      let d0h = q;
      (d0h, x_proj, y_proj)
    } else if lat < -TRANSITION_LATITUDE {
      // South polar cap, Collignon projection
      // - set the origin to (PI/4, -PI/2)
      let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 - PI_OVER_FOUR).cos(); // cos(-x) = cos(x)
      let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, sqrt_3_one_min_z);
      let d0h = q + 8;
      (d0h, x_proj, y_proj)
    } else {
      // Equatorial region, Cylindrical equal area projection
      // - set the origin to (PI/4, 0)               if q = 2
      // - set the origin to (PI/4, -PI/2)           if q = 0
      // - set the origin to (0, -TRANSITION_LAT)    if q = 3
      // - set the origin to (PI/2, -TRANSITION_LAT) if q = 1
      // let zero_or_one = (x_cea as u8) & 1;
      let y_pm1 = lat.sin() * ONE_OVER_TRANSITION_Z;
      // Inequalities have been carefully chosen so that S->E and S->W axis are part of the cell,
      // and not E->N and W->N
      
      // Version with branch
      // |\3/|
      // .2X1.
      // |/0\|
      /*let q13 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */  debug_assert!(q12 == 0 || q12 == 1);
      let q23 = (x_pm1 <=  y_pm1) as u8; /* 1/0 */  debug_assert!(q23 == 0 || q23 == 1);
      match q13 | (q23 << 1) {
        0 => ( q         , x_pm1      , y_pm1 + 2.0),
        1 => ((q + 5) & 7, x_pm1 - 1.0, y_pm1 + 1.0), // (q + 5) & 7 <=> (q + 1) | 4
        2 => ( q + 4     , x_pm1 + 1.0, y_pm1 + 1.0),
        3 => ( q + 8     , x_pm1      , y_pm1),
        _ => unreachable!(),
      }*/
      // Branch free version
      // |\2/|
      // .3X1.
      // |/0\|
      let q01 = (x_pm1 >   y_pm1) as u8;  /* 0/1 */  debug_assert!(q01 == 0 || q01 == 1);
      let q12 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */  debug_assert!(q12 == 0 || q12 == 1);
      let q1 = q01 & q12; /* = 1 if q1, 0 else */       debug_assert!( q1 == 0 ||  q1 == 1);
      let q013 = q01 + (1 - q12); // = q01 + q03; /* 0/1 + 1\0 +  */
      // x: x_pm1 + 1 if q3 | x_pm1 - 1 if q1 | x_pm1 if q0 or q2
      let x_proj = x_pm1 - ((q01 + q12) as i8 - 1) as f64;
      // y: y_pm1 + 0 if q2 | y_pm1 + 1 if q1 or q3 | y_pm1 + 2 if q0 
      let y_proj = y_pm1 + q013 as f64;
      // d0h: +8 if q0 | +4 if q3 | +5 if q1
      let d0h = (q013 << 2) + ((q + q1) & 3);
      (d0h, x_proj, y_proj)
    }
  }
  
  /// Transform the input longitude, in radians, in a value `x` in `[-1, 1[` plus a quarter in `[0, 3]`,
  /// such that `lon = (x + 1) * PI / 4 + q * PI / 2`. 
  #[inline]
  fn xpm1_and_q(lon: f64) -> (f64, u8) {
    let lon_bits = lon.to_bits();
    let lon_abs = f64::from_bits(lon_bits & F64_BUT_SIGN_BIT_MASK);
    let lon_sign = lon_bits & F64_SIGN_BIT_MASK;
    let x = lon_abs * FOUR_OVER_PI;
    let q = (x as u8 | 1_u8) & 7_u8;    debug_assert!(q < 8);
    // Remark: to avoid the branch, we could have copied lon_sign on x - q, 
    //         but I so far lack of idea to deal with q efficiently.
    //         And we are not supposed to have negative longitudes in ICRS 
    //         (the most used reference system in astronomy).
    if lon_sign == 0 { // => lon >= 0
      (x - (q as f64), q >> 1)
    } else { // case lon < 0 should be rare => few risks of branch miss-prediction
      // Since q in [0, 3]: 3 - (q >> 1)) <=> 3 & !(q >> 1)
      // WARNING: BE SURE TO HANDLE THIS CORRECTLY IN THE REMAINING OF THE CODE!
      //  - Case lon =  3/4 pi = 270 deg => x = -1, q=3
      //  - Case lon = -1/2 pi = -90 deg => x =  1, q=2
      (q as f64 - x, 3 - (q >> 1))
    }
  }
  
  
  /// Returns the cell number (hash value) associated with the given position on the unit sphere, 
  /// together with the offset `(dx, dy)` on the Euclidean plane of the projected position with
  /// respect to the origin of the cell (South vertex).
  /// # Inputs
  /// - `lon`: longitude in radians, support reasonably large positive and negative values
  ///          producing accurate results with a naive range reduction like modulo 2*pi
  ///          (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
  /// - `lat`: latitude in radians, must be in `[-pi/2, pi/2]`
  /// # Output
  /// - the cell number (hash value) associated with the given position on the unit sphere,
  ///   in `[0, 12*nside^2[`
  /// - `dx`: the positional offset $\in [0, 1[$ along the south-to-east axis
  /// - `dy`: the positional offset $\in [0, 1[$ along the south-to-west axis
  /// # Panics
  ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
  /// # Examples
  /// ```rust
  /// use cdshealpix::{nside};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 12_u8;
  /// let nside = nside(depth) as u64;
  /// let nested12: &Layer = get_or_create(depth);
  /// let h_org = nside * nside - 1;
  /// let (h_ra, h_dec) = nested12.center(h_org);
  /// let (h, dx, dy) = nested12.hash_with_dxdy(h_ra, h_dec);
  /// assert_eq!(h_org, h);
  /// assert_eq!(0.5, dx);
  /// assert_eq!(0.5, dy);
  /// ```
  pub fn hash_with_dxdy(&self, lon: f64, lat: f64) -> (u64, f64, f64) {
    let mut xy = proj(lon, lat);
    xy.0 = ensures_x_is_positive(xy.0);
    self.shift_rotate_scale(&mut xy);
    let mut ij = discretize(xy);
    let dx = xy.0 - (ij.0 as f64);
    let dy = xy.1 - (ij.1 as f64);
    let ij_d0c = self.base_cell_coos(&ij);
    let d0h_bits = self.depth0_bits(ij_d0c.0, ij_d0c.1, ij, xy/*, lon, lat*/);
    self.to_coos_in_base_cell(&mut ij);
    (self.build_hash(d0h_bits, ij.0 as u32, ij.1 as u32), dx, dy)
  }
  
  #[inline]
  fn shift_rotate_scale(&self, xy: &mut (f64, f64)) {
    let (ref mut x, ref mut y) = *xy;
    let tmp = 8.0 - *x;
    *y += 1.0;
    *x = f64::from_bits((self.time_half_nside + f64::to_bits(*x + *y) as i64) as u64);
    *y = f64::from_bits((self.time_half_nside + f64::to_bits(*y + tmp) as i64) as u64);
  }

  #[inline]
  fn base_cell_coos(&self, ij: &(u64, u64)) -> (u8, u8) {
    (self.div_by_nside_floor_u8((*ij).0), self.div_by_nside_floor_u8((*ij).1))
  }

  #[inline]
  fn div_by_nside_floor_u8(&self, val: u64) -> u8 {
    (val >> self.depth) as u8
  }

  #[inline]
  fn to_coos_in_base_cell(&self, ij: &mut (u64, u64)) {
    let (ref mut i, ref mut j) = *ij;
    *i = self.modulo_nside(*i);
    *j = self.modulo_nside(*j);
  }

  #[inline]
  fn modulo_nside(&self, val: u64) -> u64 {
    val & self.nside_remainder_mask
  }

  #[inline]
  fn build_hash_from_parts(&self, d0h: u8, i: u32, j: u32) -> u64 {
    self.build_hash((d0h as u64) << self.twice_depth, i, j)
  }

  #[inline]
  fn build_hash_from_parts_opt(&self, d0h: u8, i: u32, j: u32) -> Option<u64> {
    Some(self.build_hash_from_parts(d0h, i, j))
  }

  #[inline]
  fn build_hash(&self, d0h_bits: u64, i: u32, j: u32) -> u64 {
    debug_assert!(i < self.nside && j < self.nside);
    d0h_bits | self.z_order_curve.ij2h(i, j)
  }

  ///
  #[inline]
  fn depth0_bits(&self, i: u8, j: u8, mut ij: (u64, u64), xy: (f64, f64)/*, lon: f64, lat: f64*/) -> u64 {
    // self.base_hash_bits_lupt[i as usize][j as usize]
    let k = 5_i8 - (i + j) as i8;
    // The two branches -2 and -1 are extremely rare (north pole, NPC cells upper NE and NW borders),
    // so few risks of branch miss-prediction.
    // println!("k: {}; i: {}; j: {}; ij: {:?}; xy: {:?}", &k, &i, &j, &ij, &xy);
    match k {
      0..=2 => (((k << 2) + ( ((i as i8) + ((k - 1) >> 7)) & 3_i8)) as u64) << self.twice_depth,
      -1 => {
        if xy.0 - ij.0 as f64 > xy.1 - ij.1 as f64 {
          ((((i - 1u8) & 3u8) as u64) << self.twice_depth) | self.y_mask
        } else {
          ((((i + 2u8) & 3u8) as u64) << self.twice_depth) | self.x_mask
        }
      },
      -2 => (((i - 2u8) as u64) << self.twice_depth) | self.xy_mask,
      3 => { // rare case due to lack of numerical precision 
        let d0 = (xy.0 - ij.0 as f64).abs();
        let d1 = (xy.1 - ij.1 as f64).abs();
        if d0 < d1 {
          ij.0 += 1;
          self.depth0_bits(i + 1_u8, j, ij, xy)
        } else {
          ij.1 += 1;
          self.depth0_bits(i, j + 1_u8, ij, xy)
        }
      },
      4 => { // rare case due to lack of numerical precision 
        ij.0 += 1;
        ij.1 += 1;
       self.depth0_bits(i + 1_u8, j + 1_u8, ij, xy)
      },
      /*_ => panic!("Algorithm error: case k = {} not supported! depth: {}, lon: {}, lat: {}, x: {}, y: {}", 
                  k, self.depth, lon, lat, xy.0, xy.1),*/
      _ => panic!("Algorithm error: case k = {} not supported!", k),
    }
  }

  /// Conveniency function simply calling the [hash](fn.to_uniq.html) method with this layer *depth*.
  pub fn to_uniq(&self, hash: u64) -> u64 {
    // depth already tested, so we call the unsafe method
    nested::to_uniq_unsafe(self.depth, hash)
  }

  /// Conveniency function simply calling the [hash](fn.to_uniq_ivoa.html) method with this layer *depth*.
  pub fn to_uniq_ivoa(&self, hash: u64) -> u64 {
    // depth already tested, so we call the unsafe method
    nested::to_uniq_ivoa_unsafe(self.depth, hash)
  }
  
  /// Transforms the given NESTED hash value into the RING hash value.
  /// 
  /// # Examples
  /// 
  /// At depth 0, no differences:
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// let depth = 0;
  /// let n0 = get_or_create(depth);
  /// for h in 0..12 {
  ///   assert_eq!(n0.to_ring(h), h);
  /// }
  /// ```
  /// 
  /// At depth 1:
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// 
  /// let depth = 1;
  /// let n1 = get_or_create(depth);
  /// 
  /// assert_eq!(n1.to_ring(0),  13);
  /// assert_eq!(n1.to_ring(1),   5);
  /// assert_eq!(n1.to_ring(2),   4);
  /// assert_eq!(n1.to_ring(3),   0);
  /// assert_eq!(n1.to_ring(4),  15);
  /// assert_eq!(n1.to_ring(5),   7);
  /// assert_eq!(n1.to_ring(6),   6);
  /// assert_eq!(n1.to_ring(7),   1);
  /// assert_eq!(n1.to_ring(8),  17);
  /// assert_eq!(n1.to_ring(9),   9);
  /// assert_eq!(n1.to_ring(10),  8);
  /// assert_eq!(n1.to_ring(11),  2);
  /// assert_eq!(n1.to_ring(12), 19);
  /// assert_eq!(n1.to_ring(13), 11);
  /// assert_eq!(n1.to_ring(14), 10);
  /// assert_eq!(n1.to_ring(15),  3);
  /// assert_eq!(n1.to_ring(16), 28);
  /// assert_eq!(n1.to_ring(17), 20);
  /// assert_eq!(n1.to_ring(18), 27);
  /// assert_eq!(n1.to_ring(19), 12);
  /// assert_eq!(n1.to_ring(20), 30);
  /// assert_eq!(n1.to_ring(21), 22);
  /// assert_eq!(n1.to_ring(22), 21);
  /// assert_eq!(n1.to_ring(23), 14);
  /// assert_eq!(n1.to_ring(24), 32);
  /// assert_eq!(n1.to_ring(25), 24);
  /// assert_eq!(n1.to_ring(26), 23);
  /// assert_eq!(n1.to_ring(27), 16);
  /// assert_eq!(n1.to_ring(28), 34);
  /// assert_eq!(n1.to_ring(29), 26);
  /// assert_eq!(n1.to_ring(30), 25);
  /// assert_eq!(n1.to_ring(31), 18);
  /// assert_eq!(n1.to_ring(32), 44);
  /// assert_eq!(n1.to_ring(33), 37);
  /// assert_eq!(n1.to_ring(34), 36);
  /// assert_eq!(n1.to_ring(35), 29);
  /// assert_eq!(n1.to_ring(36), 45);
  /// assert_eq!(n1.to_ring(37), 39);
  /// assert_eq!(n1.to_ring(38), 38);
  /// assert_eq!(n1.to_ring(39), 31);
  /// assert_eq!(n1.to_ring(40), 46);
  /// assert_eq!(n1.to_ring(41), 41);
  /// assert_eq!(n1.to_ring(42), 40);
  /// assert_eq!(n1.to_ring(43), 33);
  /// assert_eq!(n1.to_ring(44), 47);
  /// assert_eq!(n1.to_ring(45), 43);
  /// assert_eq!(n1.to_ring(46), 42);
  /// assert_eq!(n1.to_ring(47), 35);
  /// ```
  /// 
  /// At depth 2 (non exhaustive test):
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// 
  /// let depth = 2;
  /// let n2 = get_or_create(depth);
  /// /// NPC
  /// assert_eq!(n2.to_ring(47),  2);
  /// assert_eq!(n2.to_ring(29),  7);
  /// assert_eq!(n2.to_ring(60), 22);
  /// // EQR
  /// assert_eq!(n2.to_ring(51),   54);
  /// assert_eq!(n2.to_ring(88),  107);
  /// assert_eq!(n2.to_ring(174), 129);
  /// // SPC
  /// assert_eq!(n2.to_ring(177), 187);
  /// assert_eq!(n2.to_ring(153), 157);
  /// assert_eq!(n2.to_ring(144), 189);
  /// ```
  pub fn to_ring(&self, hash: u64) -> u64 {
    // Number of isolatitude rings in a base cell: 
    //    nbr rings:   nr = 2 * nside - 1 
    // => index max: hmax = 2 * nside - 2
    // Index of ring at the NPC / EQR interface:
    //   i = nside - 1                                       (number of rings =     nside)
    // Index of ring at lat=0 (in the EQR): 
    //   i = nr = hmax + 1 =  = 2 * nside - 1 => always odd  (number of rings = 2 * nside)
    // Index of ring at the EQR / SPC interface:
    //   i = nr + nside = 3 * nside - 1                      (number of rings = 3 * nside)
    // We note h = x + y (<=> rotation 45 and scale sqrt(2))
    // North polar cap   base cells:   i =  hmax - h              = 2 * nside - 2 - (x + y)
    // Equatorial region base cells:   i = (hmax - h) +     nside = 3 * nside - 2 - (x + y)
    // South polar cap   base cells:   i = (hmax - h) + 2 * nside = 4 * nside - 2 - (x + y)
    let HashParts {d0h, i, j} = self.decode_hash(hash);
    let h: u64 = i as u64 + j as u64;
    let l: i64 = i as i64 - j as i64;
    let i_d0h = div4_remainder(d0h) as u64;
    let j_d0h = div4_quotient(d0h) as u64;    debug_assert!(j_d0h <= 2);
    let i_ring: u64 = self.nside_time(j_d0h + 2) - (h + 2);
    // Number of elements in isolatitude ring of index i (i in [0, 2*(2*nside - 1)]):
    // North polar cap: if (i < nside)         nj = 4 * i
    // South polar cap: if (i >= 3*nside - 1)  nj = 4 * h = 4*((4*nside-2) - i) 
    // Equatorial regi: if (ns <= i < 3*ns-1)  nj = 4 * nside
    // l = x - y; In a base cell, l in [-nside+1, nside-1] => 2*nside - 1 values
    // EQR: j = l / 2 + nside/2 (if not equatorial cell) + nside*(ipix%4) (special case if l<0 && baseCell=4)
    // NPC: j = l / 2 + (i + 1) / 2 + (i+1)*(ipix%4)
    // SPC: j = l / 2 + (h + 1) / 2 + (h+1)*(ipix%4)
    let first_isolat_index;
    let mut i_in_ring = div2_quotient(l); // Quotient such that: 1/2 = 0; -1/2 = -1
    if i_ring < self.nside as u64 { // North polar cap + NPC/EQR tansition
      // sum from i = 1 to ringIndex of 4 * i = 4 * i*(i+1)/2 = 2 * i*(i+1)
      let ip1 = i_ring + 1;
      first_isolat_index = (i_ring * ip1) << 1;
      i_in_ring += (div2_quotient(ip1) + ip1 * i_d0h) as i64;
    } else if i_ring >= self.nside_time(3) - 1 { // South polar cap
      let ip1 = h + 1;
      first_isolat_index = self.n_hash - triangular_number_x4(ip1);
      i_in_ring += (div2_quotient(ip1) + ip1 * i_d0h) as i64;
    } else { // Equatorial region
      // sum from i = 1 to nside of i
      first_isolat_index = self.first_hash_in_eqr() + self.minus_nside_x_4nside(i_ring); 
      i_in_ring += div2_quotient(self.nside_time(div2_remainder(j_d0h + 1))) as i64;
      i_in_ring += self.nside_time(if d0h == 4 && l < 0 { 4 } else { i_d0h }) as i64;
    }
    i_in_ring as u64 + first_isolat_index
  }
  
  /// Transforms the given RING hash value into the NESTED hash value.
  /// 
  /// # WARNING
  /// The RING NSIDE parameter MUST match the NESTED NSIDE parameter!! 
  /// 
  /// # Examples
  /// 
  /// At depth 0, no differences:
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// let depth = 0;
  /// let n0 = get_or_create(depth);
  /// for h in 0..12 {
  ///   assert_eq!(n0.from_ring(h), h);
  /// }
  /// ```
  /// At depth 1:
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// 
  /// let depth = 1;
  /// let n1 = get_or_create(depth);
  /// 
  /// assert_eq!( 3, n1.from_ring(0));
  /// assert_eq!( 7, n1.from_ring(1));
  /// assert_eq!(11, n1.from_ring(2));
  /// assert_eq!(15, n1.from_ring(3));
  /// assert_eq!( 2, n1.from_ring(4));
  /// assert_eq!( 1, n1.from_ring(5));
  /// assert_eq!( 6, n1.from_ring(6));
  /// assert_eq!( 5, n1.from_ring(7));
  /// assert_eq!(10, n1.from_ring(8));
  /// assert_eq!( 9, n1.from_ring(9));
  /// assert_eq!(14, n1.from_ring(10));
  /// assert_eq!(13, n1.from_ring(11));
  /// assert_eq!(19, n1.from_ring(12));  
  /// assert_eq!( 0, n1.from_ring(13));
  /// assert_eq!(23, n1.from_ring(14));
  /// assert_eq!( 4, n1.from_ring(15));
  /// assert_eq!(27, n1.from_ring(16)); 
  /// assert_eq!( 8, n1.from_ring(17));
  /// assert_eq!(31, n1.from_ring(18));
  /// assert_eq!(12, n1.from_ring(19));
  /// assert_eq!(17, n1.from_ring(20));
  /// assert_eq!(22, n1.from_ring(21));
  /// assert_eq!(21, n1.from_ring(22));
  /// assert_eq!(26, n1.from_ring(23));
  /// assert_eq!(25, n1.from_ring(24));
  /// assert_eq!(30, n1.from_ring(25));
  /// assert_eq!(29, n1.from_ring(26));
  /// assert_eq!(18, n1.from_ring(27));
  /// assert_eq!(16, n1.from_ring(28));
  /// assert_eq!(35, n1.from_ring(29));
  /// assert_eq!(20, n1.from_ring(30));
  /// assert_eq!(39, n1.from_ring(31));
  /// assert_eq!(24, n1.from_ring(32));
  /// assert_eq!(43, n1.from_ring(33));
  /// assert_eq!(28, n1.from_ring(34));
  /// assert_eq!(47, n1.from_ring(35));
  /// assert_eq!(34, n1.from_ring(36));
  /// assert_eq!(33, n1.from_ring(37));
  /// assert_eq!(38, n1.from_ring(38));
  /// assert_eq!(37, n1.from_ring(39));
  /// assert_eq!(42, n1.from_ring(40));
  /// assert_eq!(41, n1.from_ring(41));
  /// assert_eq!(46, n1.from_ring(42));
  /// assert_eq!(45, n1.from_ring(43));
  /// assert_eq!(32, n1.from_ring(44));
  /// assert_eq!(36, n1.from_ring(45));
  /// assert_eq!(40, n1.from_ring(46));
  /// assert_eq!(44, n1.from_ring(47));
  /// ```
  /// 
  /// 
  /// At depth 2 (non exhaustive test):
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  /// 
  /// let depth = 2;
  /// let n2 = get_or_create(depth);
  /// /// NPC
  /// assert_eq!(47, n2.from_ring(2));
  /// assert_eq!(29, n2.from_ring(7));
  /// assert_eq!(60, n2.from_ring(22));
  /// // EQR
  /// assert_eq!(51,  n2.from_ring(54));
  /// assert_eq!(88,  n2.from_ring(107));
  /// assert_eq!(174, n2.from_ring(129));
  /// // SPC
  /// assert_eq!(177, n2.from_ring(187));
  /// assert_eq!(153, n2.from_ring(157));
  /// assert_eq!(144, n2.from_ring(189));
  /// ```
  pub fn from_ring(&self, hash: u64) -> u64 {
    // 4*sum from i=1 to nside of i =  4 * nside(nside+1)/2 = 2*nside*(nside+1)
    let first_hash_in_eqr = self.first_hash_in_eqr();
    let first_hash_on_eqr_spc_transition = self.n_hash - first_hash_in_eqr; // == first_hash_on_eqr_spc_transition
    if hash < first_hash_in_eqr { // North polar cap
      // Solve 2*n(n+1) = x 
      //   => 2n^2+2n-x = 0 => b^2-4ac = 4+8x = 4(1+2x)
      //   => n = [-2+2*sqrt(1+2x)]/4 => n = [sqrt(1+2x) - 1] / 2
      //   => n^2+n-x/2 = 0 => b^2-4ac = 1 + 2x => n = [sqrt(1+2x) - 1] / 2
      // n - 1 = ring index
      // Here we may optimize by finding a good 'isqrt' implementation
      let i_ring: u64 = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
      let n_in_ring: u64 = i_ring + 1;
      let i_in_ring = hash - triangular_number_x4(i_ring);
      let d0h = i_in_ring / n_in_ring;
      let h = (((self.nside as u64) << 1) - 2) as i64 - i_ring as i64;
      let l = ((i_in_ring - n_in_ring * d0h) << 1) as i64 - i_ring as i64;
      self.build_hash_from_parts (
        d0h as u8,
        ((h + l) >> 1) as u32,
        ((h - l) >> 1) as u32,
      )
    } else if hash >= first_hash_on_eqr_spc_transition { // South polar cap
      let hash = self.n_hash - 1 - hash; // start counting in reverse order from south polar cap
      let i_ring = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
      let n_in_ring = i_ring + 1;
      let i_in_ring = ((n_in_ring << 2) - 1) - (hash - triangular_number_x4(i_ring));
      let d0h = i_in_ring / n_in_ring;
      let h = i_ring as i64;
      let l = ((i_in_ring - n_in_ring * d0h) << 1) as i64 - i_ring as i64;
      self.build_hash_from_parts (
        d0h as u8 + 8,
        ((h + l) >> 1) as u32,
        ((h - l) >> 1) as u32,
      )
    } else { // Equatorial region
      // Set origin of ring indexes at the center of small cell in north corner of base cell 4 (North to South direction)
      let mut i_ring = hash - first_hash_in_eqr;
      let mut i_in_ring = i_ring;
      // <=> /= 4*nside (number of hash per line) => count the number of line from first equatorial line
      i_ring >>= self.depth + 2;
      // Substract number of hash in previous rings (-= n_rings * 4*nside)
      i_in_ring -= i_ring << (self.depth + 2);
      let l = (i_in_ring << 1) + div2_remainder(i_ring);
      // Set origin of h axis at center of small cell in south corner of base cell 4 (South to North direction)
      let h = (((self.nside as u64) << 1) - 2) - i_ring;
      // Rotation of -45
      let i_in_d0c = (h + l) >> 1;
      let j_in_d0c = (h as i64 - l as i64) >> 1;
      // Offset of 4*nside in j
      let j_in_d0c = (j_in_d0c + ((self.nside as i64) << 2)) as u64;
      let i_d0c = self.div_by_nside_floor_u8(i_in_d0c);
      let j_d0c = self.div_by_nside_floor_u8(j_in_d0c);
      self.build_hash_from_parts (
        depth0_hash_unsafe(i_d0c, j_d0c),
        self.modulo_nside(i_in_d0c) as u32,
        self.modulo_nside(j_in_d0c) as u32,
      )
    }
  }
  
  /// arg * nside
  #[inline]
  fn nside_time(&self, i: u64) -> u64 {
    i << self.depth
  }
  
  /// See the same method (generalized to any NSIDE) in the RING module.
  /// Here we addiionally use the fact that NSIDE is a powrer of two.
  #[inline]
  fn first_hash_in_eqr(&self) -> u64 {
    //   2*nside*(nside + 1)
    // = 2*[nside^2 + nside]
    ((1_u64 << (self.depth << 1)) + self.nside as u64) << 1
  }

  /// (i_ring - nside) * 4 * nside
  #[inline]
  fn minus_nside_x_4nside(&self, i_ring: u64) -> u64 {
    (i_ring - self.nside as u64) << (self.depth + 2)
  }
  
  /// Compute the position on the unit sphere of the center (in the Euclidean projection plane)
  /// of the cell associated to the given hash value.
  /// 
  /// # Input
  /// - `hash`: the hash value of the cell we look for the unprojected center
  /// 
  /// # Output
  /// - `(lon, lat)` in radians, the unprojected position (on the unit sphere) of the center of 
  ///   the cell in the Euclidean plane
  ///   - `lon`, longitude in `[0, 2pi]` radians;
  ///   - `lat`, latitude in `[-pi/2, pi/2]` radians. 
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  /// 
  /// # Example
  /// ```rust
  /// use std::f64::consts::{PI};
  /// use cdshealpix::{TRANSITION_LATITUDE};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  /// 
  /// assert!(dist((PI / 4f64, TRANSITION_LATITUDE) , nested0.center(0u64)) < 1e-15);
  /// ```
  ///
  #[inline]
  pub fn center(&self, hash: u64) -> (f64, f64) {
    let (x, y) = self.center_of_projected_cell(hash);
    super::unproj(x, y)
  }

  /// Compute the position on the unit sphere of the position '(dx, dy)' from the south vertex of 
  /// the HEALPix cell associated to the given hash value.
  /// The x-axis is the South-East axis while the y-axis is the south-west axis.
  /// 
  /// # Input
  /// - `hash`: the hash value of the cell in which are defined `dx` and `dy`
  /// - `dx`: the positional offset $\in [0, 1[$ along the south-to-east axis
  /// - `dy`: the positional offset $\in [0, 1[$ along the south-to-west axis
  /// 
  /// # Output
  /// - `(lon, lat)` in radians, the unprojected position (on the unit sphere) of the given position 
  ///   inside the given cell in the Euclidean plane
  ///   - `lon`, longitude in `[0, 2pi]` radians;
  ///   - `lat`, latitude in `[-pi/2, pi/2]` radians. 
  /// 
  /// # Panics
  /// This method panics if either:
  /// - the given `hash` value is not in `[0, 12*nside^2[`, 
  /// - `dx` or `dy` is not $\in [0, 1[$
  /// 
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  /// 
  /// assert!(dist(nested0.sph_coo(0, 0.5, 0.5) , nested0.center(0)) < 1e-15);
  /// ```
  ///
  pub fn sph_coo(&self, hash: u64, dx: f64, dy: f64) -> (f64, f64) {
    assert!(0.0 <= dx && dx < 1.0);
    assert!(0.0 <= dy && dy < 1.0);
    let (mut x, mut y) = self.center_of_projected_cell(hash);
    x += (dx - dy) * self.one_over_nside;
    y += (dx + dy - 1.0) * self.one_over_nside;
    super::unproj(ensures_x_is_positive(x), y)
  }
  
  /// Computes the position on the unit sphere of the cell vertex located at the given *direction*
  /// with respect to the center of the cell.
  ///   
  /// # Input
  /// - `hash`: the hash value of the cell we look for the position of a vertex
  /// - `vertex_direction`: the direction of the wanted vertex coordiantes
  /// 
  /// # Output
  /// - `(lon, lat)` in radians, the position (on the unit sphere) of the vertex
  ///   - `lon`, longitude in `[0, 2pi]` radians;
  ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  ///
  /// # Example
  /// ```rust
  /// use std::f64::consts::{PI};
  /// use cdshealpix::compass_point::{Cardinal};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  ///
  /// assert!(dist((PI / 4f64, 0.0) , nested0.vertex(0, Cardinal::S)) < 1e-15);
  /// ```
  #[inline]
  pub fn vertex(&self, hash: u64, vertex_direction: Cardinal) -> (f64, f64) {
    let (x, y) = self.center_of_projected_cell(hash);
    self.vertex_lonlat(x, y, &vertex_direction)
  }

  /// Computes the positions on the unit sphere of the 4 vertices of the given cell.
  /// If you want to access the position for a given direction, use method 
  /// [vertices_map](#method.vertices_map).
  ///   
  /// # Input
  /// - `hash`: the hash value of the cell we look for the positions of its vertices
  /// 
  /// # Output
  /// - `[(lon_S, lat_S), (lon_E, lat_E), (lon_N, lat_N), (lon_W, lat_W)]` in radians, 
  ///   the positions (on the unit sphere) of the vertices
  ///   - `lon`, longitude in `[0, 2pi]` radians;
  ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  ///
  /// # Example
  /// ```rust
  /// use std::f64::consts::{PI};
  /// use cdshealpix::compass_point::{Cardinal};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  ///
  /// assert!(dist((PI / 4f64, 0.0) , nested0.vertices(0)[0]) < 1e-15);
  /// ```
  #[inline]
  pub fn vertices(&self, hash: u64) -> [(f64, f64); 4] {
    let (x, y) = self.center_of_projected_cell(hash);
    [
      super::unproj(x, y - self.one_over_nside), // S
      super::unproj(x + self.one_over_nside, y), // E
      super::unproj(x, y + self.one_over_nside), // N
      super::unproj(ensures_x_is_positive(x - self.one_over_nside), y)  // W
    ]
  }

  /// Computes the positions on the unit sphere of the vertices of the given cell which direction 
  /// are in the given set.
  /// If you don't care about the association between position and direction, 
  /// you should use [vertices](#method.vertices).
  ///   
  /// # Input
  /// - `hash`: the hash value of the cell we look for the positions of its vertices
  /// 
  /// # Output
  /// - the vertices position stored in a map associating each vertex direction with its position.
  ///   - `lon`, longitude in `[0, 2pi]` radians;
  ///   - `lat`, latitude in `[-pi/2, pi/2]` radians.
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  ///
  /// # Example
  /// ```rust
  /// use std::f64::consts::{PI};
  /// use cdshealpix::compass_point::{Cardinal, CardinalSet};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  ///
  /// assert!(dist((PI / 4f64, 0.0) , *nested0.vertices_map(0, CardinalSet::all()).get(Cardinal::S).unwrap()) < 1e-15);
  /// ```
  #[inline]
  pub fn vertices_map(&self, hash: u64, directions: CardinalSet) -> CardinalMap<(f64, f64)> {
    let (x, y) = self.center_of_projected_cell(hash);
    let mut result_map = CardinalMap::new();
    for direction in directions {
      let vertex = self.vertex_lonlat(x, y, &direction);
      result_map.put(direction, vertex);
    }
    result_map
  }
  
  /// Computes a list of positions on a given side of a given HEALPix cell on the unit sphere.
  /// 
  /// # Input
  /// - `hash`: the hash value of the cell we look for side path on the unit sphere.
  /// - `from_vertex`: direction (from the cell center) of the path starting vertex
  /// - `to_vertex`: direction (from the cell center) of the path ending vertex
  /// - `include_to_vertex`: if set to *false*, the result contains `n_segments` points and do
  ///                        not include the ending vertex.
  ///                        Else the result contains `n_segments + 1` points.
  /// - `n_segments`: number of segments in the path from the starting vertex to the ending vertex
  ///
  /// # Output 
  /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
  ///
  pub fn path_along_cell_side(&self, hash: u64, from_vertex:  &Cardinal, to_vertex: &Cardinal,
    include_to_vertex: bool, n_segments: u32) -> Box<[(f64, f64)]> {
    let n_points: usize = if include_to_vertex { n_segments + 1 } else { n_segments } as usize;
    let mut path_points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
    let proj_center = self.center_of_projected_cell(hash);
    self.path_along_cell_side_internal(proj_center, from_vertex, to_vertex, include_to_vertex, n_segments, &mut path_points);
    path_points.into_boxed_slice()
  }
  
  fn path_along_cell_side_internal(&self, proj_center: (f64, f64), from_vertex: &Cardinal, to_vertex: &Cardinal,
                                   include_to_vertex: bool, n_segments: u32, path_points: &mut Vec<(f64, f64)>) {
    let n_points: usize = if include_to_vertex { n_segments + 1 } else { n_segments } as usize;
    // Compute starting point offsets
    let from_offset_x = (*from_vertex).offset_we(self.one_over_nside);
    let from_offset_y = (*from_vertex).offset_sn(self.one_over_nside);
    // Compute stepX and stepY
    let step_x = ((*to_vertex).offset_we(self.one_over_nside) - from_offset_x) / (n_segments as f64);
    let step_y = ((*to_vertex).offset_sn(self.one_over_nside) - from_offset_y) / (n_segments as f64);
    // Compute intermediary vertices
    for i in 0..n_points {
      let k = i as f64;
      let x = proj_center.0 + from_offset_x + k * step_x;
      let y = proj_center.1 + from_offset_y + k * step_y;
      path_points.push(super::unproj(ensures_x_is_positive(x), y));
    }
  }

  /// Computes a list of positions on the edge of a given HEALPix cell on the unit sphere.
  /// 
  /// # Input
  /// - `hash`: the hash value of the cell we look for side path on the unit sphere.
  /// - `starting_vertex`: direction (from the cell center) of the path starting vertex
  /// - `clockwise_direction`: tells if the path is in the clockwise or anti-clockwise direction
  /// - `n_segments_by_side`: number of segments in each each side. Hence, the total number of 
  ///                         points in the path equals *4 x n_segments_by_side*.
  /// # Output 
  /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
  ///
  pub fn path_along_cell_edge(&self, hash: u64, starting_vertex: &Cardinal, clockwise_direction: bool,
                              n_segments_by_side: u32) -> Box<[(f64, f64)]> {
    // Prepare space for the result
    let mut path_points: Vec<(f64, f64)> = Vec::with_capacity((n_segments_by_side << 2) as usize);
    // Compute center
    let proj_center = self.center_of_projected_cell(hash);
    // Unrolled loop over successive sides
    // - compute vertex sequence
    let (v1, v2, v3, v4) =  if clockwise_direction {
      starting_vertex.clockwise_cycle()
    } else {
      starting_vertex.counter_clockwise_cycle()
    };
    // - make the five sides
    self.path_along_cell_side_internal(proj_center, &v1, &v2, false, n_segments_by_side, &mut path_points);
    self.path_along_cell_side_internal(proj_center, &v2, &v3, false, n_segments_by_side, &mut path_points);
    self.path_along_cell_side_internal(proj_center, &v3, &v4, false, n_segments_by_side, &mut path_points);
    self.path_along_cell_side_internal(proj_center, &v4, &v1, false, n_segments_by_side, &mut path_points);
    path_points.into_boxed_slice()
  }
  
  
  /// Computes the positions on the sky of each points located on a regular grid in the projection
  /// plane. The grid x-axis is the South-to-east axis and the y-axis is the south-to-west axis.
  /// The return array contains the square fo n_segments_by_side + 1 elements.
  /// 
  /// # Input
  /// - `hash`: the hash value of the cell we look for the grid on the unit sphere.
  /// - `n_segments_by_side`: number of segments in each each side. Hence, the total number of 
  ///                         points in the path equals *(n_segments_by_side + 1)^2*.
  /// # Output 
  /// - the list of positions on the given side of the given HEALPix cell on the unit sphere.
  ///
  /// # Motivation
  /// - to create a mesh in Unity
  pub fn grid(&self, hash: u64, n_segments_by_side: u16) -> Box<[(f64, f64)]> {
    let n_points_per_side = (n_segments_by_side as usize) + 1;
    // Prepare space for the result
    let mut grid: Vec<(f64, f64)> = Vec::with_capacity(n_points_per_side * n_points_per_side);
    // Compute center
    let proj_center = self.center_of_projected_cell(hash);
    // Compute grid
    for i in 0..n_points_per_side {
      let x = (i as f64) / (n_segments_by_side as f64);   // in [0, 1]
      for j in 0..n_points_per_side {
        let y = (j as f64) / (n_segments_by_side as f64); // in [0, 1]
        let l = x - y;
        let h = x + y - 1.0;
        grid.push(super::unproj(proj_center.0 + l * self.one_over_nside, proj_center.1 + h * self.one_over_nside));
      }
    }
    grid.into_boxed_slice()
  }
  
  //////////////////////////
  // NEIGHBOURS Functions //
  //////////////////////////

  /// Retuns the hash value of the neighbour cell of the cell of given hash, in the given direction.
  /// If the cell do not have a neighbour in the given direction (which is the case of the
  /// eastmost and westmost cells in polar caps base cells and northmost and southmost cells of the
  /// equatorial region base cells), the return Option is None.
  /// 
  /// # Input
  /// - `hash` the hash value of the cell we look for the neighbour
  /// - `direction` the direction of the neighbour we look for the cell number
  /// 
  /// # Output
  /// - the cell number (hash value) of the neighbour of the given cell in the given direction
  ///   (`None` if their is no neighbour in the given direction) .
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  ///
  /// assert_eq!(5 , nested0.neighbour(4, MainWind::E).unwrap());
  /// ```
  #[inline]
  pub fn neighbour(&self, hash: u64, direction: MainWind) -> Option<u64> {
    let h_parts: HashParts = self.decode_hash(hash);
    self.neighbour_from_parts(h_parts.d0h, h_parts.i, h_parts.j, direction)
  }

  /// Returns the hash values of all the neighbour cells of the cell of given hash.
  /// The given cell itself can be included (setting the `include_center` parameters to `true`).
  /// 
  /// # Input
  /// - `hash` the hash value of the cell we look for the neighbours
  /// - `include_center` include (or not) the input cell in the MainWind::C key of the returned map
  /// 
  /// # Output
  /// - the cell number (hash value) of the neighbour of the given cell in the given direction
  ///   (`None` if their is no neighbour in the given direction) .
  /// 
  /// # Panics
  /// If the given `hash` value is not in `[0, 12*nside^2[`, this method panics.
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 0u8;
  /// let nested0 = get_or_create(depth);
  ///
  /// assert_eq!(5 , *nested0.neighbours(4, false).get(MainWind::E).unwrap());
  /// ```
  pub fn neighbours(&self, hash: u64, include_center: bool) -> MainWindMap<u64> {
    self.check_hash(hash);
    let mut result_map = MainWindMap::new();
    if include_center {
      result_map.put(C, hash);
    }
    let h_bits: HashBits = self.pull_bits_appart(hash);
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      self.edge_cell_neighbours(hash, &mut result_map);
    } else {
      self.inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, &mut result_map);
    }
    result_map
  }
  
  /// Returns the hash values corresponding to the internal bounds of the given `hash` at 
  /// the hash depth + the given `delta_depth`:
  /// - the first quarter contains the southeast border (the z-order curve x-axis with y = 0); 
  /// - the second quarter contains the northeast border (the z-order y-axis with x = xmax - 1);
  /// - the third quarter contains the northwest border (the z-order curve x-axis with y = ymax - 1); 
  /// - the forth quarter contains the southwest border (the y-axis with x = 0).
  /// 
  /// The hashes are ordered consecutively, starting from the south (x=0, y=0) cell in the 
  /// anti-clokwise direction.
  /// 
  /// # Input
  /// - `hash ` the hash for which we look for the internal bounds
  /// - `delta_depth` difference between the depth of the edge cells and the depth of the given cell  
  ///
  /// # Output
  /// - the cell numbers (hash values) of the given hash inner edge ordered consecutively, 
  /// starting from the south (x=0, y=0) cell in the anti-clokwise direction.
  ///
  pub fn internal_edge(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
    // Compute the x and y part masks for deltaDepth.
    let zoc = get_zoc(delta_depth);
    let twice_dd = delta_depth << 1;
    let x_max_bits = x_mask(delta_depth);
    let y_max_bits = x_max_bits << 1;
    // Prepare hashes of depth of (self.depth + delta_depth), switching hash bits of 2 delta_depth to the left.
    hash <<= twice_dd;
    // Prepare filling the result.
    // am1 stands for a - 1, i.e. nSide - 1, i.e. the index of the last cell along the x or y-axis
    let am1 = (1_u32 << delta_depth) - 1; // 2^deltaDepth - 1
    let mut res: Vec<u64> = Vec::with_capacity((am1 << 2) as usize);
    // Southeast axis
    res.push(hash);
    for k in 1..am1 {
      let x_bits = zoc.i02h(k);
      res.push(hash | x_bits);
    }
    // Northeast axis
    res.push(hash | x_max_bits);
    for k in 1..am1 {
      res.push(hash | zoc.oj2h(k) | x_max_bits);
    }
    // Northwest axis
    res.push(hash | y_max_bits | x_max_bits);
    for k in 1..am1 {
      res.push(hash | y_max_bits | zoc.i02h(am1 - k));
    }
    // Southwest axis
    res.push(hash | y_max_bits);
    for k in 1..am1 {
      res.push(hash | zoc.oj2h(am1 - k));
    }
    res.into_boxed_slice()
  }

  /// Same as method [internal_edge](#method.internal_edge) except that the returned array is sorted.
  pub fn internal_edge_sorted(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
    // Compute the x and y part masks for deltaDepth.
    let zoc = get_zoc(delta_depth);
    let twice_dd = delta_depth << 1;
    let x_max_bits = x_mask(delta_depth);
    let y_max_bits = x_max_bits << 1;
    // Prepare hashes of depth of (self.depth + delta_depth), switching hash bits of 2 delta_depth to the left.
    hash <<= twice_dd;
    // Set grid size (nSide inside the cell of depth this.depth)
    let nside = 1 << delta_depth;
    let am1 = nside - 1;
    let n_half_side = nside >> 1;
    // South sub-square (dividing in 4 sub-squares)
    let mut x = 1_u32;
    let mut lim = 2_u32;
    let mut k0 = 1_usize; 
    let mut k1 = 2_usize;
    let mut k2 = (am1 + n_half_side) as usize;
    let mut k3 = ((am1 << 1) + n_half_side) as usize;
    let size = (am1 << 2) as usize;
    let mut result = vec![0_u64; size];
    // Set South corner (first element)
    result[0] = hash;
    // Set east corner
    result[k2 - 1] = hash | x_max_bits;
    // Set west corner 
    result[k3 - 1] = hash | y_max_bits;
    // Set north corner (last element)
    result[size - k0] = hash | y_max_bits | x_max_bits;
    while x < n_half_side { // while (k < nHalfSize)
      x += 1;
      let xn = zoc.ij2h(x, nside - x); // x shuffled
      let xs = xn & x_max_bits;
      let xn = xn & y_max_bits;
      // South square, south east part
      result[k0] = hash | xs;
      k0 += 1;
      // South square, south west part
      result[k1] = hash | (xs << 1);
      k1 += 1;
      // East square, north east part
      result[k2] = hash | (xs << 1) | x_max_bits;
      k2 += 1;
      // West square, north west part
      result[k3] = hash | y_max_bits | xs;
      k3 += 1;
      // North square, north west
      result[size - k0] = hash | y_max_bits | (xn >> 1);
      // North square, north east
      result[size - k1] = hash | xn | x_max_bits;
      // West square, north west part
      result[size - k2] = hash | xn;
      // East square, nort east part
      result[size - k3] = hash | (xn >> 1);
      // Change k0, k1 and limit if x== limit.
      // The following lines of code are equivalent to:
      /* if (x == lim) {
              k0 = k1;
              k1 += lim; // +2 +4 +8
              lim <<= 1; // 4 8 16 32 ...
          } */
      // To be tested if they are faster (since no risk of branch miss-prediction):
      // probably true for small deltaDepth but not for large deltaDepth.
      let mut tmp = x & lim;  debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == x));
      k0 += (tmp >> 1) as usize;
      k1 += tmp as usize;
      tmp -= x;           debug_assert!((x < lim            ) || (x == lim && tmp == 0));
      tmp = 1 >> tmp;     debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == 1));
      lim <<= tmp;
    }
    result.into_boxed_slice()
  }
  
  /// Similar to [external_edge](#method.external_edge) except that the returned structure allow
  /// to access each elements of the external edge: the 4 corners plus the 4 edges.
  pub fn external_edge_struct(&self, hash: u64, delta_depth: u8) -> ExternalEdge {
    self.check_hash(hash);
    let mut res = ExternalEdge::new_empty();
    let h_bits: HashBits = self.pull_bits_appart(hash);
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      // Not easy: opposite directions depends on base cell neighbours
      let mut neighbours = MainWindMap::new();
      self.edge_cell_neighbours(hash, &mut neighbours);
      let h_parts: HashParts = self.decode_hash(hash);
      for (direction, hash_value) in neighbours.entries_vec().drain(..) {
        
       let dir_from_neig = if h_parts.d0h == self.h_2_d0h(hash_value) {
          direction.opposite()
        } else if self.depth == 0 {
	  direction_from_neighbour(h_parts.d0h, &direction)
        } else {
          let dir_in_basce_cell_border = self.direction_in_base_cell_border(h_bits.i, h_bits.j);
          // println!("B: {:?}, {}, {:?}", &direction, &hash_value, &dir_in_basce_cell_border);
          edge_cell_direction_from_neighbour(h_parts.d0h, &dir_in_basce_cell_border, &direction)
        };
        // println!("{:?}, {}, {:?}", &direction, &hash_value, &dir_from_neig);
        add_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, &direction,&mut res);
      }
    } else {
      // Easy: always use the opposite direction
      let mut neighbours = MainWindMap::new();
      self.inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, &mut neighbours);
      for (direction, hash_value) in neighbours.entries_vec().drain(..) {
        let dir_from_neig = direction.opposite();
        add_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, &direction,&mut res);
      }
    }
    res
  }
  
  /// Provides the list of all cells of depth this layer depth + the given `delta_depth`
  /// surrounding the cell of given hash value.  
  /// 
  /// Here the result of both following codes:
  /// ![External edge depth 1, cells 10 and 11, delta_depth = +2](external_edge.png)
  /// 
  /// ```rust
  /// use cdshealpix::nested::{external_edge_sorted};
  /// 
  /// let depth = 1;
  /// let delta_depth = 2;
  /// 
  /// let hash = 10;
  /// let actual_res = external_edge_sorted(depth, hash, delta_depth);
  /// let expected_res: [u64; 19] = [85, 87, 93, 95, 117, 138, 139, 142, 143, 154, 176, 178, 184, 186, 415, 437, 439, 445, 447];
  /// for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
  ///   assert_eq!(h1, h2);
  /// }
  /// assert_eq!(expected_res.len(), actual_res.len());
  /// 
  /// let hash = 11;
  /// let actual_res = external_edge_sorted(depth, hash, delta_depth);
  /// let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
  /// for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
  ///   assert_eq!(h1, h2);
  /// }
  /// assert_eq!(expected_res.len(), actual_res.len());
  /// 
  /// ```
  pub fn external_edge(&self, hash: u64, delta_depth: u8) -> Box<[u64]> {
    self.external_edge_generic(hash, delta_depth, false)
  }
  
  /// Similar to [external_edge](#method.external_edge) except that the returned list of cells is ordered.
  pub fn external_edge_sorted(&self, hash: u64, delta_depth: u8) -> Box<[u64]> {
    self.external_edge_generic(hash, delta_depth, true)
  }

  fn external_edge_generic(&self, hash: u64, delta_depth: u8, sorted: bool) -> Box<[u64]> {
    self.check_hash(hash);
    let mut edge = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
    let h_bits: HashBits = self.pull_bits_appart(hash);
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      // Not easy: opposite directions depends on base cell neighbours
      let mut neighbours = MainWindMap::new();
      self.edge_cell_neighbours(hash, &mut neighbours);
      let mut neighbours = if sorted { neighbours.sorted_entries_vec() } else { neighbours.entries_vec() };
      let h_parts: HashParts = self.decode_hash(hash);
      for (direction, hash_value) in neighbours.drain(..) {
        let dir_from_neig = if h_parts.d0h == self.h_2_d0h(hash_value) {
          direction.opposite()
        } else if self.depth == 0 {
          direction_from_neighbour(h_parts.d0h, &direction)
        }  else {
          edge_cell_direction_from_neighbour(h_parts.d0h, &self.direction_in_base_cell_border(h_bits.i, h_bits.j), &direction)
        };
        append_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, &mut edge);
      }
    } else {
      // Easy: always use the opposite direction
      let mut neighbours = MainWindMap::new();
      self.inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, &mut neighbours);
      let mut neighbours = if sorted { neighbours.sorted_entries_vec() } else { neighbours.entries_vec() };
      for (direction, hash_value) in neighbours.drain(..) {
        append_sorted_internal_edge_element(hash_value, delta_depth,  direction.opposite(), &mut edge);
      }
    }
    edge.into_boxed_slice()
  }
  
  #[inline]
  fn is_in_base_cell_border(&self, i_in_base_cell_bits: u64, j_in_base_cell_bits: u64) -> bool {
    0_u64 == i_in_base_cell_bits || i_in_base_cell_bits == self.x_mask
      || 0_u64 == j_in_base_cell_bits || j_in_base_cell_bits == self.y_mask
  }

  #[inline]
  fn direction_in_base_cell_border(&self, i_in_base_cell_bits: u64, j_in_base_cell_bits: u64) -> MainWind {
    let i = if 0_u64 == i_in_base_cell_bits {
      0
    } else if i_in_base_cell_bits == self.x_mask {
      2
    } else {
      1
    };
    let j = if 0_u64 == j_in_base_cell_bits {
      0
    } else if j_in_base_cell_bits == self.y_mask {
      2
    } else {
      1
    };
    MainWind::from_index(3 * j + i)
  }

  fn inner_cell_neighbours(&self, d0h_bits: u64, i_in_d0h_bits: u64, j_in_d0h_bits: u64,
                           result_map: &mut MainWindMap<u64>) {
    let ij = self.z_order_curve.h2ij(i_in_d0h_bits | j_in_d0h_bits);
    let i = self.z_order_curve.ij2i(ij);
    let j = self.z_order_curve.ij2j(ij);
    // Compute i-1 and j-1 bits.
    // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
    // in making fewer operations
    let ij = self.z_order_curve.ij2h(i - 1, j - 1);
    let im1_bits = ij & self.x_mask;
    let jm1_bits = ij & self.y_mask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    let ij = self.z_order_curve.ij2h(i + 1, j + 1);
    let ip1_bits = ij & self.x_mask;
    let jp1_bits = ij & self.y_mask;
    // Unrolled for loop an MainWind enumset
    result_map.put(S, bits_2_hash(d0h_bits, im1_bits, jm1_bits));
    result_map.put(SE, bits_2_hash(d0h_bits, i_in_d0h_bits, jm1_bits));
    result_map.put(E, bits_2_hash(d0h_bits, ip1_bits, jm1_bits));
    result_map.put(SW, bits_2_hash(d0h_bits, im1_bits, j_in_d0h_bits));
    result_map.put(NE, bits_2_hash(d0h_bits, ip1_bits, j_in_d0h_bits));
    result_map.put(W, bits_2_hash(d0h_bits, im1_bits, jp1_bits));
    result_map.put(NW, bits_2_hash(d0h_bits, i_in_d0h_bits, jp1_bits));
    result_map.put(N, bits_2_hash(d0h_bits, ip1_bits, jp1_bits));
  }

  fn edge_cell_neighbours(&self, hash: u64, result_map: &mut MainWindMap<u64>) {
    // Could have simply been edgeCellNeighbours(hash, EnumSet.allOf(MainWind.class) result)
    // but we prefered to unroll the for loop.
    let h_parts: HashParts = self.decode_hash(hash);
    let d0h = h_parts.d0h;
    let i = h_parts.i;
    let j = h_parts.j;
    result_map.put_opt(S, self.neighbour_from_parts(d0h, i, j, S));
    result_map.put_opt(SE, self.neighbour_from_parts(d0h, i, j, SE));
    result_map.put_opt(E, self.neighbour_from_parts(d0h, i, j, E));
    result_map.put_opt(SW, self.neighbour_from_parts(d0h, i, j, SW));
    result_map.put_opt(NE, self.neighbour_from_parts(d0h, i, j, NE));
    result_map.put_opt(W, self.neighbour_from_parts(d0h, i, j, W));
    result_map.put_opt(NW, self.neighbour_from_parts(d0h, i, j, NW));
    result_map.put_opt(N, self.neighbour_from_parts(d0h, i, j, N));
  }

  fn neighbour_from_parts(&self, d0h: u8, i: u32, j: u32, dir: MainWind) -> Option<u64> {
    let i = (i as i32) + (dir.offset_se() as i32);
    let j = (j as i32) + (dir.offset_sw() as i32);
    let d0_neighbour_dir = MainWind::from_offsets(
      self.neighbour_base_cell_offset(i),
      self.neighbour_base_cell_offset(j));
    self.neighbour_from_shifted_coos(d0h, i as u32, j as u32, d0_neighbour_dir)
  }

  /*fn neighbour_base_cell_offset(&self, offset: i8, coo: u32) -> i8 {
    if coo == 0_u32 && offset == -1_i8 {
      -1_i8
    } else if coo == self.nside_minus_1 && offset == 1_i8 {
      1_i8
    } else {
      0_i8
    }
  }*/
  /// This method has a single input parameters `coo` which must be in `[-1, nside]`, and returns:
  /// - -1 if `coo` == -1
  /// -  0 if `coo` in `[0, nside[`
  /// -  1 if `coo` == nside
  #[inline]
  fn neighbour_base_cell_offset(&self, coo_in_base_cell: i32) -> i8 {
    debug_assert!(-1_i32 <= coo_in_base_cell && coo_in_base_cell <= (self.nside as i32));
    let offset = (coo_in_base_cell >> 31 | coo_in_base_cell >> self.depth) as i8;
    debug_assert!(
      (coo_in_base_cell == -1_i32 && offset == -1_i8)
        || (coo_in_base_cell == (self.nside as i32) && offset == 1_i8)
        || (-1_i32 < coo_in_base_cell && coo_in_base_cell < (self.nside as i32) && offset == 0_i8)
    );
    offset
  }

  #[inline]
  fn neighbour_from_shifted_coos(&self, d0h: u8, i: u32, j: u32, base_cell_neighbour_dir: MainWind) -> Option<u64> {
    if base_cell_neighbour_dir == MainWind::C {
      debug_assert!(i < self.nside && j < self.nside);
      self.build_hash_from_parts_opt(d0h, i, j)
    } else {
      let d0h_mod_4 = d0h & 3_u8;  // <=> base_cell modulo 4
      match div4_quotient(d0h) {
        // <=> base_cell / 4
        0 => self.ncp_neighbour(d0h_mod_4, i, j, base_cell_neighbour_dir),
        1 => self.eqr_neighbour(d0h_mod_4, i, j, base_cell_neighbour_dir),
        2 => self.spc_neighbour(d0h_mod_4, i, j, base_cell_neighbour_dir),
        _ => panic!("Base cell must be in [0, 12["),
      }
    }
  }

  #[inline]
  fn ncp_neighbour(&self, d0h_mod_4: u8, i: u32, j: u32, base_cell_neighbour_dir: MainWind) -> Option<u64> {
    let m = self.nside_minus_1;
    match base_cell_neighbour_dir {
      S => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 2), m, m),
      SE => self.build_hash_from_parts_opt(base_cell(next(d0h_mod_4), 1), i, m),
      SW => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 1), m, j),
      NE => self.build_hash_from_parts_opt(base_cell(next(d0h_mod_4), 0), j, m),
      NW => self.build_hash_from_parts_opt(base_cell(prev(d0h_mod_4), 0), m, i),
      N => self.build_hash_from_parts_opt(base_cell(oppo(d0h_mod_4), 0), m, m),
      _ => None,
    }
  }

  #[inline]
  fn eqr_neighbour(&self, d0h_mod_4: u8, i: u32, j: u32, base_cell_neighbour_dir: MainWind) -> Option<u64> {
    let m = self.nside_minus_1;
    match base_cell_neighbour_dir {
      SE => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 2), i, m),
      E => self.build_hash_from_parts_opt(base_cell(next(d0h_mod_4), 1), 0, m),
      SW => self.build_hash_from_parts_opt(base_cell(prev(d0h_mod_4), 2), m, j),
      NE => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 0), 0, j),
      W => self.build_hash_from_parts_opt(base_cell(prev(d0h_mod_4), 1), m, 0),
      NW => self.build_hash_from_parts_opt(base_cell(prev(d0h_mod_4), 0), i, 0),
      _ => None,
    }
  }

  #[inline]
  fn spc_neighbour(&self, d0h_mod_4: u8, i: u32, j: u32, base_cell_neighbour_dir: MainWind) -> Option<u64> {
    match base_cell_neighbour_dir {
      S => self.build_hash_from_parts_opt(base_cell(oppo(d0h_mod_4), 2), 0, 0),
      SE => self.build_hash_from_parts_opt(base_cell(next(d0h_mod_4), 2), 0, i),
      SW => self.build_hash_from_parts_opt(base_cell(prev(d0h_mod_4), 2), j, 0),
      NE => self.build_hash_from_parts_opt(base_cell(next(d0h_mod_4), 1), 0, j),
      NW => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 1), i, 0),
      N => self.build_hash_from_parts_opt(base_cell(iden(d0h_mod_4), 0), 0, 0),
      _ => None,
    }
  }
  
  /// Center of the given cell in the Euclidean projection space.
  /// # Output
  /// - `(x, y)` coordinates such that $x \in [0, 8[$ and $y \in [-2, 2]$. 
  pub fn center_of_projected_cell(&self, hash: u64) -> (f64, f64) {
    self.check_hash(hash);
    let h_parts: HashParts = self.decode_hash(hash);
    let mut hl: (i32, i32) = rotate45_scale2(h_parts.i, h_parts.j);
    self.shift_from_small_cell_center_to_base_cell_center(&mut hl);
    let mut xy: (f64, f64) = self.scale_to_proj_dividing_by_nside(hl);
    let (offset_x, offset_y) = compute_base_cell_center_offsets_in_8x3_grid(h_parts.d0h);
    apply_base_cell_center_offsets(&mut xy, offset_x, offset_y);
    xy
  }

  // Computes the position on the unit sphere of the vertex, located at the given direction,
  // of the cell of given center coordinate on the projection plane.
  #[inline]
  fn vertex_lonlat(&self, center_x: f64, center_y: f64, vertex_direction: &Cardinal) -> (f64, f64) {
    let x = center_x + (*vertex_direction).offset_we(self.one_over_nside);
    let y = center_y + (*vertex_direction).offset_sn(self.one_over_nside);
    super::unproj(ensures_x_is_positive(x), y)
  }

  #[inline]
  fn is_hash(&self, hash: u64) -> bool { hash < self.n_hash }

  #[inline]
  fn check_hash(&self, hash: u64) { assert!(self.is_hash(hash), "Wrong hash value: too large."); }

  #[inline]
  fn h_2_d0h(&self, hash: u64) -> u8 {
      (hash >> self.twice_depth) as u8
  }
  
  #[inline]
  fn decode_hash(&self, hash: u64) -> HashParts {
    let ij: u64 = self.z_order_curve.h2ij(hash & self.xy_mask);
    HashParts {
      d0h: self.h_2_d0h(hash),
      i: self.z_order_curve.ij2i(ij),
      j: self.z_order_curve.ij2j(ij),
    }
  }

  #[inline]
  fn pull_bits_appart(&self, hash: u64) -> HashBits {
    HashBits {
      d0h: hash & self.d0h_mask,
      i: hash & self.x_mask,
      j: hash & self.y_mask,
    }
  }

  #[inline]
  fn shift_from_small_cell_center_to_base_cell_center(&self, ij: &mut (i32, i32)) {
    let (ref mut _i, ref mut j) = *ij;
    *j -= self.nside_remainder_mask as i32; // nside_remainder_mask == nside - 1
  }

  #[inline]
  fn scale_to_proj_dividing_by_nside(&self, (x, y): (i32, i32)) -> (f64, f64) {
    (x as f64 * self.one_over_nside, y as f64 * self.one_over_nside)
  }
  
  #[inline]
  fn h_and_shs_to_lower_h(&self, deeper_depth: u8) -> impl Fn((u64, f64)) -> u64 {
    let twice_depth_diff = (deeper_depth - self.depth) << 1;
    move |(h, _)| h >> twice_depth_diff
  }

  ////////////////////////////
  // Bilinear interpolation //
  ////////////////////////////
  
  /// See [wikipeida](https://en.wikipedia.org/wiki/Bilinear_interpolation) about bilinear interpolation.
  /// The main difficulty here are the corners of base cells for which the number of neighbours is not
  /// equals to 8.
  /// In the normal case we have:
  /// ```math
  /// f(x, y) = f(0, 0) (1 - x) (1 - y) 
  ///         + f(1, 0) x (1 - y)         
  ///         + f(0, 1) (1 - x) y
  ///         + f(1, 1) x y
  /// ```
  /// If a neighbour is missing, we share equally its contribution between the 2 cells that do not
  /// contains the given coordinate, and we fill the array with the cell of the given coordinate
  /// with a weight of 0.  
  /// # Output
  /// - `[(cell, weigth), (cell, weigth), (cell, weigth), (cell, weigth)]` the cell number
  ///    together with their weight
  pub fn bilinear_interpolation(&self, lon: f64, lat: f64) -> [(u64, f64); 4] {
    let (h, dx, dy) = self.hash_with_dxdy(lon, lat);
    // We can probably optimize here since we are interested in only 3 neighbours
    let neigbours_map = self.neighbours(h, true);
    // Look at the four pixels
    let xcoo = (dx > 0.5) as u8;
    let ycoo = (dy > 0.5) as u8;
    let quarter: u8 = (ycoo << 1) + xcoo;
    match quarter {
      0 => { // => S => (dx + 0.5, dy + 0.5, S, SE, SW, C)
        match neigbours_map.get(S) {
          Some(nh) => [
            (*nh, (0.5 - dx) * (0.5 - dy)), 
            (*neigbours_map.get(SE).unwrap(), (0.5 + dx) * (0.5 - dy)), 
            (*neigbours_map.get(SW).unwrap(), (0.5 - dx) * (0.5 + dy)), 
            (h, (0.5 + dx) * (0.5 + dy) )
          ],
          None => [
            (h, 0.0),
            (*neigbours_map.get(SE).unwrap(), (0.5 - dy) * (0.75 + 0.5 * dx)),
            (*neigbours_map.get(SW).unwrap(), (0.5 - dx) * (0.75 + 0.5 * dy)),
            (h, (0.5 + dx) * (0.5 + dy))
          ],
        }
      }
      1 => // => E => (dx - 0.5, dy + 0.5, SE, E, C, NE)
        match neigbours_map.get(E) {
          Some(nh) => [
            (*neigbours_map.get(SE).unwrap(), (1.5 - dx) * (0.5 - dy)),
            (*nh, (dx - 0.5) * (0.5 - dy)),
            (h, (1.5 - dx) * (0.5 + dy)),
            (*neigbours_map.get(NE).unwrap(), (dx - 0.5) * (0.5 + dy))
          ],
          None => [
            (*neigbours_map.get(SE).unwrap(), (0.5 - dy) * (1.25 - 0.5 * dx)),
            (h, 0.0),
            (h, (1.5 - dx) * (0.5 + dy)),
            (*neigbours_map.get(NE).unwrap(), (dx - 0.5) * (0.75 + 0.5 * dy))
          ],
        }
      2 => // => W => (dx + 0.5, dy - 0.5, SW, C, W, NW)
        match neigbours_map.get(W) {
          Some(nh) => [
            (*neigbours_map.get(SW).unwrap(), (0.5 - dx) * (1.5 - dy)),
            (h, (dx + 0.5) * (1.5 - dy)),
            (*nh, (0.5 - dx) * (dy - 0.5)),
            (*neigbours_map.get(NW).unwrap(), (0.5 + dx) * (dy - 0.5))
          ],
          None => [
            (*neigbours_map.get(SW).unwrap(), (0.5 - dx) * (1.25 - 0.5 * dy)),
            (h, (dx + 0.5) * (1.5 - dy)),
            (h, 0.0),
            (*neigbours_map.get(NW).unwrap(), (dy - 0.5) * (0.5 * dx + 0.75))
          ],
        }
      3 => // => N => (dx - 0.5, dy - 0.5, C, NE, NW, N)
        match neigbours_map.get(N) {
          Some(nh) => [
            (h, (1.5 - dx) * (1.5 - dy)),
            (*neigbours_map.get(NE).unwrap(), (dx - 0.5) * (1.5 - dy)),
            (*neigbours_map.get(NW).unwrap(), (1.5 - dx) * (dy - 0.5)),
            (*nh, (dx - 0.5) * (dy - 0.5))
          ],
          None => [
            (h, (1.5 - dx) * (1.5 - dy)),
            (*neigbours_map.get(NE).unwrap(), (dx - 0.5) * (1.25 - 0.5 * dy)),
            (*neigbours_map.get(NW).unwrap(), (1.25 - 0.5 * dx) * (dy - 0.5)),
            (h, 0.0)
          ],
        }
      _ => unreachable!(),
    }
  }
  
  //////////////////////
  // Coverage methods //
  //////////////////////
  
  
  /// Returns a hierarchical view of the list of cells overlapped by the given cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
  /// The algorithm is fast but approximated: it may return false positive, 
  /// i.e. cells which are near from the cone but do not overlap it.
  /// To control the approximation, see the method 
  /// [cone_coverage_approx_custom](#method.cone_coverage_approx_custom)
  /// 
  /// # Input
  /// - `cone_lon` the longitude of the center of the cone, in radians
  /// - `cone_lat` the latitude of the center of the cone, in radians
  /// - `cone_radius` the radius of the cone, in radians
  /// 
  /// # Output
  /// - the list of cells overlapped by the given cone, in a BMOC (hierarchical view also telling
  ///   if a cell is fully or partially covered).
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get_or_create(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  /// 
  /// let actual_res = nested3.cone_coverage_approx(lon, lat, radius);
  /// let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///      assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_coverage_approx(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    self.cone_coverage_approx_internal(cone_lon, cone_lat, cone_radius).to_bmoc_packing()
  }
  
  /// Returns a hierarchical view of the list of cells overlapped by the given cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
  /// The algorithm is fast but approximated: it may return false positive, 
  /// i.e. cells which are near from the cone but do not overlap it.
  /// To control the approximation, you can choose to perform the computations at a deeper depth
  /// using the `delta_depth` parameter.
  /// 
  /// # Input
  /// - `delta_depth` the difference between this Layer depth and the depth at which the computations
  ///   are made (should remain quite small).
  /// - `cone_lon` the longitude of the center of the cone, in radians
  /// - `cone_lat` the latitude of the center of the cone, in radians
  /// - `cone_radius` the radius of the cone, in radians
  /// 
  /// # Output
  /// - the list of cells overlapped by the given cone, in a BMOC (hierarchical view also telling
  ///   if a cell is fully or partially covered).
  ///
  /// # Panics
  /// If this layer depth + `delta_depth` > the max depth (i.e. 29)
  /// 
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get_or_create(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  /// 
  /// let actual_res = nested3.cone_coverage_approx_custom(2, lon, lat, radius);
  /// let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_coverage_approx_custom(&self, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    if delta_depth == 0 {
      self.cone_coverage_approx(cone_lon, cone_lat, cone_radius)
    } else {
      // TODO: change the algo not to put all cell in the MOC and pruning it
      get_or_create(self.depth + delta_depth)
        .cone_coverage_approx_internal(cone_lon, cone_lat, cone_radius)
        //.to_lower_depth_bmoc(self.depth)
        .to_lower_depth_bmoc_packing(self.depth)
    }
  }
  
  fn cone_coverage_approx_internal(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOCBuilderUnsafe {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      return self.allsky_bmoc_builder();
    }
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius) {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, distances);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for h in 0..12 {
        self.cone_coverage_approx_recur(0, h,
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                                &minmax_array, 0, &mut bmoc_builder);
      }
      return bmoc_builder; //.to_bmoc_packing();
    }
    // Normal case
    let depth_start = best_starting_depth(cone_radius);
    if depth_start >= self.depth {
      let shs_max = to_squared_half_segment(cone_radius
        + largest_center_to_vertex_distance_with_radius(depth_start, cone_lon, cone_lat, cone_radius));
      let root_layer = get_or_create(depth_start);
      let mut neigs: Vec<u64> = root_layer.neighbours(root_layer.hash(cone_lon, cone_lat), true)
        .values_vec().iter()
        .map(h_to_h_and_shs(cone_lon, cone_lat, cos_cone_lat, root_layer))
        .filter(shs_lower_than(shs_max))
        .map(self.h_and_shs_to_lower_h(depth_start))
        .collect();
      neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
      neigs.dedup();         // remove duplicates (vector must be sorted first)
      // return BMOC::create_unsafe(self.depth, neigs.into_boxed_slice());
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, false); 
      }
      return bmoc_builder;
    } else {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, distances);
      let root_layer = get_or_create(depth_start);
      let root_center_hash = root_layer.hash(cone_lon, cone_lat);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for &root_hash in neigs.sorted_values().into_iter() {
        self.cone_coverage_approx_recur(depth_start, root_hash,
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                                &minmax_array, 0, &mut bmoc_builder);
      }
      return bmoc_builder; //.to_bmoc_packing();
    }
  }
  fn cone_coverage_approx_recur<F>(&self, depth: u8, hash: u64, shs_computer: &F, shs_minmax: &[MinMax],
                           recur_depth: u8, bmoc_builder: &mut BMOCBuilderUnsafe)
    where F: Fn((f64, f64)) -> f64  {
    let center = get_or_create(depth).center(hash);
    let shs = shs_computer(center);
    let MinMax{min, max} = shs_minmax[recur_depth as usize];
    if shs <= min {
      bmoc_builder.push(depth, hash, true);
    } else if shs <= max {
      if depth == self.depth {
        bmoc_builder.push(depth, hash, false);
      } else {
        let hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.cone_coverage_approx_recur(depth, hash, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_approx_recur(depth, hash | 1_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_approx_recur(depth, hash | 2_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_approx_recur(depth, hash | 3_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
      }
    }
  }
  
  /// cone_radius in radians
  /// TODO: find a better function!!
  #[inline]
  fn ncell_in_cone_upper_bound(&self, cone_radius: f64) -> usize {
    // cell_area = 4 pi / ncell
    // cone_area = pi r^2
    // => area_ratio = cone_area / cell_area = r^2 * ncell / 4 = r^2 * (ncell >> 2)
    let mut area_ratio = pow2(cone_radius) * ((self.n_hash >> 2) as f64);
    area_ratio += 1.0_f64; // to be sure the minimum = 1
    // We want:
    //  - ratio = 1    --> n = 9
    //  - ratio = +inf --> n = 1.25 (125% => + 25%) 
    //  - function shape = a/x + b
    // => a / 1 + b = 9  and a / +inf + b = 1.25 
    // => b = 1.25 and a = 7.75
    let correction_factor = 1.25_f64 + 7.75_f64 / area_ratio; 
    (correction_factor * area_ratio) as usize
  }

  #[inline]
  fn n_moc_cell_in_cone_upper_bound(&self, cone_radius: f64) -> usize {
    const TWICE_SQRT_3: f64 = 2.0_f64 * 1.73205080756887729352_f64; // sqrt(3)
    // cell_area = 4 * pi / ncell = 4 * pi / (3 * 4 * nside^2) = pi / (3 * nside^2) =  pi * r^2
    // cell_radius = r = 1 / (sqrt(3) * nside)
    // As a very simple and naive rule, we take 4x the number of cells needed to cover
    // the cone external annulus
    // Annulus area = 4 pi ((R + r)^2 - R^2) = 4 pi (r^2 + 2rR)
    // N cells = 4 pi (r^2 + 2rR) / 4 pi r^2 = 1 + 2 R/r = 1 + 2 * sqrt(3) * nside * R
    4_usize * (1_usize + (self.nside as f64 * TWICE_SQRT_3 * cone_radius + 0.99_f64) as usize)
  }



  /*pub fn cone_coverage(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      return self.allsky_bmoc();
    }
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius) {
      
    }
    // Normal case
    
  }*/
  
  /// Returns a hierarchical view of the list of cells overlapped by the given elliptical cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the elliptical cone.
  /// The algorithm is approximated: it may return false positive, 
  /// i.e. cells which are near from the elliptical cone but do not overlap it.
  /// To control the approximation, see the method 
  /// [cone_coverage_approx_custom](#method.elliptical_cone_coverage_custom)
  /// 
  /// # Input
  /// - `lon` the longitude of the center of the elliptical cone, in radians
  /// - `lat` the latitude of the center of the elliptical cone, in radians
  /// - `a` the semi-major axis of the elliptical cone, in radians
  /// - `b` the semi-minor axis of the elliptical cone, in radians
  /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
  /// 
  /// # Output
  /// - the list of cells overlapped by the given elliptical cone, in a BMOC 
  ///   (hierarchical view also telling if a cell is fully or partially covered).
  ///
  /// # Panics
  /// - if the semi-major axis is > PI/2
  /// 
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get_or_create(depth);
  ///
  /// let lon = 36.80105218_f64.to_radians();
  /// let lat = 56.78028536_f64.to_radians();
  /// let a = 14.93_f64.to_radians();
  /// let b = 4.93_f64.to_radians();
  /// let pa = 75.0_f64.to_radians();
  /// 
  /// let actual_res = nested3.elliptical_cone_coverage(lon, lat, a, b, pa);
  /// let expected_res: [u64; 16] = [27, 30, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57];
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn elliptical_cone_coverage(&self, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
    self.elliptical_cone_coverage_internal(lon, lat, a, b, pa).to_bmoc_packing()
  }
  
  /// Returns a hierarchical view of the list of cells overlapped by the given elliptical cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the elliptical cone.
  /// The algorithm is approximated: it may return false positive, 
  /// i.e. cells which are near from the cone but do not overlap it.
  /// To control the approximation, you can choose to perform the computations at a deeper depth
  /// using the `delta_depth` parameter.
  /// 
  /// # Input
  /// - `delta_depth` the difference between this Layer depth and the depth at which the computations
  ///   are made (should remain quite small).
  /// - `lon` the longitude of the center of the elliptical cone, in radians
  /// - `lat` the latitude of the center of the elliptical cone, in radians
  /// - `a` the semi-major axis of the elliptical cone, in radians
  /// - `b` the semi-minor axis of the elliptical cone, in radians
  /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
  /// 
  /// # Output
  /// - the list of cells overlapped by the given elliptical cone, in a BMOC 
  ///   (hierarchical view also telling if a cell is fully or partially covered).
  ///
  /// # Panics
  /// - if the semi-major axis is > PI/2
  /// - if this layer depth + `delta_depth` > the max depth (i.e. 29)
  ///
  pub fn elliptical_cone_coverage_custom(&self, delta_depth: u8, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
    if delta_depth == 0 {
      self.elliptical_cone_coverage(lon, lat, a, b, pa)
    } else {
      // TODO: change the algo not to put all cell in the MOC and pruning it
      get_or_create(self.depth + delta_depth)
        .elliptical_cone_coverage_internal(lon, lat, a, b, pa)
        .to_lower_depth_bmoc_packing(self.depth)
    }
  }

  pub fn elliptical_cone_coverage_internal(&self, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOCBuilderUnsafe {
    if a >= HALF_PI {
      panic!("Unable to handle ellipses with a semi-major axis > PI/2");
    }
    // Special case: the full sky is covered
    if b >= PI {
      return self.allsky_bmoc_builder();
    }
    // Common variable
    let sph_ellipse = EllipticalCone::new(lon, lat, a, b, pa);
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(a) {
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(a));
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, lon, lat, a);
      for h in 0..12 {
        self.elliptical_cone_coverage_recur(0, h ,&sph_ellipse, &distances, 0, &mut bmoc_builder);
      }
      return bmoc_builder;
    }
    // Normal case
    let depth_start = best_starting_depth(a);
    let root_layer = get_or_create(depth_start);
    let root_center_hash = root_layer.hash(lon, lat);
    // Small ellipse case
    if depth_start >= self.depth {
      let distance = largest_center_to_vertex_distance_with_radius(depth_start, lon, lat, a);
      let mut neigs: Vec<u64> = root_layer.neighbours(root_center_hash, true)
        .values_vec()
        .into_iter()
        .filter(|h| {
          let (l, b) = root_layer.center(*h);
          sph_ellipse.contains(l, b) || sph_ellipse.overlap_cone(l, b, distance)
        })
        .map(|h| h >> ((depth_start - &self.depth) << 1)) // h_to_lower_depth
        .collect();
      neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
      neigs.dedup();         // remove duplicates (vector must be sorted first)
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, false);
      }
      return bmoc_builder;
    } else {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, lon, lat, a);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(a));
      for &root_hash in neigs.sorted_values().into_iter() {
        self.elliptical_cone_coverage_recur(depth_start, root_hash, &sph_ellipse, &distances, 0, &mut bmoc_builder);
      }
      return bmoc_builder;
    }
  }
  
  fn elliptical_cone_coverage_recur(&self, depth: u8, hash: u64, 
                                       ellipse: &EllipticalCone, distances: &[f64],
                                       recur_depth: u8, bmoc_builder: &mut BMOCBuilderUnsafe) {
    let (lon, lat) = get_or_create(depth).center(hash);
    let distance = distances[recur_depth as usize];
    /*eprintln!("d: {}; h: {}; lon: {}, lat: {}; dist: {}; contains: {}; overlap: {}", 
             &depth, &hash, &lon.to_degrees(), &lat.to_degrees(), &distance.to_degrees(),
             &ellipse.contains_cone(lon, lat, distance),
             &ellipse.overlap_cone(lon, lat, distance));*/
    if ellipse.contains_cone(lon, lat, distance) {
      bmoc_builder.push(depth, hash, true);
    } else if ellipse.contains(lon, lat) || ellipse.overlap_cone(lon, lat, distance) {
      if depth == self.depth {
        let mut is_full = true;
        for (lon, lat) in self.vertices(hash).iter() {
          is_full &= ellipse.contains(*lon, *lat); // Not sure computation not done if is_full==false, to be verfied
        }
        bmoc_builder.push(depth, hash, is_full);
      } else {
        let hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.elliptical_cone_coverage_recur(depth,      hash        , ellipse, distances, recur_depth, bmoc_builder);
        self.elliptical_cone_coverage_recur(depth, hash | 1_u64, ellipse, distances, recur_depth, bmoc_builder);
        self.elliptical_cone_coverage_recur(depth, hash | 2_u64, ellipse, distances, recur_depth, bmoc_builder);
        self.elliptical_cone_coverage_recur(depth, hash | 3_u64, ellipse, distances, recur_depth, bmoc_builder);
      }
    }
  }
  
  /// Returns a hierarchical view of the list of cells overlapped by the given polygon.
  /// Self intersecting polygons are supported.
  /// The BMOC also tells if the cell if fully or partially overlapped by the polygon.
  /// 
  /// If you want the complementary solution, apply the NOT operator on the BMOC.
  /// 
  /// This method supports both an *exact* (need more tests) and an *approximated* solution.
  /// The second one being faster (TODO: measure and provided perf differences as a function of the
  /// number of vertices in the polygon).
  /// 
  /// The approximation is the following one: 
  /// > when testing the intersection between a polygon segment and an HEALPix cell edge
  /// > we consider that each edge of the HEALPix cell is on a great-circle arc (which is not
  /// > true, especially a low resolutions).
  /// 
  /// For the exact solution:
  /// > for each polygon segment, we first test if the segment contains a 'special point', 
  /// > if it is the case, we add it to the list of cell number computed from each polygin vertex
  /// > A 'special point' is a point such that in the HEALPix projection Euclidean plane
  /// ```math
  /// \mathrm{d}\DeltaX(z) / \mathrm{d}Y = \pm 1
  /// ``` 
  ///
  /// # Input
  /// - `vertices` the list of vertices (in a slice) coordinates, in radians
  ///              `[(lon, lat), (lon, lat), ..., (lon, lat)]`
  /// - `exact_solution` if set 
  /// 
  /// # Output
  /// - the list of cells overlapped by the given polygon, in a BMOC (hierarchical view also telling
  ///   if a cell is fully or partially covered).
  /// 
  /// # Example
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// use cdshealpix::nested::{get_or_create, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get_or_create(depth);
  /// 
  /// let actual_res =nested3.polygon_coverage(&[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)], false);
  /// let expected_res: [u64; 8] = [304, 305, 306, 307, 308, 310, 313, 316];
  /// 
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn polygon_coverage(&self, vertices: &[(f64, f64)], exact_solution: bool) -> BMOC {
    let poly = Polygon::new(
      vertices.iter().map(|(lon, lat)| LonLat { lon: *lon, lat: *lat} )
        .collect::<Vec<LonLat>>().into_boxed_slice()
    );
    let bounding_cone: Cone = Cone::bounding_cone(poly.vertices());
    let mut depth_start = 0;
    let neigs: Vec<u64> = if !has_best_starting_depth(bounding_cone.radius()) {
      (0..12).collect()
    } else {
      depth_start = best_starting_depth(bounding_cone.radius()).min(self.depth);
      let root_layer = get_or_create(depth_start);
      let LonLat{lon, lat} = bounding_cone.center().lonlat();
      let center_hash = root_layer.hash(lon, lat);
      let mut neigs: Vec<u64> = root_layer.neighbours(center_hash, true).values_vec();
      neigs.sort_unstable();
      neigs
    };
    // Compute and sort the list of cells containing at least one polygon vertex
    let mut sorted_poly_vertices_hash = self.hashs_vec(poly.vertices());
    // Special treatment for the exact solution
    if exact_solution {
      let vertices= poly.vertices();
      let mut left= &vertices[vertices.len() - 1];
      for right in vertices {
        let special_lonlats = arc_special_points(left, right, 1.0e-14, 20);
        // println!("special_lonlats: {:?}", &special_lonlats);
        sorted_poly_vertices_hash.append(&mut self.hashs_vec(&special_lonlats));
        left = right;
      }
    }
    // Back to general case
    sorted_poly_vertices_hash.sort_unstable();
    sorted_poly_vertices_hash.dedup();
    let sorted_poly_vertices_hash = sorted_poly_vertices_hash.into_boxed_slice();
    // Build the list (removing duplicated) for all deltaDepth?
    let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 10_000_usize);
    for root_hash in neigs { 
      self.polygon_coverage_recur(&mut bmoc_builder, depth_start, root_hash, &poly, &sorted_poly_vertices_hash);
    }
    bmoc_builder.to_bmoc()
  }

  fn polygon_coverage_recur(&self, moc_builder: &mut BMOCBuilderUnsafe, depth: u8, hash: u64,
  poly: &Polygon, sorted_poly_vertices_hash: &[u64]) {
    if is_in_list(depth, hash, self.depth, sorted_poly_vertices_hash) {
      if depth == self.depth {
        moc_builder.push(depth, hash, false);
      } else {
        let hash = hash << 2;
        let depth = depth + 1;
        self.polygon_coverage_recur(moc_builder, depth, hash    , poly, sorted_poly_vertices_hash);
        self.polygon_coverage_recur(moc_builder, depth, hash | 1, poly, sorted_poly_vertices_hash);
        self.polygon_coverage_recur(moc_builder, depth, hash | 2, poly, sorted_poly_vertices_hash);
        self.polygon_coverage_recur(moc_builder, depth, hash | 3, poly, sorted_poly_vertices_hash);
      }
    } else {
      let (n_vertices_in_poly, poly_vertices) = n_vertices_in_poly(depth, hash, poly);
      if n_vertices_in_poly == 4 {
        moc_builder.push(depth, hash, true);
      } else if n_vertices_in_poly > 0 || has_intersection(poly, poly_vertices) {
        if depth == self.depth {
          moc_builder.push(depth, hash, false);
        } else {
          // I known, I don't like this code repetition, TODO: see how to remove it
          let hash = hash << 2;
          let depth = depth + 1;
          self.polygon_coverage_recur(moc_builder, depth, hash    , poly, sorted_poly_vertices_hash);
          self.polygon_coverage_recur(moc_builder, depth, hash | 1, poly, sorted_poly_vertices_hash);
          self.polygon_coverage_recur(moc_builder, depth, hash | 2, poly, sorted_poly_vertices_hash);
          self.polygon_coverage_recur(moc_builder, depth, hash | 3, poly, sorted_poly_vertices_hash);
        }
      }
    }
  }

  fn hashs<T: LonLatT>(&self, poly_vertices: &[T]) -> Box<[u64]> {
    poly_vertices.iter().map(|coo| self.hash(coo.lon(), coo.lat()))
      .collect::<Vec<u64>>().into_boxed_slice()
  }

  fn hashs_vec<T: LonLatT>(&self, poly_vertices: &[T]) -> Vec<u64> {
    poly_vertices.iter().map(|coo| self.hash(coo.lon(), coo.lat()))
      .collect::<Vec<u64>>()
  }
  
  fn allsky_bmoc(&self) -> BMOC {
    self.allsky_bmoc_builder().to_bmoc()
  }

  fn allsky_bmoc_builder(&self) -> BMOCBuilderUnsafe {
    let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 12);
    bmoc_builder.push_all(0_u8, 0_u64, 12_u64, true);
    bmoc_builder
  }
}

/// Returns the hash value of the base cell (the depth 0 cell hash value) given its '(i, j)'
/// coordinates in the projected, rotated, scaled plane.
/// Here we suppose that we get `(i, j)` from the center of a cell, so that we do not have to take
/// care of border effects (see `depth0_bits` accounting for border effects).
#[inline]
fn depth0_hash_unsafe(i: u8, j: u8) -> u8 {
  let k = 5_i8 - (i + j) as i8;
  (((k << 2) + ( ((i as i8) + ((k - 1) >> 7)) & 3_i8)) as u8)
}

/// Returns the hash value of the cell of depth this layer depth + the given `delta_depth`
/// located in the corner of given direction in the given cell.
/// ```rust
/// use cdshealpix::compass_point::{Cardinal};
/// use cdshealpix::nested::{internal_corner};
///
/// let delta_depth = 1;
/// assert_eq!(0 , internal_corner(0, delta_depth, &Cardinal::S));
/// assert_eq!(1 , internal_corner(0, delta_depth, &Cardinal::E));
/// assert_eq!(2 , internal_corner(0, delta_depth, &Cardinal::W));
/// assert_eq!(3 , internal_corner(0, delta_depth, &Cardinal::N));
/// ```
pub fn internal_corner(hash: u64, delta_depth: u8, direction: &Cardinal) -> u64 {
  match *direction {
    Cardinal::S => internal_corner_south(hash, delta_depth),
    Cardinal::E => internal_corner_east(hash, delta_depth),
    Cardinal::N => internal_corner_north(hash, delta_depth),
    Cardinal::W => internal_corner_west(hash, delta_depth),
  }
}

/// Returns the hash value of the cell of depth the given hash depth + the given `delta_depth`
/// located in the south corner of the given cell.
pub fn internal_corner_south(hash: u64, delta_depth: u8) -> u64 {
  hash << (delta_depth << 1)
}

/// Returns the hash value of the cell of depth the given hash depth + the given `delta_depth`
/// located in the east corner of the given cell.
pub fn internal_corner_east(hash: u64, delta_depth: u8) -> u64 {
  (hash << (delta_depth << 1)) | x_mask(delta_depth)
}

/// Returns the hash value of the cell of depth the given hash depth + the given `delta_depth`
/// located in the west corner of the given cell.
pub fn internal_corner_west(hash: u64, delta_depth: u8) -> u64 {
  (hash << (delta_depth << 1)) | y_mask(delta_depth)
}

/// Returns the hash value of the cell of depth the given hash depth + the given `delta_depth`
/// located in the north corner of the given cell.
pub fn internal_corner_north(hash: u64, delta_depth: u8) -> u64 {
  (hash << (delta_depth << 1)) | xy_mask(delta_depth)
}


/// Returns the hash values of the cells of depth this layer depth + the given `delta_depth`
/// located in the internal edge of given direction in the given cell.
/// 
/// # Info
/// The returned vector is sorted.
/// 
/// ```rust
/// use cdshealpix::compass_point::{Ordinal};
/// use cdshealpix::nested::{internal_edge_part};
///
/// 
/// let delta_depth = 1;
/// assert_eq!(vec![0, 1].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::SE));
/// assert_eq!(vec![0, 2].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::SW));
/// assert_eq!(vec![1, 3].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::NE));
/// assert_eq!(vec![2, 3].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::NW));
/// 
/// let delta_depth = 2;
/// assert_eq!(vec![ 0,  1,  4,  5].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::SE));
/// assert_eq!(vec![ 0,  2,  8, 10].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::SW));
/// assert_eq!(vec![ 5,  7, 13, 15].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::NE));
/// assert_eq!(vec![10, 11, 14, 15].into_boxed_slice() , internal_edge_part(0, delta_depth, &Ordinal::NW));
/// ```
pub fn internal_edge_part(hash: u64, delta_depth: u8, direction: &Ordinal) -> Box<[u64]> {
  match *direction {
    Ordinal::SE => internal_edge_southeast(hash, delta_depth),
    Ordinal::SW => internal_edge_southwest(hash, delta_depth),
    Ordinal::NE => internal_edge_northeast(hash, delta_depth),
    Ordinal::NW => internal_edge_northwest(hash, delta_depth),
  }
}

/// Returns the hash values of the cells of depth the given hash depth + the given `delta_depth`
/// located in the southeast internal edge of the given cell.
/// # Info
/// The returned vector is sorted.
pub fn internal_edge_southeast(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let mut v: Vec<u64> = Vec::with_capacity(nside as usize);
  hash <<= delta_depth << 1;
  for x in 0..nside {
    v.push(hash | get_zoc(delta_depth).i02h(x));
  }
  v.into_boxed_slice()
}

/// Returns the hash values of the cells of depth the given hash depth + the given `delta_depth`
/// located in the southwest internal edge of the given cell.
/// # Info
/// The returned vector is sorted.
pub fn internal_edge_southwest(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let mut v: Vec<u64> = Vec::with_capacity(nside as usize);
  hash <<= delta_depth << 1;
  for y in 0..nside {
    v.push(hash | get_zoc(delta_depth).oj2h(y));
  }
  v.into_boxed_slice()
}

/// Returns the hash values of the cells of depth the given hash depth + the given `delta_depth`
/// located in the northeast internal edge of the given cell.
/// # Info
/// The returned vector is sorted.
pub fn internal_edge_northeast(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let mut v: Vec<u64> = Vec::with_capacity(nside as usize);
  hash <<= delta_depth << 1;
  let x_bits = get_zoc(delta_depth).i02h(nside - 1);
  for y in 0..nside {
    v.push(hash | get_zoc(delta_depth).oj2h(y) | x_bits);
  }
  v.into_boxed_slice()
}

/// Returns the hash values of the cells of depth the given hash depth + the given `delta_depth`
/// located in the northwest internal edge of the given cell.
/// # Info
/// The returned vector is sorted.
pub fn internal_edge_northwest(mut hash: u64, delta_depth: u8) -> Box<[u64]> {
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let mut v: Vec<u64> = Vec::with_capacity(nside as usize);
  hash <<= delta_depth << 1;
  let y_bits = get_zoc(delta_depth).oj2h(nside - 1);
  for x in 0..nside {
    v.push(hash | get_zoc(delta_depth).i02h(x) | y_bits);
  }
  v.into_boxed_slice()
}

/// Same as [internal_edge_part](fn.internal_edge_part.html) except that the result is appended
/// to the given vec.
pub fn append_internal_edge_part(hash: u64, delta_depth: u8, direction: &Ordinal, result: &mut Vec<u64>) {
  match *direction {
    Ordinal::SE => append_internal_edge_southeast(hash, delta_depth, result),
    Ordinal::SW => append_internal_edge_southwest(hash, delta_depth, result),
    Ordinal::NE => append_internal_edge_northeast(hash, delta_depth, result),
    Ordinal::NW => append_internal_edge_northwest(hash, delta_depth, result),
  }
}

/// Same as [internal_edge_southeast](fn.internal_edge_southeast.html) except that
/// the result is appended to the given vec.
pub fn append_internal_edge_southeast(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
  hash <<= delta_depth << 1;
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  for x in 0..nside {
    result.push(hash | get_zoc(delta_depth).i02h(x));
  }
}

/// Same as [internal_edge_southwest](fn.internal_edge_southwest.html) except that
/// the result is appended to the given vec.
pub fn append_internal_edge_southwest(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
  hash <<= delta_depth << 1;
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  for y in 0..nside {
    result.push(hash | get_zoc(delta_depth).oj2h(y));
  }
}

/// Same as [internal_edge_northeast](fn.internal_edge_northeast.html) except that
/// the result is appended to the given vec.
pub fn append_internal_edge_northeast(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
  hash <<= delta_depth << 1;
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let x_bits = get_zoc(delta_depth).i02h(nside - 1);
  for y in 0..nside {
    result.push(hash | get_zoc(delta_depth).oj2h(y) | x_bits);
  }
}

/// Same as [internal_edge_northwest](fn.internal_edge_northwest.html) except that
/// the result is appended to the given vec.
pub fn append_internal_edge_northwest(mut hash: u64, delta_depth: u8, result: &mut Vec<u64>) {
  hash <<= delta_depth << 1;
  let nside = 1_u32 << delta_depth; // 2^deltaDepth
  let y_bits = get_zoc(delta_depth).oj2h(nside - 1);
  for x in 0..nside {
    result.push(hash | get_zoc(delta_depth).i02h(x) | y_bits);
  }
}


/// # Panics
/// If the given Main Wind is the Center (i.e. is neither Ordinal nor Cardinal).
fn add_sorted_internal_edge_element(hash: u64, delta_depth: u8, direction: MainWind, ext_direction: &MainWind, result: &mut ExternalEdge) {
  if direction.is_cardinal() {
    let cardinal = direction.to_cardinal();
    result.set_corner(&ext_direction.to_cardinal(), 
                      internal_corner(hash, delta_depth, &cardinal));
  } else if direction.is_ordinal() {
    let ordinal = direction.to_ordinal();
    result.set_edge(&ext_direction.to_ordinal(), 
                    internal_edge_part(hash, delta_depth, &ordinal));
  } else {
    panic!("Main wind {:?} is neither ordinal not cardinal", &direction);
  }
}

/// # Panics
/// If the given Main Wind is the Center (i.e. is neither Ordinal nor Cardinal).
fn append_sorted_internal_edge_element(hash: u64, delta_depth: u8, direction: MainWind, result: &mut Vec<u64>) {
  if direction.is_cardinal() {
    result.push(internal_corner(hash, delta_depth, &direction.to_cardinal()));
  } else if direction.is_ordinal() {
    append_internal_edge_part(hash, delta_depth, &direction.to_ordinal(), result);
  } else {
    panic!("Main wind {:?} is neither ordinal not cardinal", &direction);
  }
}





fn is_in_list(depth: u8, hash: u64, depth_hashs: u8, sorted_hashs: &[u64]) -> bool {
  let twice_delta_depth = (depth_hashs - depth) << 1;
  let hash_at_depth_max = hash << twice_delta_depth;
  match sorted_hashs.binary_search(&hash_at_depth_max) {
    Ok(_) => true,
    Err(i) => {
           (i < sorted_hashs.len() && (sorted_hashs[i] >> twice_delta_depth) == hash)
        || (i > 0_usize && (sorted_hashs[i - 1] >> twice_delta_depth) == hash)
    },
  }
}

fn n_vertices_in_poly(depth: u8, hash: u64, poly: &Polygon) -> (u8, [Coo3D; 4]) {
  let [(l_south, b_south), (l_east, b_east), (l_north, b_north), (l_west, b_west)] = vertices(depth, hash);
  let vertices = [
    Coo3D::from_sph_coo(l_south, b_south),
    Coo3D::from_sph_coo(l_east, b_east),
    Coo3D::from_sph_coo(l_north, b_north),
    Coo3D::from_sph_coo(l_west, b_west)
  ];
  let n_vertices_in_poly = 
      (poly.contains(&vertices[0]) as u8)
    + (poly.contains(&vertices[1]) as u8)
    + (poly.contains(&vertices[2]) as u8)
    + (poly.contains(&vertices[3]) as u8);
  (n_vertices_in_poly, vertices)
}

fn has_intersection(poly: &Polygon, vertices: [Coo3D; 4]) -> bool {
       poly.intersect_great_circle_arc(&vertices[2], &vertices[1]) // N vs E
    || poly.intersect_great_circle_arc(&vertices[0], &vertices[1]) // S vs E
    || poly.intersect_great_circle_arc(&vertices[3], &vertices[2]) // W vs N
    || poly.intersect_great_circle_arc(&vertices[3], &vertices[0]) // W vs S
}


struct MinMax {
  min: f64,
  max: f64,
}

struct HashParts {
  d0h: u8, // base cell number (depth 0 hash value)
  i: u32, // in the base cell, z-order curve coordinate along the x-axis
  j: u32, // in the base cell, z-order curve coordinate along the x-axis
}

struct HashBits {
  d0h: u64, // base cell number (depth 0 hash value) bits
  i: u64,   // in the base cell, z-order curve coordinate along the x-axis bits
  j: u64,   // in the base cell, z-order curve coordinate along the y-axis bits
}

#[inline]
const fn discretize(xy: (f64, f64)) -> (u64, u64) {
  (xy.0 as u64, xy.1 as u64)
}

#[inline]
const fn rotate45_scale2(i_in_d0h: u32, j_in_d0h: u32) -> (i32, i32) {
  (i_in_d0h as i32 - j_in_d0h as i32, (i_in_d0h + j_in_d0h) as i32)
}

/// offset_x in [0, 7], odd for polar caps, even for equatorial region
/// offset_y in [-1, 1], -1 or 1 for polar caps, 0 for equatorial region
#[inline]
const fn compute_base_cell_center_offsets_in_8x3_grid(d0h: u8) -> (u8, i8) { 
  let offset_y = 1 - div4_quotient(d0h) as i8;
  let mut offset_x = (div4_remainder(d0h)) << 1u8;
  // +1 if the base cell is not equatorial
  offset_x |= (offset_y & 1_i8) as  u8;
  (offset_x, offset_y)
}

#[inline]
fn apply_base_cell_center_offsets(xy: &mut (f64, f64), offset_x: u8, offset_y: i8) {
  let (ref mut x, ref mut y) = *xy;
  *x += offset_x as f64;
  *y += offset_y as f64;
  // If x < 0, then x += 8; (happens only in case of base cell 4)
  *x += ((f64::to_bits(*x) & super::F64_SIGN_BIT_MASK) >> 60) as f64;
}

#[inline]
const fn bits_2_hash(d0h_bits: u64, i_in_d0h_bits: u64, j_in_d0h_bits: u64) -> u64 {
  d0h_bits | i_in_d0h_bits | j_in_d0h_bits
}

/// x / 2
#[inline]
fn div2_quotient<T: Shr<u8, Output=T>>(x: T) -> T {
  // x >> 1
  x.shr(1)
}

/// x modulo 2
#[inline]
const fn div2_remainder(x: u64) -> u64 {
  x & 1
}

/// x / 4
#[inline]
const fn div4_quotient(x: u8) -> u8 {
  x >> 2
}

/// x modulo 4
#[inline]
const fn div4_remainder(x: u8) -> u8 {
  x & 3
}

/// mask ...010101
/// ```rust
/// use cdshealpix::nested::{x_mask};
/// assert_eq!(x_mask(3), 0b00010101);
/// ```
#[inline]
pub const fn x_mask(depth: u8) -> u64 {
  0x5555555555555555_u64 >> (64 - (depth << 1))
}

/// mask ...101010
/// ```rust
/// use cdshealpix::nested::{y_mask, x_mask};
/// assert_eq!(y_mask(3), 0b00101010);
/// assert_eq!(y_mask(3), x_mask(3) << 1);
/// ```
#[inline]
pub const fn y_mask(depth: u8) -> u64{
  0xAAAAAAAAAAAAAAAA_u64 >> (64 - (depth << 1))
}

/// mask ...111111
/// ```rust
/// use cdshealpix::nested::{xy_mask};
/// assert_eq!(xy_mask(3), 0b00111111);
/// let depth = 3_u8;
/// assert_eq!(0xFFFFFFFFFFFFFFFF_u64 >> (64 - (depth << 1)), (1_u64 << (depth << 1)) - 1_u64);
/// ```
#[inline]
pub const fn xy_mask(depth: u8) -> u64 {
  0xFFFFFFFFFFFFFFFF_u64 >> (64 - (depth << 1))
}

#[inline]
fn to_shs_min_max_array(cone_radius: f64, distances: Box<[f64]>) -> Box<[MinMax]> {
  let v: Vec<MinMax> = distances.iter()
    .map(|d| to_shs_min_max(cone_radius, *d))
    .collect();
  v.into_boxed_slice()
}
#[inline]
fn to_shs_min_max(cone_radius: f64, distance: f64) -> MinMax {
  MinMax {
    min: if cone_radius < distance { 
      0.0
    } else {
      to_squared_half_segment(cone_radius - distance)
    },
    max: to_squared_half_segment(cone_radius + distance),
  }
}

/// - `h` stands for `hash`
/// - `shs` stands for `squared half segment`
#[inline]
fn h_to_h_and_shs(cone_lon: f64, cone_lat: f64, cos_cone_lat:f64, layer: &'static Layer)
                  -> impl FnMut(&u64) -> (u64, f64) {
  move |&hash| {
    let (lon, lat) = layer.center(hash);
    (hash, squared_half_segment(lon - cone_lon, lat - cone_lat, lon.cos(), cos_cone_lat))
  }
}

/// - `shs` stands for `squared half segment`
#[inline]
fn shs_lower_than(shs_max: f64) -> impl FnMut(&(u64, f64)) -> bool {
  move |&(_, shs)| shs <= shs_max
}

/// The returned closure computes the 
#[inline]
fn shs_computer(cone_lon: f64, cone_lat: f64, cos_cone_lat: f64) -> impl Fn((f64, f64)) -> f64 {
  move |(lon, lat)| {
    squared_half_segment(lon - cone_lon, lat - cone_lat, lat.cos(), cos_cone_lat)
  }
}

/* Commented because not used so far
/// The returned closure computes the squared half segment between 
/// the center of a given hash (the closure parameter) ,at the depth of the given layer,
/// and the given point on the unit sphere.
/// We ask in argument the cosine of the latitude since it is a costly operation that may have already
/// been performed.
#[inline]
fn shs_computer_4layer(lon: f64, lat: f64, cos_lat: f64, layer: &'static Layer) -> impl Fn(u64) -> (f64) {
  debug_assert!(cos_lat == lat.cos());
  let func = shs_computer(lon, lat, cos_lat);
  move |hash| {
    func(layer.center(hash))
  }
}*/

#[cfg(test)]
mod tests {
  use super::*;

  fn build_lupt(depth: u8) -> [[u64; 6]; 6] {
    let depth_x2: u8= depth << 1u8;
    let mut y_mask = 0u64;
    let mut xy_mask = 0u64;
    if depth > 0 {
      xy_mask = (1u64 << depth_x2) - 1u64;
      y_mask = 0x5555555555555555u64 >> (64 - depth_x2); // ...0101
      y_mask <<= 1;                                      // ...1010
    }
    let null = 0xFFFFFFFFFFFFFFFFu64;
    let bc00 =  0u64 << depth_x2;
    let bc01 =  1u64 << depth_x2;
    let bc02 =  2u64 << depth_x2;
    let bc03 =  3u64 << depth_x2;
    let bc04 =  4u64 << depth_x2;
    let bc05 =  5u64 << depth_x2;
    let bc06 =  6u64 << depth_x2;
    let bc07 =  7u64 << depth_x2;
    let bc08 =  8u64 << depth_x2;
    let bc09 =  9u64 << depth_x2;
    let bc10 = 10u64 << depth_x2;
    let bc11 = 11u64 << depth_x2;
    let mg0y = (0u64 << depth_x2) | y_mask;
    let mg1y = (1u64 << depth_x2) | y_mask;
    let mg2y = (2u64 << depth_x2) | y_mask;
    let mg3y = (3u64 << depth_x2) | y_mask;
    let m0xy = (0u64 << depth_x2) | xy_mask;
    let m1xy = (1u64 << depth_x2) | xy_mask;
    let m2xy = (2u64 << depth_x2) | xy_mask;
    let m3xy = (3u64 << depth_x2) | xy_mask;
    [
      [null, null, null, bc08, bc04, null], //   ----> y-axis
      [null, null, bc09, bc05, bc00, mg0y], //  |
      [null, bc10, bc06, bc01, mg1y, m0xy], //  |
      [bc11, bc07, bc02, mg2y, m1xy, null], //  v
      [bc04, bc03, mg3y, m2xy, null, null], // x-axis
      [null, mg0y, m3xy, null, null, null]
    ]
  }
  
  /*#[test]
  fn testok_d0h_bits() {
    let depth = 0_u8;
    let layer = get_or_create(depth);
    let lupm = build_lupt(depth);
    let ijs = [
      (0_u8, 3_u8), (0_u8, 4_u8),
      (1_u8, 2_u8), (1_u8, 3_u8), (1_u8, 4_u8), (1_u8, 5_u8),
      (2_u8, 1_u8), (2_u8, 2_u8), (2_u8, 3_u8), (2_u8, 4_u8), (2_u8, 5_u8),
      (3_u8, 0_u8), (3_u8, 1_u8), (3_u8, 2_u8), (3_u8, 3_u8), (3_u8, 4_u8),
      (4_u8, 0_u8), (4_u8, 1_u8), (4_u8, 2_u8), (4_u8, 3_u8),
      (5_u8, 1_u8), (5_u8, 2_u8),
    ];
    for (i, j) in ijs.iter() {
      assert_eq!(lupm[*i as usize][*j as usize], layer.depth0_bits(*i, *j/*, &mut (0_u64, 0_u64)*/, (*i as f64, *j as f64 + 1e-14)/*, 0.0, 0.0*/));
    }
  }*/
  
  #[test]
  fn testok_hash_d0() {
    let layer = get_or_create(0);
    assert_eq!(4_u64, layer.hash(0.0_f64.to_radians(), 0.0_f64.to_radians()));
    assert_eq!(5_u64, layer.hash(90.0_f64.to_radians(), 0.0_f64.to_radians()));
    assert_eq!(6_u64, layer.hash(180.0_f64.to_radians(), 0.0_f64.to_radians()));
    assert_eq!(7_u64, layer.hash(270.0_f64.to_radians(), 0.0_f64.to_radians()));
    
    assert_eq!(0_u64, layer.hash(45.0_f64.to_radians(), 41.0_f64.to_radians()));
    assert_eq!(1_u64, layer.hash(135.0_f64.to_radians(), 41.0_f64.to_radians()));
    assert_eq!(2_u64, layer.hash(225.0_f64.to_radians(), 41.0_f64.to_radians()));
    assert_eq!(3_u64, layer.hash(315.0_f64.to_radians(), 41.0_f64.to_radians()));

    assert_eq!(8_u64, layer.hash(45.0_f64.to_radians(), -41.0_f64.to_radians()));
    assert_eq!(9_u64, layer.hash(135.0_f64.to_radians(), -41.0_f64.to_radians()));
    assert_eq!(10_u64, layer.hash(225.0_f64.to_radians(), -41.0_f64.to_radians()));
    assert_eq!(11_u64, layer.hash(315.0_f64.to_radians(), -41.0_f64.to_radians()));
  }
  
  #[test]
  fn testok_hash() {
    let layer = get_or_create(3);
    let hash = layer.hash(333.5982493968911_f64.to_radians(), -25.919634217871433_f64.to_radians());
    assert_eq!(735_u64, hash);
  }

  #[test]
  fn testok_hash_2() {
    let layer = get_or_create(0);
    // ra = 179.99999999999998633839 deg
    // de = -48.13786699999999889561 deg
    let hash = layer.hash(3.141592653589793, -0.8401642740371252);
    // The difference comes from the fact that:
    //   ra * 4/pi = 4.0 instead of 3.99999999999999969640 due to numerical approximations.
    // In v1, it is a particular case (k=3) in depth0_bits, but we do not handle all of them
    if hash == 9_u64 {
      assert_eq!(9_u64, hash);   // with hash_v1
    } else {
      assert_eq!(10_u64, hash);  // with hash_v2
    }
  }

  #[test]
  fn testok_hash_3() {
    let layer = get_or_create(0);
    // ra = 89.99999999999999889877 deg
    // de = -42.68491599999999973256 deg
    let hash = layer.hash(1.5707963267948966, -0.7449923251372079);
    // The difference comes from the fact that:
    //   ra * 4/pi = 2.0 instead of 1.99999999999999997552 due to numerical approximations.
    // In v1, it is a particular case (k=3) in depth0_bits, but we do not handle all of them
    if hash == 8_u64 {
      assert_eq!(8_u64, hash);
    } else {
      assert_eq!(9_u64, hash);
    }
  }

  #[test]
  fn testok_hash_4() {
    let layer = get_or_create(3);
    let ra = 180.0_f64;
    let dec = -45.85_f64;
    let hash = layer.hash(ra.to_radians(), dec.to_radians());
    assert_eq!(682_u64, hash);
  }
  
  #[test]
  fn testok_neighbour_d0() {
    let layer = get_or_create(0);
    // North polar cap
    // - 0
    assert_eq!(0_u64, layer.neighbour(0,  C).unwrap());
    assert_eq!(2_u64, layer.neighbour(0,  N).unwrap());
    assert_eq!(1_u64, layer.neighbour(0, NE).unwrap());
    assert_eq!(3_u64, layer.neighbour(0, NW).unwrap());
    assert_eq!(8_u64, layer.neighbour(0,  S).unwrap());
    assert_eq!(5_u64, layer.neighbour(0, SE).unwrap());
    assert_eq!(4_u64, layer.neighbour(0, SW).unwrap());
    assert_eq!( None, layer.neighbour(0,  E));
    assert_eq!( None, layer.neighbour(0,  W));
    // - 1
    assert_eq!(1_u64, layer.neighbour(1,  C).unwrap());
    assert_eq!(3_u64, layer.neighbour(1,  N).unwrap());
    assert_eq!(2_u64, layer.neighbour(1, NE).unwrap());
    assert_eq!(0_u64, layer.neighbour(1, NW).unwrap());
    assert_eq!(9_u64, layer.neighbour(1,  S).unwrap());
    assert_eq!(6_u64, layer.neighbour(1, SE).unwrap());
    assert_eq!(5_u64, layer.neighbour(1, SW).unwrap());
    assert_eq!( None, layer.neighbour(1,  E));
    assert_eq!( None, layer.neighbour(1,  W));
    // - 2
    assert_eq!(2_u64, layer.neighbour(2,  C).unwrap());
    assert_eq!(0_u64, layer.neighbour(2,  N).unwrap());
    assert_eq!(3_u64, layer.neighbour(2, NE).unwrap());
    assert_eq!(1_u64, layer.neighbour(2, NW).unwrap());
    assert_eq!(10_u64, layer.neighbour(2,  S).unwrap());
    assert_eq!(7_u64, layer.neighbour(2, SE).unwrap());
    assert_eq!(6_u64, layer.neighbour(2, SW).unwrap());
    assert_eq!( None, layer.neighbour(2,  E));
    assert_eq!( None, layer.neighbour(2,  W));
    // - 3
    assert_eq!( 3_u64, layer.neighbour(3,  C).unwrap());
    assert_eq!( 1_u64, layer.neighbour(3,  N).unwrap());
    assert_eq!( 0_u64, layer.neighbour(3, NE).unwrap());
    assert_eq!( 2_u64, layer.neighbour(3, NW).unwrap());
    assert_eq!(11_u64, layer.neighbour(3,  S).unwrap());
    assert_eq!( 4_u64, layer.neighbour(3, SE).unwrap());
    assert_eq!( 7_u64, layer.neighbour(3, SW).unwrap());
    assert_eq!(  None, layer.neighbour(3,  E));
    assert_eq!(  None, layer.neighbour(3,  W));
    // Equatorial region
    // - 4
    assert_eq!( 4_u64, layer.neighbour(4,  C).unwrap());
    assert_eq!(  None, layer.neighbour(4,  N));
    assert_eq!( 0_u64, layer.neighbour(4, NE).unwrap());
    assert_eq!( 3_u64, layer.neighbour(4, NW).unwrap());
    assert_eq!(  None, layer.neighbour(4,  S));
    assert_eq!( 8_u64, layer.neighbour(4, SE).unwrap());
    assert_eq!(11_u64, layer.neighbour(4, SW).unwrap());
    assert_eq!( 5_u64, layer.neighbour(4,  E).unwrap());
    assert_eq!( 7_u64, layer.neighbour(4,  W).unwrap());
    // - 5
    assert_eq!(5_u64, layer.neighbour(5,  C).unwrap());
    assert_eq!(  None, layer.neighbour(5,  N));
    assert_eq!(1_u64, layer.neighbour(5, NE).unwrap());
    assert_eq!(0_u64, layer.neighbour(5, NW).unwrap());
    assert_eq!(  None, layer.neighbour(5,  S));
    assert_eq!(9_u64, layer.neighbour(5, SE).unwrap());
    assert_eq!(8_u64, layer.neighbour(5, SW).unwrap());
    assert_eq!(6_u64, layer.neighbour(5,  E).unwrap());
    assert_eq!(4_u64, layer.neighbour(5,  W).unwrap());
    // - 6
    assert_eq!( 6_u64, layer.neighbour(6,  C).unwrap());
    assert_eq!(  None, layer.neighbour(6,  N));
    assert_eq!( 2_u64, layer.neighbour(6, NE).unwrap());
    assert_eq!( 1_u64, layer.neighbour(6, NW).unwrap());
    assert_eq!(  None, layer.neighbour(6,  S));
    assert_eq!(10_u64, layer.neighbour(6, SE).unwrap());
    assert_eq!( 9_u64, layer.neighbour(6, SW).unwrap());
    assert_eq!( 7_u64, layer.neighbour(6,  E).unwrap());
    assert_eq!( 5_u64, layer.neighbour(6,  W).unwrap());
    // - 7
    assert_eq!( 7_u64, layer.neighbour(7,  C).unwrap());
    assert_eq!(  None, layer.neighbour(7,  N));
    assert_eq!( 3_u64, layer.neighbour(7, NE).unwrap());
    assert_eq!( 2_u64, layer.neighbour(7, NW).unwrap());
    assert_eq!(  None, layer.neighbour(7,  S));
    assert_eq!(11_u64, layer.neighbour(7, SE).unwrap());
    assert_eq!(10_u64, layer.neighbour(7, SW).unwrap());
    assert_eq!( 4_u64, layer.neighbour(7,  E).unwrap());
    assert_eq!( 6_u64, layer.neighbour(7,  W).unwrap());
    //  South polar cap
    // - 8
    assert_eq!( 8_u64, layer.neighbour(8,  C).unwrap());
    assert_eq!( 0_u64, layer.neighbour(8,  N).unwrap());
    assert_eq!( 5_u64, layer.neighbour(8, NE).unwrap());
    assert_eq!( 4_u64, layer.neighbour(8, NW).unwrap());
    assert_eq!(10_u64, layer.neighbour(8,  S).unwrap());
    assert_eq!( 9_u64, layer.neighbour(8, SE).unwrap());
    assert_eq!(11_u64, layer.neighbour(8, SW).unwrap());
    assert_eq!(  None, layer.neighbour(8,  E));
    assert_eq!(  None, layer.neighbour(8,  W));
    // - 9
    assert_eq!( 9_u64, layer.neighbour(9,  C).unwrap());
    assert_eq!( 1_u64, layer.neighbour(9,  N).unwrap());
    assert_eq!( 6_u64, layer.neighbour(9, NE).unwrap());
    assert_eq!( 5_u64, layer.neighbour(9, NW).unwrap());
    assert_eq!(11_u64, layer.neighbour(9,  S).unwrap());
    assert_eq!(10_u64, layer.neighbour(9, SE).unwrap());
    assert_eq!( 8_u64, layer.neighbour(9, SW).unwrap());
    assert_eq!(  None, layer.neighbour(9,  E));
    assert_eq!(  None, layer.neighbour(9,  W));
    // - 10
    assert_eq!(10_u64, layer.neighbour(10,  C).unwrap());
    assert_eq!( 2_u64, layer.neighbour(10,  N).unwrap());
    assert_eq!( 7_u64, layer.neighbour(10, NE).unwrap());
    assert_eq!( 6_u64, layer.neighbour(10, NW).unwrap());
    assert_eq!( 8_u64, layer.neighbour(10,  S).unwrap());
    assert_eq!(11_u64, layer.neighbour(10, SE).unwrap());
    assert_eq!( 9_u64, layer.neighbour(10, SW).unwrap());
    assert_eq!(  None, layer.neighbour(10,  E));
    assert_eq!(  None, layer.neighbour(10,  W));
    // - 11
    assert_eq!(11_u64, layer.neighbour(11,  C).unwrap());
    assert_eq!( 3_u64, layer.neighbour(11,  N).unwrap());
    assert_eq!( 4_u64, layer.neighbour(11, NE).unwrap());
    assert_eq!( 7_u64, layer.neighbour(11, NW).unwrap());
    assert_eq!( 9_u64, layer.neighbour(11,  S).unwrap());
    assert_eq!( 8_u64, layer.neighbour(11, SE).unwrap());
    assert_eq!(10_u64, layer.neighbour(11, SW).unwrap());
    assert_eq!(  None, layer.neighbour(11,  E));
    assert_eq!(  None, layer.neighbour(11,  W));
  }

  #[test]
  fn testok_neighbours_d0() {
    let layer = get_or_create(0);
    // North polar cap
    check_equals(layer.neighbours(0_u64, false), [1, 2, 3, 4, 5, 8]);
    check_equals(layer.neighbours(1_u64, false), [0, 2, 3, 5, 6, 9]);
    check_equals(layer.neighbours(2_u64, false), [0, 1, 3, 6, 7, 10]);
    check_equals(layer.neighbours(3_u64, false), [0, 1, 2, 4, 7, 11]);
    // Equatorial region
    check_equals(layer.neighbours(4_u64, false), [0, 3, 5, 7, 8, 11]);
    check_equals(layer.neighbours(5_u64, false), [0, 1, 4, 6, 8, 9]);
    check_equals(layer.neighbours(6_u64, false), [1, 2, 5, 7, 9, 10]);
    check_equals(layer.neighbours(7_u64, false), [2, 3, 4, 6, 10, 11]);
    //  South polar cap
    check_equals(layer.neighbours(8_u64, false), [0, 4, 5, 9, 10, 11]);
    check_equals(layer.neighbours(9_u64, false), [1, 5, 6, 8, 10, 11]);
    check_equals(layer.neighbours(10_u64, false), [2, 6, 7, 8, 9, 11]);
    check_equals(layer.neighbours(11_u64, false), [3, 4, 7, 8, 9, 10]);
  }
  
  #[test]
  fn test_ok_neighbours_t1() {
    let depth = 2_u8;
    let hash = 130;
    
    let layer = get_or_create(depth);
    check_equals_all(layer.neighbours(hash, true), [128, 129, 130, 131, 136, 137, 176, 177, 180]);
  }
  
  fn check_equals(map: MainWindMap<u64>, array: [u64; 6]) {
    array.iter().zip(map.sorted_values_vec().iter())
      .for_each(|(h1, h2)| assert_eq!(h1, h2) );
  }

  fn check_equals_all(map: MainWindMap<u64>, array: [u64; 9]) {
    array.iter().zip(map.sorted_values_vec().iter())
      .for_each(|(h1, h2)| assert_eq!(h1, h2) );
  }
  
  /*#[test]
  fn testok_cone_flat() {
    let res = cone_overlap_flat(3, 13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
    // println!("@@@@@@@@@@@@@@ {:?}", &res);
  }*/
  
  
  #[test]
  fn testok_external_edge_struct() {
    let depth = 1;
    let hash = 10;
    let delta_depth = 2;
    // draw moc 3/117,138,139,142,143,437,439,445,447,176, 178, 184, 186,85, 87, 93, 95,415,154
    // let e = external_edge_struct(depth, hash, delta_depth);
    // println!("{:?}", &e);
    let actual_res = external_edge_sorted(depth, hash, delta_depth);
    let expected_res: [u64; 19] = [85, 87, 93, 95, 117, 138, 139, 142, 143, 154, 176, 178, 184, 186, 415, 437, 439, 445, 447];
    // println!("{:?}", &actual_res);
    for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
      assert_eq!(h1, h2);
    }
    assert_eq!(expected_res.len(), actual_res.len());
  }

  #[test]
  fn testok_external_edge_struct_v2() {
    let depth = 1;
    let hash = 11;
    let delta_depth = 2;
    // draw moc 3/63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255
    let actual_res = external_edge_sorted(depth, hash, delta_depth);
    let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
    // println!("{:?}", &actual_res);
    for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
      assert_eq!(h1, h2);
    }
    assert_eq!(expected_res.len(), actual_res.len());
  }
 
  #[test]
  fn testok_external_edge_struct_v3() {
    let depth = 0;
    let hash = 0;
    let delta_depth = 2;
    // draw moc 3/63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255
    let actual_res = external_edge_sorted(depth, hash, delta_depth);
    //let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
    println!("{:?}", &actual_res);
    /*for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
      assert_eq!(h1, h2);
    }
    assert_eq!(expected_res.len(), actual_res.len());*/
  }
 
  #[test]
  fn testok_cone_approx_bmoc() {
    // let res = cone_overlap_approx(5, 0.01, 0.02, 0.05);
    // let res = cone_overlap_approx(6, 160.771389_f64.to_radians(), 64.3813_f64.to_radians(), 0.8962_f64.to_radians());
    let actual_res = cone_coverage_approx(3, 13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
    let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    /*for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
  }

  #[test]
  fn testok_cone_approx_custom_bmoc_2() {
    // let res = cone_overlap_approx(5, 0.01, 0.02, 0.05);
    // let res = cone_overlap_approx(6, 160.771389_f64.to_radians(), 64.3813_f64.to_radians(), 0.8962_f64.to_radians());
    let actual_res = cone_coverage_approx_custom(3, 2,36.80105218_f64.to_radians(), 56.78028536_f64.to_radians(), 14.93_f64.to_radians());
    let expected_res: [u64; 22] = [26, 27, 30, 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57, 58, 59, 60];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    /*println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
  }
  
  #[test]
  fn testok_cone_approx_custom_bmoc() {
    let actual_res = cone_coverage_approx_custom(3, 2,13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
    let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }
  
  #[test]
  fn testok_cone_approx_custom_bmoc_dbg() {
    let actual_res = cone_coverage_approx_custom(2, 1,20_f64.to_radians(), 0.0_f64.to_radians(), 50.0_f64.to_radians());
    let expected_res: [u64; 50] = [0, 1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 49, 52, 53, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 82, 88, 89, 90, 91, 94, 131, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 181, 183, 189];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }
  
  
  #[test]
  fn testok_elliptical_cone() {
    let lon = 36.80105218_f64.to_radians();
    let lat = 56.78028536_f64.to_radians();
    let a = 14.93_f64.to_radians();
    let b = 4.93_f64.to_radians();
    let pa = 75.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(3, lon, lat, a, b, pa);
    let expected_res: [u64; 16] = [27, 30, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57];
    to_aladin_moc(&actual_res);
    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    to_aladin_moc(&actual_res);
    println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) { 
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_2() {
    let lon = 0.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 30.0_f64.to_radians();
    let b = 5.0_f64.to_radians();
    let pa = 30.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(4, lon, lat, a, b, pa);
    let expected_res: [u64; 112] = [130, 136, 138, 1035, 1038, 1039, 1050, 1051, 1056, 1057, 1058,
      1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 
      1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1120,
      1122, 1123, 1126, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1144, 1146, 1147, 1150,
      1153, 1156, 1157, 1159, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1177 ,1180, 1181,
      1183, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225 ,1226, 1227 ,1228 ,1229,
      1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244,
      1245, 1246, 1247, 1252, 1253, 1264, 1265, 1268, 2933, 2935, 2941];
    /*76 [130, 136, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064,
      1065, 1067, 1068, 1069, 1070, 1071, 1072, 1074, 1075, 1078, 1079, 1080, 1081, 1082, 1083, 
      1084, 1085, 1086, 1087, 1128, 1129, 1130, 1131, 1132, 1134, 1135, 1146, 1157, 1168, 1169,
      1171, 1172, 1173, 1174, 1175, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 
      1228, 1229, 1231, 1232, 1233, 1234, 1235, 1236, 1238, 1239, 1240, 1241, 1242, 1243, 1244,
      1245, 1246, 1247, 2935, 2941];*/
    // to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }
  
  #[test]
  fn testok_elliptical_cone_3() {
    let lon = 0.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 50.0_f64.to_radians();
    let b = 5.0_f64.to_radians();
    let pa = 30.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(3, lon, lat, a, b, pa);
    let expected_res: [u64; 70] = [10, 11, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 256, 257, 258,
      259, 260, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 274, 280, 281, 282, 283, 284, 
      286, 287, 288, 289, 291, 292, 293, 294, 295, 301, 304, 305, 306, 307, 308, 309, 310, 311, 312,
      313, 315, 316, 317, 318, 319, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 756, 757];
    /*40 [32, 33, 34, 35, 36, 38, 40, 258, 264, 265, 266, 267, 268, 269,
      270, 271, 280, 282, 283, 286, 289, 292, 293, 295, 304, 305, 306, 307, 308, 309, 310, 311, 317,
      727, 729, 731, 732, 733, 734, 735];*/
    // to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_4() {
    let lon = 0.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 50.0_f64.to_radians();
    let b = 5.0_f64.to_radians();
    let pa = 20.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(3, lon, lat, a, b, pa);
    let expected_res: [u64; 62] = [32, 34, 35, 38, 40, 41, 42, 43, 44, 256, 257, 258, 259, 260,
      262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 274, 280, 281, 282, 283, 286, 287, 288, 289,
      292, 293, 294, 295, 301, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 315, 316 ,317, 318,
      319, 723, 724, 725, 726, 727, 729, 732, 733, 735];
    /*40 [34, 35, 38, 40, 41, 44, 258, 259, 262, 264, 265, 266, 267, 268, 
      269, 270, 271, 280, 282, 283, 292, 293, 295, 304, 305, 306, 307, 308, 309, 310, 311, 313, 
      316, 317, 723, 726, 727, 729, 732, 733];*/
    // to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_5() {
    let lon = 0.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 42.0_f64.to_radians();
    let b = 10.0_f64.to_radians();
    let pa = 20.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(2, lon, lat, a, b, pa);
    let expected_res: [u64; 26] = [8, 9, 10, 11, 52, 53, 64, 65, 66, 67, 68, 70, 71, 72, 73, 75, 
      76, 77, 78, 79, 138, 139, 180, 181, 182, 183];
    /*20 [8, 9, 10, 64, 65, 66, 67, 68, 70, 71, 72, 73, 75, 76, 77, 78, 79,
      181 ,182 ,183];*/
    to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_6() {
    let lon = 360.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 50.0_f64.to_radians();
    let b = 5.0_f64.to_radians();
    let pa = 0.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(2, lon, lat, a, b, pa);
    let expected_res: [u64; 18] = [10, 11, 53, 55, 64, 65, 66, 67, 70, 73, 76, 77, 78, 79, 136, 
      138, 180, 181];
    // to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_7() {
    let lon = 0.0_f64.to_radians();
    let lat = 0.0_f64.to_radians();
    let a = 50.0_f64.to_radians();
    let b = 5.0_f64.to_radians();
    let pa = 90.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(3, lon, lat, a, b, pa);
    let expected_res: [u64; 40] = [0, 192, 269, 270, 271, 273, 274, 275, 276, 277, 278, 279, 280, 
      281, 282, 283, 284, 285, 286, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301,
      302, 304, 305, 306, 362, 469, 575, 767];
    // to_aladin_moc(&actual_res);
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());
  }

  #[test]
  fn testok_elliptical_cone_8() {
    let lon = 50.0_f64.to_radians();
    let lat = 50.0_f64.to_radians();
    let a = 60.0_f64.to_radians();
    let b = 2.0_f64.to_radians();
    let pa = 35.0_f64.to_radians();
    let actual_res = elliptical_cone_coverage(6, lon, lat, a, b, pa);
    /*let expected_res: [u64; 40] = [0, 192, 269, 270, 271, 273, 274, 275, 276, 277, 278, 279, 280,
      281, 282, 283, 284, 285, 286, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301,
      302, 304, 305, 306, 362, 469, 575, 767];*/
    to_aladin_moc(&actual_res);
    /*for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());*/
  }
  
  fn to_aladin_moc(bmoc: &BMOC) {
    print!("draw moc {}/", bmoc.get_depth_max());
    for cell in bmoc.flat_iter() {
      print!("{},", cell);
    }
  }
  
  #[test]
  fn testok_polygone_approx() {
    let actual_res = polygon_coverage(3, &[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)], false);
    /*let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }*/
    println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
  }

  #[test]
  fn testok_polygone_exact_npc() {
    // Aladin: draw polygon(65.11781779000003, 85.012424, 89.70533626000001, 87.06130188, 60.23667431000001, 85.609882)
    let depth = 6;
    let mut vertices = [(65.11781779000003, 85.012424), (89.70533626000001, 87.06130188), (60.23667431000001, 85.609882)];
    let expected_res_approx: [u64; 2] = [4062, 4063];
    let expected_res_exact: [u64; 4] = [4062, 4063, 4084, 4085];

    to_radians(&mut vertices);
    
    let actual_res_approx = polygon_coverage(depth, &vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }    
    
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_npc_2() {
    // Aladin: draw polygon(359.70533626,+87.06130188, 330.23667431,+85.60988200, 335.11781779,+85.01242400)
    let depth = 6;
    let mut vertices = [(359.70533626, 87.06130188), (330.23667431, 85.60988200), (335.11781779, 85.01242400)];
    let expected_res_approx: [u64; 2] = [16350, 16351];
    let expected_res_exact: [u64; 4] = [16350, 16351, 16372, 16373];
    
    to_radians(&mut vertices);

    let actual_res_approx = polygon_coverage(depth, &vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }

    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_npc_3() {
    // Aladin:  draw polygon(224.86211710,+78.10924662, 176.91129363 +83.92878811, 135.81578643,+78.24840426, 200.73574863,+73.58038790)
    let depth = 3;
    let mut vertices = [(224.86211710, 78.10924662), (176.91129363, 83.92878811), (135.81578643, 78.24840426), (200.73574863, 73.58038790)];
    let expected_res_approx: [u64; 5] = [119, 125, 187, 188, 190];
    let expected_res_exact: [u64; 7] =  [119, 125, 127, 187, 188, 190, 191];

    // Pb pour cell 127:
    // from (180, +83) --> (224, +78)=> 191
    // from (135, +78) --> (176, +83) => 127
    
    to_radians(&mut vertices);

    let actual_res_approx = polygon_coverage(depth, &vertices, false);
    println!("draw moc 3/ {:?}", actual_res_approx.to_flat_array());
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }
    
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    println!("draw moc 3/ {:?}", actual_res_exact.to_flat_array());
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_spc() {
    // Aladin: draw polygon(359.70533626,-87.06130188, 330.23667431,-85.60988200, 335.11781779,-85.01242400)
    let depth = 6;
    let mut vertices = [(359.70533626, -87.06130188), (330.23667431, -85.60988200), (335.11781779, -85.01242400)];
    let expected_res_approx: [u64; 2] = [45072, 45074];
    let expected_res_exact: [u64; 4] = [45061, 45063, 45072, 45074];

    to_radians(&mut vertices);

    let actual_res_approx = polygon_coverage(depth, &vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }

    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_eqr() {
    // In Aladin: draw polygon(180.08758393,-41.43289179, 191.00310758,-29.99207687, 181.59160475,-34.21976170)
    let depth = 3;
    let mut vertices = [(180.08758393,-41.43289179), (191.00310758,-29.99207687), (181.59160475,-34.219761700)];
    let expected_res_approx: [u64; 2] = [384, 385];
    let expected_res_exact: [u64; 4] = [384, 385, 682, 683];

    to_radians(&mut vertices);

    let actual_res_approx = polygon_coverage(depth, &vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }

    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_2() {
    // In Aladin: draw polygon(174.75937396073138, -49.16744206799886, 185.24062603926856, -49.16744206799887, 184.63292896369916, -42.32049830486584, 175.3670710363009, -42.32049830486584)
    let depth = 10;
    
    let mut vertices = [(174.75937396073138, -49.16744206799886), 
      (185.24062603926856, -49.16744206799887), 
      (184.63292896369916, -42.32049830486584),
      (175.3670710363009, -42.32049830486584)];
    //let expected_res_approx: [u64; 2] = [384, 385];
    //let expected_res_exact: [u64; 4] = [384, 385, 682, 683];

    to_radians(&mut vertices);

    /*let actual_res_approx = polygon_coverage(depth, &vertices, false);
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }*/

    let actual_res_exact = polygon_coverage(depth, &vertices, false);

    println!("@@@@@ FLAT VIEW");
    for cell in actual_res_exact.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    
    assert!(actual_res_exact.deep_size() > 0);
    
    /*assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }*/
  }
  
  fn to_radians(lonlats: &mut [(f64, f64)]) {
    for (lon, lat) in lonlats.iter_mut() {
      *lon = lon.to_radians();
      *lat = lat.to_radians();
    }
  }
  
  /*fn test_polygone(depth: u8, lonlats: &[(f64, f64)]) {
    let bmoc = polygon_coverage(depth, lonlats, true);
    
  }*/
  
  #[test]
  fn testok_bmoc_not() {
    let actual_res = cone_coverage_approx_custom(3, 4, 36.80105218_f64.to_radians(), 56.78028536_f64.to_radians(), 14.93_f64.to_radians());
    println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    println!("@@@@@ HIERARCH VIEW, NOT");
    let complement: BMOC = actual_res.not();
    for cell in complement.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    let org = complement.not();
    assert!(actual_res.equals(&org));
    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
  }
  
  #[test]
  fn test_prec_1() {
    let lon_deg = 179.99999999999997_f64;
    let lat_deg = 41.813964843754924_f64;
    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get_or_create(0);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    let layer_0 = get_or_create(0);
    assert_eq!(1, layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_prec_2() {
    let lon_deg = 359.99999999999994_f64;
    let lat_deg = 41.81031502783791_f64;
    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_1 = get_or_create(1);
    layer_1.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_1.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);*/
    let layer_1 = get_or_create(1);
    assert_eq!(13, layer_1.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_prec_3() {
    let lon_deg = 359.99999999999994_f64; //359.99999999999994_f64;
    let lat_deg = 41.81031489577861_f64; //41.81031489577857_f64;

    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get_or_create(6);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));*/
    
    let layer_6 = get_or_create(6);
    assert_eq!(13653, layer_6.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_prec_4() {
    let lon_deg = 292.49999999999994_f64; //359.99999999999994_f64;
    let lat_deg = 41.810314895778546_f64; //41.81031489577857_f64;

    let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    println!("proj_x: {:.17}, proj_y: {:.17}", xy.0, xy.1);
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get_or_create(0);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));

    let layer_2 = get_or_create(2);
    assert_eq!(56, layer_2.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_ring() {
    for depth in 0..10 {
      let layer = get_or_create(depth);
      for h in 0..layer.n_hash {
        assert_eq!(layer.from_ring(layer.to_ring(h)), h);
      }
    }
  }

  #[test]
  fn test_bilinear_interpolation() {
    let lon_deg = 89.18473162_f64; // 322.99297784_f64;// 324.8778822_f64
    let lat_deg = -28.04159707_f64;// 39.9302924_f64;// -41.08635508_f64
    let res = bilinear_interpolation(1, lon_deg.to_radians(), lat_deg.to_radians());
    // println!("{:?}", res);
    assert_eq!(res, [
      (20, 0.0), 
      (38, 0.1661686383097217), 
      (33, 0.2024027885319438), 
      (20, 0.6314285731583344)
    ]);
  }
  
  #[test]
  fn test_bilinear_interpolation_2() {
    let lon_deg = 83.633478_f64;
    let lat_deg = 22.015110_f64;
    let res = bilinear_interpolation(18, lon_deg.to_radians(), lat_deg.to_radians());
    // println!("{:?}", res);
    assert_eq!(res, [
      (405766747916, 0.5757471135241182), 
      (405766747917, 0.3604806280107034), 
      (405766747918, 0.039217694696856834), 
      (405766747919, 0.024554563768321474)
    ]);
  }

  #[test]
  fn test_gen_file() -> std::io::Result<()> {
    use std::fs::File;
    use std::io::prelude::*;
    
    let depth = 8;
    let layer = get_or_create(depth);
    let n_cells = layer.n_hash();
    
    let mut s = String::with_capacity(8 * 1024);
    s.push_str("i,ra,dec\n");
    for h in 0..n_cells {
      let (lon, lat) = layer.center(h);
      s.push_str(&format!("{},{},{}\n", &h, &lon.to_degrees(), &lat.to_degrees()));
    }
    
    let mut file = std::fs::File::create(format!("hpx.{}.csv", depth))?;
    file.write_all(s.as_bytes())?;
    Ok(())
  }


  #[test]
  fn test_to_uniq() {
    for depth in 0..8 {
      for idx in 0..n_hash(depth) {
        assert_eq!((depth, idx), from_uniq(to_uniq(depth, idx)));
      }
    }
  }

  #[test]
  fn test_to_uniq_ivoa() {
    for depth in 0..8 {
      for idx in 0..n_hash(depth) {
        assert_eq!((depth, idx), from_uniq_ivoa(to_uniq_ivoa(depth, idx)));
      }
    }
  }

  #[test]
  fn test_xpm1_and_q () {
    let a = Layer::d0h_lh_in_d0c(std::f64::NAN, 0.0);
    eprintln!("{:?}", &a);
    let a = Layer::d0h_lh_in_d0c(std::f64::INFINITY, 0.0);
    eprintln!("{:?}", &a);
  }

}
