use std::sync::Once;
use super::compass_point::*;
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

pub fn to_range(hash: u64, delta_depth: u8) -> std::ops::Range<u64> {
  let twice_delta_depth = delta_depth << 1;
  (hash << twice_delta_depth)..((hash + 1) << twice_delta_depth)
}

/// Conveniency function simply calling the [hash](struct.Layer.html#method.hash) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn hash(depth: u8, lon: f64, lat: f64) -> u64 {
  get_or_create(depth).hash(lon, lat)
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

/// Conveniency function simply calling the [neighbours](struct.Layer.html#method.neighbours) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn neighbours(depth: u8, hash: u64, include_center: bool) -> MainWindMap<u64> {
  get_or_create(depth).neighbours(hash, include_center)
}

/// Conveniency function simply calling the [cone_overlap_approx](struct.Layer.html#method.cone_overlap_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_overlap_approx(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get_or_create(depth).cone_overlap_approx(cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_overlap_flat](struct.Layer.html#method.cone_overlap_approx) method
/// of the [Layer] of the given *depth* and retrieving a flat array.
#[inline]
pub fn cone_overlap_approx_flat(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> Box<[u64]> {
  get_or_create(depth).cone_overlap_approx(cone_lon, cone_lat, cone_radius).to_flat_array()
}

/// Conveniency function simply calling the [cone_overlap_approx_custom](struct.Layer.html#method.cone_overlap_approx_custom) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_overlap_approx_custom(depth: u8, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get_or_create(depth).cone_overlap_approx_custom(delta_depth, cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [polygon_overlap_approx](struct.Layer.html#method.polygon_overlap_approx) method
/// of the [Layer] of the given *depth*.
pub fn polygon_overlap_approx(depth: u8, vertices: &[(f64, f64)]) -> BMOC {
  get_or_create(depth).polygon_overlap_approx(vertices)
}


pub mod bmoc;
mod zordercurve;

use self::zordercurve::{ZOrderCurve, get_zoc};
use self::bmoc::*;
use super::sph_geom::coo3d::*;
use super::sph_geom::{Polygon};
use super::sph_geom::cone::{Cone};
use super::sph_geom::coo3d::LonLat;
use super::{proj};
use super::compass_point::{Cardinal, CardinalSet, CardinalMap};
use super::compass_point::MainWind::{S, SE, E, SW, C, NE, W, NW, N};

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
  z_order_curve: &'static ZOrderCurve,
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
    let mut xy = proj(lon, lat);
    xy.0 = ensures_x_is_positive(xy.0);
    self.shift_rotate_scale(&mut xy);
    let mut ij = discretize(xy);
    let ij_d0c = self.base_cell_coos(&ij);
    let d0h_bits = self.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    self.to_coos_in_base_cell(&mut ij);
    self.build_hash(d0h_bits, ij.0 as u32, ij.1 as u32)
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
  fn depth0_bits(&self, i: u8, j: u8/*, ij: &mut (u64, u64)*, xy: (f64, f64), lon: f64, lat: f64*/) -> u64 {
    // self.base_hash_bits_lupt[i as usize][j as usize]
    // Useful for quick tests: https://play.rust-lang.org
    let k = 5_i8 - (i + j) as i8;
    // The two branches -1 and -2 are extremely rare (north pole, NPC cells upper border),
    // so few risks of branch miss-prediction.
    match k {
      0 ... 2 => (((k << 2) + ( ((i as i8) + ((k - 1) >> 7)) & 3_i8)) as u64) << self.twice_depth,
      -1 => ((((i - 1u8) & 3u8) as u64) << self.twice_depth) | self.y_mask,
      -2 => (((i - 2u8) as u64) << self.twice_depth) | self.xy_mask,
      /* I wrote this code when I forget to handle negative longitude (and I found catalogues containing
         small negative longitudes).
      3 => { // rare case
        let d0 = (xy.0 - ij.0 as f64).abs();
        let d1 = (xy.1 - ij.1 as f64).abs();
        if d0 < d1 {
          ij.0 += 1;
          return self.depth0_bits(i + 1_u8, j, ij, xy, lon, lat);
        } else {
          ij.1 += 1;
          return self.depth0_bits(i, j + 1_u8, ij, xy, lon, lat);
        }
      },
      4 => { // rare case
        ij.0 += 1;
        ij.1 += 1;
       return self.depth0_bits(i + 1_u8, j + 1_u8, ij, xy, lon, lat);
      },
      _ => panic!("Algorithm error: case k = {} not supported! depth: {}, lon: {}, lat: {}, x: {}, y: {}", 
                  k, self.depth, lon, lat, xy.0, xy.1),*/
      _ => panic!("Algorithm error: case k = {} not supported!"),
    }
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
      result_map.get(C);
    }
    let h_bits: HashBits = self.pull_bits_appart(hash);
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      self.edge_cell_neighbours(hash, &mut result_map);
    } else {
      self.inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, &mut result_map);
    }
    result_map
  }

  #[inline]
  fn is_in_base_cell_border(&self, i_in_base_cell_bits: u64, j_in_base_cell_bits: u64) -> bool {
    0_u64 == i_in_base_cell_bits || i_in_base_cell_bits == self.x_mask
      || 0_u64 == j_in_base_cell_bits || j_in_base_cell_bits == self.y_mask
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
      match d0h >> 2 {
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

  // Center of a cell, on the projection pace (i.e. in the Euclidean plane).
  fn center_of_projected_cell(&self, hash: u64) -> (f64, f64) {
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
  fn decode_hash(&self, hash: u64) -> HashParts {
    let ij: u64 = self.z_order_curve.h2ij(hash & self.xy_mask);
    HashParts {
      d0h: (hash >> self.twice_depth) as u8,
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
  fn h_and_shs_to_lower_h(&self, deeper_depth: u8) 
    -> impl Fn((u64, f64)) -> (u64) {
    let twice_depth_diff = (deeper_depth - self.depth) << 1;
    move |(h, _)| h >> twice_depth_diff
  }

  /// Returns a hierarchical view of the list of cells overlapped by the given cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
  /// The algorithm is fast but approximated: it may return false positive, 
  /// i.e. cells which are near from the cone but do not overlap it.
  /// To control the approximation, see the method 
  /// [cone_overlap_approx_custom](#method.cone_overlap_approx_custom)
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
  /// let actual_res = nested3.cone_overlap_approx(lon, lat, radius);
  /// let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///      assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_overlap_approx(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    self.cone_overlap_approx_internal(cone_lon, cone_lat, cone_radius).to_bmoc_packing()
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
  /// let actual_res = nested3.cone_overlap_approx_custom(2, lon, lat, radius);
  /// let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_overlap_approx_custom(&self, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    get_or_create(self.depth + delta_depth)
      .cone_overlap_approx_internal(cone_lon, cone_lat, cone_radius)
      .to_lower_depth_bmoc_packing(self.depth)
  }
  
  fn cone_overlap_approx_internal(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOCBuilderUnsafe {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 12);
      bmoc_builder.push_all(0_u8, 0_u64, 12_u64, true);
      return bmoc_builder; //.to_bmoc();
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
        self.cone_overlap_approx_recur(0, h,
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
        self.cone_overlap_approx_recur(depth_start, root_hash,
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                                &minmax_array, 0, &mut bmoc_builder);
      }
      return bmoc_builder; //.to_bmoc_packing();
    }
  }
  fn cone_overlap_approx_recur<F>(&self, depth: u8, hash: u64, shs_computer: &F, shs_minmax: &[MinMax],
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
        self.cone_overlap_approx_recur(depth, hash, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_overlap_approx_recur(depth, hash | 1_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_overlap_approx_recur(depth, hash | 2_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_overlap_approx_recur(depth, hash | 3_u64, shs_computer, shs_minmax, recur_depth, bmoc_builder);
      }
    }
  }
  /*
  /// Returns the list of cells overlapping the given cone
  /// REPLACE THIS BY cone_overlap_approx AND CALLING A METHOD ON THE RESULTING BMOC
  pub fn cone_overlap_flat(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> Box<[u64]> {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      let v: Vec<u64> = (0..self.n_hash).collect();
      return v.into_boxed_slice();
    }
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius) {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, distances);
      let mut res: Vec<u64> = Vec::with_capacity(self.ncell_in_cone_upper_bound(cone_radius));
      for h in 0..12 {
        self.cone_overlap_recur(0, h, 
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat), 
                                &minmax_array, 0, &mut res);
      }
      return res.into_boxed_slice();
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
      return neigs.into_boxed_slice();
    } else {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, distances);
      let root_layer = get_or_create(depth_start);
      let root_center_hash = root_layer.hash(cone_lon, cone_lat);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut res: Vec<u64> = Vec::with_capacity(self.ncell_in_cone_upper_bound(cone_radius));
      for &root_hash in neigs.sorted_values().into_iter() {
        self.cone_overlap_recur(depth_start, root_hash, 
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat), 
                                &minmax_array, 0, &mut res);
      }
      return res.into_boxed_slice();
    }
  }
  
  fn cone_overlap_recur<F>(&self, depth: u8, hash: u64, shs_computer: &F, shs_minmax: &[MinMax], 
                           recur_depth: u8, result: &mut Vec<u64>)
  where F: Fn((f64, f64)) -> f64  {
    let shs = shs_computer(get_or_create(depth).center(hash));
    let MinMax{min, max} = shs_minmax[recur_depth as usize];
    if shs <= min {
      let twice_depth_diff = (self.depth - depth) << 1;
      let h_min = hash << twice_depth_diff;
      let h_max = (hash + 1) << twice_depth_diff;
      // result.extend((h_min..h_max).into_iter()); // do not work in webassembly??
      for h in (h_min..h_max).into_iter() {
        result.push(h);
      }
    } else if shs <= max {
      if depth == self.depth {
        result.push(hash);
      } else {
        let deeper_hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.cone_overlap_recur(depth, deeper_hash, shs_computer, shs_minmax, recur_depth, result);
        self.cone_overlap_recur(depth, deeper_hash | 1_u64, shs_computer, shs_minmax, recur_depth, result);
        self.cone_overlap_recur(depth, deeper_hash | 2_u64, shs_computer, shs_minmax, recur_depth, result);
        self.cone_overlap_recur(depth, deeper_hash | 3_u64, shs_computer, shs_minmax, recur_depth, result);
      }
    }
  }*/
  
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
    4_usize * (1_usize + (self.nside as f64 * TWICE_SQRT_3 * cone_radius + 0.99f64) as usize)
  }


  /// Returns a hierarchical view of the list of cells overlapped by the given polygon.
  /// Self intersecting polygons are supported.
  /// The BMOC also tells if the cell if fully or partially overlapped by the polygon.
  /// 
  /// If you want the complementary solution, apply the NOT operator on the BMOC (to be implemented).
  /// 
  /// Here we make the following approximation when testing the intersection between a polygon segment
  /// and an HEALPix cell edge: we consider that each edge of the HEALPix cell is on a great-circle arc.*
  /// We plan to provide an exact solution by first testing for each polygon segment if it
  /// contains a 'special point', like for the exact cone solution (implemented in Java but not yet in Rust).
  ///
  /// 
  /// # Input
  /// - `vertices` the list of vertices (in a slice) coordinates, in radians
  ///              `[(lon, lat), (lon, lat), ..., (lon, lat)]`
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
  /// let actual_res =nested3.polygon_overlap_approx(&[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)]);
  /// let expected_res: [u64; 8] = [304, 305, 306, 307, 308, 310, 313, 316];
  /// 
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn polygon_overlap_approx(&self, vertices: &[(f64, f64)]) -> BMOC {
    let poly = Polygon::new(
      vertices.iter().map(|(lon, lat)| LonLat { lon: *lon, lat: *lat} )
        .collect::<Vec<LonLat>>().into_boxed_slice()
    );
    let bounding_cone: Cone = Cone::bounding_cone(poly.vertices());
    let mut depth_start = 0;
    let neigs: Vec<u64> = if !has_best_starting_depth(bounding_cone.radius()) {
      (0..12).collect()
    } else {
      depth_start = best_starting_depth(bounding_cone.radius());
      let root_layer = get_or_create(depth_start);
      let LonLat{lon, lat} = bounding_cone.center().lonlat();
      let center_hash = root_layer.hash(lon, lat);
      let mut neigs: Vec<u64> = root_layer.neighbours(center_hash, true).values_vec();
      neigs.sort_unstable();
      neigs
    };
    // Compute and sort the list of cell containing at least one polygon vertex
    let mut sorted_poly_vertices_hash = self.hashs(poly.vertices());
    sorted_poly_vertices_hash.sort_unstable();
    // Build the list (removing duplicated) for all deltaDepth?
    let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 10_000_usize);
    for root_hash in neigs { 
      self.polygon_overlap_recur(&mut bmoc_builder, depth_start, root_hash, &poly, &sorted_poly_vertices_hash);
    }
    bmoc_builder.to_bmoc()
    // return HealpixNestedBMOC.createUnsafe(this.depthMax, moc.a, moc.size());
    // BMOC::create_unsafe(self.depth, Vec::with_capacity(0).into_boxed_slice())
  }

  fn polygon_overlap_recur(&self, moc_builder: &mut BMOCBuilderUnsafe, depth: u8, hash: u64,
  poly: &Polygon, sorted_poly_vertices_hash: &[u64]) {
    if is_in_list(depth, hash, self.depth, sorted_poly_vertices_hash) {
      if depth == self.depth {
        moc_builder.push(depth, hash, false);
      } else {
        let hash = hash << 2;
        let depth = depth + 1;
        self.polygon_overlap_recur(moc_builder, depth, hash    , poly, sorted_poly_vertices_hash);
        self.polygon_overlap_recur(moc_builder, depth, hash | 1, poly, sorted_poly_vertices_hash);
        self.polygon_overlap_recur(moc_builder, depth, hash | 2, poly, sorted_poly_vertices_hash);
        self.polygon_overlap_recur(moc_builder, depth, hash | 3, poly, sorted_poly_vertices_hash);
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
          self.polygon_overlap_recur(moc_builder, depth, hash    , poly, sorted_poly_vertices_hash);
          self.polygon_overlap_recur(moc_builder, depth, hash | 1, poly, sorted_poly_vertices_hash);
          self.polygon_overlap_recur(moc_builder, depth, hash | 2, poly, sorted_poly_vertices_hash);
          self.polygon_overlap_recur(moc_builder, depth, hash | 3, poly, sorted_poly_vertices_hash);
        }
      }
    }
  }

  fn hashs<T: LonLatT>(&self, poly_vertices: &[T]) -> Box<[u64]> {
    poly_vertices.iter().map(|coo| self.hash(coo.lon(), coo.lat()))
      .collect::<Vec<u64>>().into_boxed_slice()
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
  j: u64,   // in the base cell, z-order curve coordinate along the x-axis bits
}

#[inline]
fn discretize(xy: (f64, f64)) -> (u64, u64) {
  (xy.0 as u64, xy.1 as u64)
}

#[inline]
fn rotate45_scale2(i_in_d0h: u32, j_in_d0h: u32) -> (i32, i32) {
  (i_in_d0h as i32 - j_in_d0h as i32, (i_in_d0h + j_in_d0h) as i32)
}

/// offset_x in [0, 7], odd for polar caps, even for equatorial region
/// offset_y in [-1, 1], -1 or 1 for polar caps, 0 for equatorial region
#[inline]
fn compute_base_cell_center_offsets_in_8x3_grid(d0h: u8) -> (u8, i8) { 
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
fn bits_2_hash(d0h_bits: u64, i_in_d0h_bits: u64, j_in_d0h_bits: u64) -> u64 {
  d0h_bits | i_in_d0h_bits | j_in_d0h_bits
}

#[inline]
fn ensures_x_is_positive(x: f64) -> f64 {
  if x < 0.0 { x + 8.0 } else { x }
}

/// x / 4
#[inline]
fn div4_quotient(x: u8) -> u8 {
  x >> 2
}

/// x modulo 4
#[inline]
fn div4_remainder(x: u8) -> u8 {
  x & 3
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
fn shs_lower_than(shs_max: f64)
                  -> impl FnMut(&(u64, f64)) -> (bool) {
  move |&(_, shs)| shs <= shs_max
}

/// The returned closure computes the 
#[inline]
fn shs_computer(cone_lon: f64, cone_lat: f64, cos_cone_lat: f64) -> impl Fn((f64, f64)) -> (f64) {
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
  
  #[test]
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
      assert_eq!(lupm[*i as usize][*j as usize], layer.depth0_bits(*i, *j/*, &mut (0_u64, 0_u64), (0.0, 0.0), 0.0, 0.0*/));
    }
  }
  
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
  fn testok_cone_approx_bmoc() {
    // let res = cone_overlap_approx(5, 0.01, 0.02, 0.05);
    // let res = cone_overlap_approx(6, 160.771389_f64.to_radians(), 64.3813_f64.to_radians(), 0.8962_f64.to_radians());
    let actual_res = cone_overlap_approx(3, 13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
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
    let actual_res = cone_overlap_approx_custom(3, 2,36.80105218_f64.to_radians(), 56.78028536_f64.to_radians(), 14.93_f64.to_radians());
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
    let actual_res = cone_overlap_approx_custom(3, 2,13.158329_f64.to_radians(), -72.80028_f64.to_radians(), 5.64323_f64.to_radians());
    let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    /*println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
  }

  #[test]
  fn testok_polygone_approx() {
    let actual_res = polygon_overlap_approx(3, &[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)]);
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
  fn testok_bmoc_not() {
    let actual_res = cone_overlap_approx_custom(3, 4, 36.80105218_f64.to_radians(), 56.78028536_f64.to_radians(), 14.93_f64.to_radians());
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
  
}