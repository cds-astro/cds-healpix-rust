use std::ops::Shr;
use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};

use super::compass_point::*;
use super::external_edge::*;
use super::*;

pub mod zordercurve;

use self::zordercurve::{ZOC, ZOrderCurve, get_zoc};

// We use an array here since operations on f64 are not yet stable for `const fn` :o/
// Run the following code on e.g. https://play.rust-lang.org/
// fn main() {
//   for d in 0..=29 {
//     let nside = (1 as u32) << d;
//     let time_nside = (d as u64) << 52;
//     let one_over_nside_v1 = 1_f64 / nside as f64;
//     let one_over_nside_v2 = f64::from_bits(1_f64.to_bits() - time_nside);
//     assert_eq!(one_over_nside_v1, one_over_nside_v2);
//     println!("  {},", one_over_nside_v2);
//   }
// }
const ONE_OVER_NSIDE: [f64; 30] = [
  1.0,
  0.5,
  0.25,
  0.125,
  0.0625,
  0.03125,
  0.015625,
  0.0078125,
  0.00390625,
  0.001953125,
  0.0009765625,
  0.00048828125,
  0.000244140625,
  0.0001220703125,
  0.00006103515625,
  0.000030517578125,
  0.0000152587890625,
  0.00000762939453125,
  0.000003814697265625,
  0.0000019073486328125,
  0.00000095367431640625,
  0.000000476837158203125,
  0.0000002384185791015625,
  0.00000011920928955078125,
  0.00000005960464477539063,
  0.000000029802322387695313,
  0.000000014901161193847656,
  0.000000007450580596923828,
  0.000000003725290298461914,
  0.000000001862645149230957,
];

/*
// Wait for stabilization: https://github.com/rust-lang/rust/issues/72447
const fn one_over_nside(depth: u8) -> f64 {
  let time_nside = (depth as u64) << 52;
  f64::from_bits(1_f64.to_bits() - time_nside)
}*/

/// Array storing pre-computed values for each of the 30 possible depth (from 0 to 29)
const LAYERS: [Layer; 30] =  [
  Layer::new(0, ZOC::EMPTY),
  Layer::new(1, ZOC::SMALL),
  Layer::new(2, ZOC::SMALL),
  Layer::new(3, ZOC::SMALL),
  Layer::new(4, ZOC::SMALL),
  Layer::new(5, ZOC::SMALL),
  Layer::new(6, ZOC::SMALL),
  Layer::new(7, ZOC::SMALL),
  Layer::new(8, ZOC::SMALL),
  Layer::new(9, ZOC::MEDIUM),
  Layer::new(10, ZOC::MEDIUM),
  Layer::new(11, ZOC::MEDIUM),
  Layer::new(12, ZOC::MEDIUM),
  Layer::new(13, ZOC::MEDIUM),
  Layer::new(14, ZOC::MEDIUM),
  Layer::new(15, ZOC::MEDIUM),
  Layer::new(16, ZOC::MEDIUM),
  Layer::new(17, ZOC::LARGE),
  Layer::new(18, ZOC::LARGE),
  Layer::new(19, ZOC::LARGE),
  Layer::new(20, ZOC::LARGE),
  Layer::new(21, ZOC::LARGE),
  Layer::new(22, ZOC::LARGE),
  Layer::new(23, ZOC::LARGE),
  Layer::new(24, ZOC::LARGE),
  Layer::new(25, ZOC::LARGE),
  Layer::new(26, ZOC::LARGE),
  Layer::new(27, ZOC::LARGE),
  Layer::new(28, ZOC::LARGE),
  Layer::new(29, ZOC::LARGE)
];

pub fn get(depth: u8) -> &'static Layer {
  &LAYERS[depth as usize]
}


pub mod moc;

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
/// use cdshealpix::nested::{get, Layer};
/// let l0 = get(0);
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

#[allow(dead_code)]
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
  get(depth).n_hash
}

/// Conveniency function simply calling the [hash](struct.Layer.html#method.hash) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn hash(depth: u8, lon: f64, lat: f64) -> u64 {
  get(depth).hash(lon, lat)
}

/// Conveniency function simply calling the [hash](struct.Layer.html#method.hash_with_dxdy) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn hash_with_dxdy(depth: u8, lon: f64, lat: f64) -> (u64, f64, f64) {
  get(depth).hash_with_dxdy(lon, lat)
}

/// Conveniency function simply calling the [center](struct.Layer.html#method.sph_coo) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn sph_coo(depth: u8, hash: u64, dx: f64, dy: f64) -> (f64, f64) { 
  get(depth).sph_coo(hash, dx, dy)
}

/// Conveniency function simply calling the [center](struct.Layer.html#method.center) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn center(depth: u8, hash: u64) -> (f64, f64) {
  get(depth).center(hash)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.vertices) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn vertices(depth: u8, hash: u64) -> [(f64, f64); 4] {
  get(depth).vertices(hash)
}


/// Conveniency function simply calling the [vertices](struct.Layer.html#method.path_along_cell_side) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn path_along_cell_side(depth: u8, hash: u64, from_vertex:  &Cardinal, to_vertex: &Cardinal, include_to_vertex: bool, n_segments: u32) -> Box<[(f64, f64)]> {
  get(depth).path_along_cell_side(hash, from_vertex, to_vertex, include_to_vertex, n_segments)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.path_along_cell_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn path_along_cell_edge(depth: u8, hash: u64, starting_vertex: &Cardinal, clockwise_direction: bool, n_segments_by_side: u32) -> Box<[(f64, f64)]> {
  get(depth).path_along_cell_edge(hash, starting_vertex, clockwise_direction, n_segments_by_side)
}

/// Conveniency function simply calling the [vertices](struct.Layer.html#method.grid) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn grid(depth: u8, hash: u64, n_segments_by_side: u16) -> Box<[(f64, f64)]> {
  get(depth).grid(hash, n_segments_by_side)
}

/// Conveniency function simply calling the [neighbours](struct.Layer.html#method.neighbours) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn neighbours(depth: u8, hash: u64, include_center: bool) -> MainWindMap<u64> {
  get(depth).neighbours(hash, include_center)
}

/// Conveniency function simply calling the [append_bulk_neighbours](struct.Layer.html#method.append_bulk_neighbours) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn append_bulk_neighbours(depth: u8, hash: u64, dest: &mut Vec<u64>) {
  get(depth).append_bulk_neighbours(hash, dest);
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
pub fn internal_edge_sorted(depth: u8, hash: u64, delta_depth: u8) -> Vec<u64> {
  assert!(depth + delta_depth < DEPTH_MAX);
  Layer::internal_edge_sorted(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge_struct](struct.Layer.html#method.external_edge_struct) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge_struct(depth: u8, hash: u64, delta_depth: u8) -> ExternalEdge {
  get(depth).external_edge_struct(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge](struct.Layer.html#method.external_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  get(depth).external_edge(hash, delta_depth)
}

/// Conveniency function simply calling the [external_edge_sorted](struct.Layer.html#method.external_edge_sorted) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn external_edge_sorted(depth: u8, hash: u64, delta_depth: u8) -> Box<[u64]> {
  get(depth).external_edge_sorted(hash, delta_depth)
}

/// Conveniency function simply calling the [append_external_edge](struct.Layer.html#method.append_external_edge) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn append_external_edge(depth: u8, hash: u64, delta_depth: u8, dest: &mut Vec<u64>) {
  get(depth).append_external_edge(hash, delta_depth, dest)
}

/// Conveniency function simply calling the [append_external_edge_sorted](struct.Layer.html#method.append_external_edge_sorted) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn append_external_edge_sorted(depth: u8, hash: u64, delta_depth: u8, dest: &mut Vec<u64>) {
  get(depth).append_external_edge_sorted(hash, delta_depth, dest)
}

/// Conveniency function simply calling the [bilinear_interpolation](struct.Layer.html#method.bilinear_interpolation) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn bilinear_interpolation(depth: u8, lon: f64, lat: f64) -> [(u64, f64); 4] {
  get(depth).bilinear_interpolation(lon, lat)
}

/// Conveniency function simply calling the [cone_coverage_approx](struct.Layer.html#method.cone_coverage_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_approx(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get(depth).cone_coverage_approx(cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_coverage_flat](struct.Layer.html#method.cone_coverage_approx) method
/// of the [Layer] of the given *depth* and retrieving a flat array.
#[inline]
pub fn cone_coverage_approx_flat(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> Box<[u64]> {
  get(depth).cone_coverage_approx(cone_lon, cone_lat, cone_radius).to_flat_array()
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.cone_coverage_approx_custom) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_approx_custom(depth: u8, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get(depth).cone_coverage_approx_custom(delta_depth, cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.cone_coverage_centers) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_centers(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get(depth).cone_coverage_centers(cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.cone_coverage_fullin) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn cone_coverage_fullin(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
  get(depth).cone_coverage_fullin(cone_lon, cone_lat, cone_radius)
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.ring_coverage_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn ring_coverage_approx(depth: u8, cone_lon: f64, cone_lat: f64, cone_radius_int: f64, cone_radius_ext: f64) -> BMOC {
  get(depth).ring_coverage_approx(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
}

/// Conveniency function simply calling the [cone_coverage_approx_custom](struct.Layer.html#method.ring_coverage_approx) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn ring_coverage_approx_custom(depth: u8, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius_int: f64, cone_radius_ext: f64) -> BMOC {
  get(depth).ring_coverage_approx_custom(delta_depth, cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
}

/// Conveniency function simply calling the [elliptical_cone_coverage](struct.Layer.html#method.elliptical_cone_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn elliptical_cone_coverage(depth: u8, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
  get(depth).elliptical_cone_coverage( lon, lat, a, b, pa)
}

/// Conveniency function simply calling the [elliptical_cone_coverage_custom](struct.Layer.html#method.elliptical_cone_coverage_custom) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn elliptical_cone_coverage_custom(depth: u8, delta_depth: u8, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
  get(depth).elliptical_cone_coverage_custom(delta_depth, lon, lat, a, b, pa)
}

/// Conveniency function simply calling the [box_coverage](struct.Layer.html#method.box_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn box_coverage(depth: u8,  lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
  get(depth).box_coverage(lon, lat, a, b, pa)
}

/// Conveniency function simply calling the [polygon_coverage](struct.Layer.html#method.polygon_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn polygon_coverage(depth: u8, vertices: &[(f64, f64)], exact_solution: bool) -> BMOC {
  get(depth).polygon_coverage(vertices, exact_solution)
}

/// Conveniency function simply calling the [polygon_coverage](struct.Layer.html#method.polygon_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn custom_polygon_coverage(
  depth: u8,
  vertices: &[(f64, f64)],
  south_pole_method: &ContainsSouthPoleMethod,
  exact_solution: bool) -> BMOC
{
  get(depth).custom_polygon_coverage(vertices, south_pole_method, exact_solution)
}

/// Conveniency function simply calling the [zone_coverage](struct.Layer.html#method.zone_coverage) method
/// of the [Layer] of the given *depth*.
#[inline]
pub fn zone_coverage(depth: u8, lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> BMOC {
  get(depth).zone_coverage(lon_min, lat_min, lon_max, lat_max)
}

pub mod bmoc;
pub mod gpu;

use self::bmoc::*;
use super::ring::{triangular_number_x4};
use super::sph_geom::coo3d::*;
use super::sph_geom::{Polygon, ContainsSouthPoleMethod};
use super::sph_geom::cone::{Cone};
use super::sph_geom::zone::{Zone};
use crate::sph_geom::frame::RefToLocalRotMatrix;
use super::sph_geom::elliptical_cone::EllipticalCone;
use super::{direction_from_neighbour, edge_cell_direction_from_neighbour};
use super::compass_point::{Cardinal, CardinalSet, CardinalMap};
use super::compass_point::MainWind::{S, SE, E, SW, C, NE, W, NW, N};
use super::special_points_finder::{arc_special_points};


/// Defines an HEALPix layer in the NESTED scheme.
/// A layer is simply an utility structure containing all constants and methods related
/// to a given depth.
pub struct Layer {
  depth: u8,
  nside: u32,
  nside_minus_1: u32,
  n_hash: u64,
  twice_depth: u8,
  d0h_mask: u64,
  x_mask: u64,
  y_mask: u64,
  xy_mask: u64,
  time_half_nside: i64,
  one_over_nside: f64,
  // z_order_curve: &'static dyn ZOrderCurve,
  z_order_curve: ZOC,
}

impl Layer {
  
  const fn new(depth: u8, z_order_curve: ZOC) -> Layer {
    // const onef64_to_bits: u64 = 4607182418800017408_u64; // 1_f64.to_bits()
    // let time_nside = (depth as u64) << 52;
    let twice_depth: u8 = depth << 1u8;
    let nside: u32 = super::nside_unsafe(depth);
    Layer {
      depth,
      nside,
      nside_minus_1: nside - 1,
      n_hash: super::n_hash_unsafe(depth),
      twice_depth,
      d0h_mask: 15_u64 << twice_depth,
      x_mask: nested::x_mask(depth),
      y_mask: nested::y_mask(depth), //x_mask << 1,
      xy_mask: nested::xy_mask(depth),
      time_half_nside: (depth as i64 - 1) << 52,
      one_over_nside: ONE_OVER_NSIDE[depth as usize], // f64::from_bits(onef64_to_bits - time_nside),
      z_order_curve
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 12_u8;
  /// let nside = nside(depth) as u64;
  /// let nested12: &Layer = get(depth);
  /// assert_eq!(nside * nside - 1, nested12.hash(12.5_f64.to_radians(), 89.99999_f64.to_radians()));
  /// ```
  pub fn hash(&self, lon: f64, lat: f64) -> u64 {
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

  pub fn hash_checked(&self, lon: f64, lat: f64) -> Result<u64, String> {
    check_lat_res(lat)?;
    Ok(self.hash(lon, lat))
  }
  
  /* The clean way to do, but changes the API...
  pub fn hash(&self, lon: f64, lat: f64) -> Result<u64, Error> {
    check lon, lat
    hash_uncheked(lon, lat)
  }
  pub fn hash_uncheked(&self, lon: f64, lat: f64) -> u64 {
    ...
  }*/

  #[inline]
  fn d0h_lh_in_d0c(lon: f64, lat: f64) -> (u8, f64, f64) {
    let (x_pm1, q) = Layer::xpm1_and_q(lon);
    if lat > TRANSITION_LATITUDE {
      // North polar cap, Collignon projection.
      // - set the origin to (PI/4, 0)
      let sqrt_3_one_min_z = SQRT6 * (HALF * lat + FRAC_PI_4).cos();
      let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
      let d0h = q;
      (d0h, x_proj, y_proj)
    } else if lat < -TRANSITION_LATITUDE {
      // South polar cap, Collignon projection
      // - set the origin to (PI/4, -PI/2)
      let sqrt_3_one_min_z = SQRT6 * (HALF * lat - FRAC_PI_4).cos(); // cos(-x) = cos(x)
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
      let q01 = (x_pm1 >   y_pm1) as u8; /* 0/1 */  debug_assert!(q01 == 0 || q01 == 1);
      let q12 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */  debug_assert!(q12 == 0 || q12 == 1);
      let q1 = q01 & q12; /* = 1 if q1, 0 else */   debug_assert!( q1 == 0 ||  q1 == 1);
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
    let q = x as u8 | 1_u8;
    // Remark: to avoid the branch, we could have copied lon_sign on x - q, 
    //         but I so far lack of idea to deal with q efficiently.
    //         And we are not supposed to have negative longitudes in ICRS 
    //         (the most used reference system in astronomy).
    if lon_sign == 0 { // => lon >= 0
      (x - (q as f64), (q & 7_u8) >> 1)
    } else { // case lon < 0 should be rare => few risks of branch miss-prediction
      // Since q in [0, 3]: 3 - (q >> 1)) <=> 3 & !(q >> 1)
      // WARNING: BE SURE TO HANDLE THIS CORRECTLY IN THE REMAINING OF THE CODE!
      //  - Case lon =  3/4 pi = 270 deg => x = -1, q=3
      //  - Case lon = -1/2 pi = -90 deg => x =  1, q=2
      (q as f64 - x, 3 - ((q & 7_u8) >> 1))
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 12_u8;
  /// let nside = nside(depth) as u64;
  /// let nested12: &Layer = get(depth);
  /// let h_org = nside * nside - 1;
  /// let (h_ra, h_dec) = nested12.center(h_org);
  /// let (h, dx, dy) = nested12.hash_with_dxdy(h_ra, h_dec);
  /// assert_eq!(h_org, h);
  /// // A precision of 1e-12 in a cell of depth 12 (side of ~51.5 arcsec)
  /// // leads to an absolute precision of ~0.05 nanoarcsec
  /// assert!((dx - 0.5).abs() < 1e-12_f64);
  /// assert!((dy - 0.5).abs() < 1e-12_f64);
  /// ```
  pub fn hash_with_dxdy(&self, lon: f64, lat: f64) -> (u64, f64, f64) {
    check_lat(lat);
    let (d0h, l_in_d0c, h_in_d0c) = Layer::d0h_lh_in_d0c(lon, lat);
    let x = f64::from_bits((self.time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64);
    debug_assert!(-1e-14 < x || x < self.nside as f64 * (1.0_f64 + 1e-14), "x: {}, x_proj: {}; y_proj: {}", &x, &h_in_d0c, &l_in_d0c);
    let y = f64::from_bits((self.time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64);
    debug_assert!(-1e-14 < y || y < self.nside as f64 * (1.0_f64 + 1e-14), "y: {}, x_proj: {}; y_proj: {}", &y, &h_in_d0c, &l_in_d0c);
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
  fn div_by_nside_floor_u8(&self, val: u64) -> u8 {
    (val >> self.depth) as u8
  }
  
  #[inline]
  fn modulo_nside_u32(&self, val: u32) -> u32 {
    val & self.nside_minus_1
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
    debug_assert!(i < self.nside && j < self.nside, "nside: {}; i: {}, j: {}", self.nside, i, j);
    d0h_bits | self.z_order_curve.ij2h(i, j)
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
  /// use cdshealpix::nested::{get, Layer};
  /// let depth = 0;
  /// let n0 = get(depth);
  /// for h in 0..12 {
  ///   assert_eq!(n0.to_ring(h), h);
  /// }
  /// ```
  /// 
  /// At depth 1:
  /// ```rust
  /// use cdshealpix::nested::{get, Layer};
  /// 
  /// let depth = 1;
  /// let n1 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  /// 
  /// let depth = 2;
  /// let n2 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  /// let depth = 0;
  /// let n0 = get(depth);
  /// for h in 0..12 {
  ///   assert_eq!(n0.from_ring(h), h);
  /// }
  /// ```
  /// At depth 1:
  /// ```rust
  /// use cdshealpix::nested::{get, Layer};
  /// 
  /// let depth = 1;
  /// let n1 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  /// 
  /// let depth = 2;
  /// let n2 = get(depth);
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
        self.modulo_nside_u32(i_in_d0c as u32),
        self.modulo_nside_u32(j_in_d0c as u32) 
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get(depth);
  /// 
  /// assert!(dist(nested0.sph_coo(0, 0.5, 0.5) , nested0.center(0)) < 1e-15);
  /// ```
  ///
  pub fn sph_coo(&self, hash: u64, dx: f64, dy: f64) -> (f64, f64) {
    assert!((0.0..1.0).contains(&dx));
    assert!((0.0..1.0).contains(&dy));
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// fn dist(p1: (f64, f64), p2: (f64, f64)) -> f64 {
  ///   let sindlon = f64::sin(0.5 * (p2.0 - p1.0));
  ///   let sindlat = f64::sin(0.5 * (p2.1 - p1.1));
  ///   2f64 * f64::asin(f64::sqrt(sindlat * sindlat + p1.1.cos() * p2.1.cos() * sindlon * sindlon))
  /// }
  /// 
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 0u8;
  /// let nested0 = get(depth);
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
  
  /// Same method as [neighbours](#method.neighbours) except that neighbours are appended
  /// to the given vector.
  pub fn append_bulk_neighbours(&self, hash: u64, dest: &mut Vec<u64>) {
    self.check_hash(hash);
    let h_bits: HashBits = self.pull_bits_appart(hash);
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      self.append_bulk_edge_cell_neighbours(hash, dest);
    } else {
      self.append_bulk_inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, dest);
    }
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
  pub fn internal_edge_sorted(mut hash: u64, delta_depth: u8) -> Vec<u64> {
    if delta_depth == 0 {
      return vec![hash; 1];
    } else if delta_depth == 1 {
      hash <<= 2;
      return vec![hash, hash | 1, hash | 2, hash | 3];
    }
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
      if x == lim {
        k0 = k1;
        k1 += lim as usize; // +2 +4 +8
        lim <<= 1; // 4 8 16 32 ...
      }
      // To be tested if they are faster (since no risk of branch miss-prediction):
      // probably true for small deltaDepth but not for large deltaDepth.
     /* let mut tmp = x & lim;  debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == x));
      k0 += (tmp >> 1) as usize;
      k1 += tmp as usize;
      tmp -= x;           debug_assert!((x < lim            ) || (x == lim && tmp == 0));
      tmp = 1 >> tmp;     debug_assert!((x < lim && tmp == 0) || (x == lim && tmp == 1));
      lim <<= tmp;*/
    }
    result
  }
  
  /// Similar to [external_edge](#method.external_edge) except that the returned structure allow
  /// to access each elements of the external edge: the 4 corners plus the 4 edges.
  pub fn external_edge_struct(&self, hash: u64, delta_depth: u8) -> ExternalEdge {
    self.check_hash(hash);
    let mut res = ExternalEdge::new_empty();
    let h_bits: HashBits = self.pull_bits_appart(hash);
    let mut neighbours = MainWindMap::new();
    if self.is_in_base_cell_border(h_bits.i, h_bits.j) {
      // Not easy: opposite directions depends on base cell neighbours
      self.edge_cell_neighbours(hash, &mut neighbours);
      let h_parts: HashParts = self.decode_hash(hash);
      for (direction, hash_value) in neighbours.entries_vec().drain(..) {
        
       let dir_from_neig = if h_parts.d0h == self.h_2_d0h(hash_value) {
          direction.opposite()
        } else if self.depth == 0 {
	        direction_from_neighbour(h_parts.d0h, &direction)
        } else {
          let dir_in_basce_cell_border = self.direction_in_base_cell_border(h_bits.i, h_bits.j);
          edge_cell_direction_from_neighbour(h_parts.d0h, &dir_in_basce_cell_border, &direction)
        };
        add_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, &direction,&mut res);
      }
    } else {
      // Easy: always use the opposite direction
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
    let mut dest = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
    self.external_edge_generic(hash, delta_depth, false, &mut dest);
    dest.into_boxed_slice()
  }
  
  /// Similar to [external_edge](#method.external_edge) except that the returned list of cells is ordered.
  pub fn external_edge_sorted(&self, hash: u64, delta_depth: u8) -> Box<[u64]> {
    let mut dest = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
    self.external_edge_generic(hash, delta_depth, true, &mut dest);
    dest.into_boxed_slice()
  }

  /// Similar to [external_edge](#method.external_edge) except that the result is appended to the
  /// provided  `dest` list.
  pub fn append_external_edge(&self, hash: u64, delta_depth: u8, dest: &mut  Vec<u64>) {
    self.external_edge_generic(hash, delta_depth, false, dest);
  }

  /// Similar to [external_edge_sorted](#method.external_edge_sorted) except that the result is appended to the
  /// provided  `dest` list.
  pub fn append_external_edge_sorted(&self, hash: u64, delta_depth: u8,  dest: &mut Vec<u64>) {
    self.external_edge_generic(hash, delta_depth, true, dest);
  }

  fn external_edge_generic(&self, hash: u64, delta_depth: u8, sorted: bool, dest: &mut Vec<u64>) { //} -> Box<[u64]> {
    self.check_hash(hash);
    if delta_depth == 0 {
      if sorted {
        let mut tmp: Vec<u64> = Vec::with_capacity(8);
        self.append_bulk_neighbours(hash, &mut tmp);
        tmp.sort_unstable();
        dest.append(&mut tmp);
      } else {
        self.append_bulk_neighbours(hash, dest);
      }
      return;
    }
    
    // let mut edge = Vec::with_capacity((4 + (self.nside << 2)) as usize); // 4 borders (nside) + 4 corners (1)
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
        append_sorted_internal_edge_element(hash_value, delta_depth, dir_from_neig, dest);
      }
    } else {
      // Easy: always use the opposite direction
      let mut neighbours = MainWindMap::new();
      self.inner_cell_neighbours(h_bits.d0h, h_bits.i, h_bits.j, &mut neighbours);
      let mut neighbours = if sorted { neighbours.sorted_entries_vec() } else { neighbours.entries_vec() };
      for (direction, hash_value) in neighbours.drain(..) {
        append_sorted_internal_edge_element(hash_value, delta_depth,  direction.opposite(), dest);
      }
    }
    // edge.into_boxed_slice()
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

  fn append_bulk_inner_cell_neighbours(&self, d0h_bits: u64, i_in_d0h_bits: u64, j_in_d0h_bits: u64,
                           dest: &mut Vec<u64>) {
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
    dest.push(bits_2_hash(d0h_bits, im1_bits, jm1_bits));
    dest.push(bits_2_hash(d0h_bits, i_in_d0h_bits, jm1_bits));
    dest.push(bits_2_hash(d0h_bits, ip1_bits, jm1_bits));
    dest.push(bits_2_hash(d0h_bits, im1_bits, j_in_d0h_bits));
    dest.push(bits_2_hash(d0h_bits, ip1_bits, j_in_d0h_bits));
    dest.push(bits_2_hash(d0h_bits, im1_bits, jp1_bits));
    dest.push(bits_2_hash(d0h_bits, i_in_d0h_bits, jp1_bits));
    dest.push(bits_2_hash(d0h_bits, ip1_bits, jp1_bits));
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

  fn append_bulk_edge_cell_neighbours(&self, hash: u64, dest: &mut Vec<u64>) {
    let h_parts: HashParts = self.decode_hash(hash);
    let d0h = h_parts.d0h;
    let i = h_parts.i;
    let j = h_parts.j;
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, S) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, SE) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, E) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, SW) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, NE) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, W) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, NW) { dest.push(n); }
    if let Some(n) = self.neighbour_from_parts(d0h, i, j, N) { dest.push(n); }
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
    *j -= self.nside_minus_1 as i32;
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
  /// # Panics
  ///   If `lat` **not in** `[-pi/2, pi/2]`, this method panics.
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

  /// Returns a hierarchical view of the list of cells overlapped by the given zone.
  ///
  /// TODO: tests!!!!
  ///
  /// # Input
  /// - `lon_min` the longitude of the bottom left corner
  /// - `lat_min` the latitude of the bottom left corner
  /// - `lon_max` the longitude of the upper left corner
  /// - `lat_max` the latitude of the upper left corner
  ///
  /// # Output
  /// - the list of cells overlapped by the given zone, in a BMOC (hierarchical view also telling
  ///   if a cell is fully or partially covered).
  ///
  /// # Remark
  /// - If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
  /// - The north pole is included only if `lon_min == 0 && lat_max == pi/2`
  ///
  /// # Panics
  /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
  /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
  /// * `lat_min >= lat_max`.
  pub fn zone_coverage(&self, lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> BMOC {
    let zone = Zone::new(lon_min, lat_min, lon_max, lat_max);
    let (depth_start, hashs_start) = zone.smallest_enclosing_cone()
      .map(|mec| {
        let center = mec.center().lonlat();
        if zone.contains(center.lon(), center.lat()) && has_best_starting_depth(mec.radius()) {
          let d = best_starting_depth(mec.radius()).min(self.depth);
          let h = hash(d, center.lon(), center.lat());
          (d, neighbours(d, h, true).sorted_values_vec())
        } else {
          (0, (0..12).into_iter().collect::<Vec<u64>>())
        }
      }).unwrap_or((0, (0..12).into_iter().collect::<Vec<u64>>()));
    let zone_perimeter = 2.0 * (zone.dlon() + zone.dlat());
    let aprox_cell_side_at_depth_max = 2.0 * (PI / n_hash(self.depth) as f64).sqrt();
    let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, 6 * (2 + (zone_perimeter / aprox_cell_side_at_depth_max) as usize));
    let zone_vertices = zone.vertices();
    let zone_vertices_hashs: [u64; 4] = [
      self.hash(zone_vertices[0].0, zone_vertices[0].1),
      self.hash(zone_vertices[1].0, zone_vertices[1].1),
      self.hash(zone_vertices[2].0, zone_vertices[2].1),
      self.hash(zone_vertices[3].0, zone_vertices[3].1)
    ];
    for h in hashs_start {
      self.zone_coverage_recur(depth_start, h, &zone, &zone_vertices_hashs, &mut bmoc_builder);
    }
    bmoc_builder.to_bmoc()
  }

  fn zone_coverage_recur(
    &self,
    depth: u8, hash: u64,
    zone: &Zone, zone_vertices_hashs: &[u64; 4],
    bmoc_builder: &mut BMOCBuilderUnsafe)
  {
    let shift = (self.depth - depth) << 1; // <=> shift = 2 * (self.depth - depth)
    let [(l_s, b_s), (l_e, b_e), (l_n, b_n), (l_w, b_w)] = vertices(depth, hash);
    let n_in = zone.contains(l_s, b_s) as u8
      + zone.contains(l_e, b_e) as u8
      + zone.contains(l_n, b_n) as u8
      + zone.contains(l_w, b_w) as u8;
    // A cell may intersect a zone without having a vertex inside:
    // * either a vertex of the zone is in the cell: we can easily test this
    // * or no vertex of the zone is int the cell: we must detect either
    //     + the cell NS  vertical line intersect the zone => (N_lat > lat_max && S_lat < lat_min) && NS_lon_NS in dlon
    //     + the cell EW horizontal line intersect the zone => (EW_lat in dlat) && E_lon in ...
    // Remark: in the projected plane:
    // * Equatorial region: the zone is a rectangle
    // * Polar caps: the zone is trapeze
    if n_in == 4 {
      bmoc_builder.push(depth, hash, true);
    } else if n_in > 0
      || zone_vertices_hashs[0] >> shift == hash
      || zone_vertices_hashs[1] >> shift == hash
      || zone_vertices_hashs[2] >> shift == hash
      || zone_vertices_hashs[3] >> shift == hash
      || zone.crossed_vertically(l_n, b_n, b_s) // we may have false positive in the polar caps around n*PI/2
      || zone.crossed_vertically(l_s, b_n, b_s) // we may have false positive in the polar caps around n*PI/2
      || zone.crossed_horizontally(l_w, l_e, b_w) // we may have false positive in the polar caps around n*PI/2
      || zone.crossed_horizontally(l_w, l_e, b_e) // we may have false positive in the polar caps around n*PI/2
      // We should probably consider the zone in the projection plane, to be done later if necesary
    {
      if depth == self.depth {
        bmoc_builder.push(depth, hash, false);
      } else {
        let hash = hash << 2;
        let depth = depth + 1;
        self.zone_coverage_recur(depth, hash, zone, zone_vertices_hashs, bmoc_builder);
        self.zone_coverage_recur(depth, hash | 1_u64, zone, zone_vertices_hashs, bmoc_builder);
        self.zone_coverage_recur(depth, hash | 2_u64, zone, zone_vertices_hashs, bmoc_builder);
        self.zone_coverage_recur(depth, hash | 3_u64, zone, zone_vertices_hashs, bmoc_builder);
      }
    }
  }

  /// Returns a hierarchical view of the list of cells having their center in the given cone.
  /// All BMOC flags are set to 1.
  ///
  /// # Input
  /// - `cone_lon` the longitude of the center of the cone, in radians
  /// - `cone_lat` the latitude of the center of the cone, in radians
  /// - `cone_radius` the radius of the cone, in radians
  ///
  /// # Output
  /// - the list of cells having their center in the given cone, in a BMOC (hierarchical view with
  ///   all flags set to 1).
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 4_u8;
  /// let nested = get(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  ///
  /// let actual_res = nested.cone_coverage_centers(lon, lat, radius);
  /// let expected_res: [u64; 7] = [2058, 2059, 2080, 2081, 2082, 2083, 2088];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///      assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_coverage_centers(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      return self.allsky_bmoc_builder().to_bmoc();
    }
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    let shs_cone_radius = dbg!(to_squared_half_segment(cone_radius));
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius) {
      let distances = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for h in 0..12 {
        self.cone_coverage_centers_recur(
          0, h, shs_cone_radius,
          &shs_computer(cone_lon, cone_lat, cos_cone_lat),
           &minmax_array, 0, &mut bmoc_builder
        );
      }
      return bmoc_builder.to_bmoc_packing();
    }
    // Normal case
    let depth_start = best_starting_depth(cone_radius);
    if depth_start >= self.depth {
      let neigs: Vec<u64> = self.neighbours(self.hash(cone_lon, cone_lat), true)
        .values_vec().into_iter()
        .filter(|neigh| {
          let (lon, lat) = self.center(*neigh);
          squared_half_segment(lon - cone_lon, lat - cone_lat, lon.cos(), cos_cone_lat) <= shs_cone_radius
        }).collect();
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, true);
      }
      bmoc_builder.to_bmoc()
    } else {
        let distances = largest_center_to_vertex_distances_with_radius(
          depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius);
        let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
        let root_layer = get(depth_start);
        let root_center_hash = root_layer.hash(cone_lon, cone_lat);
        let neigs = root_layer.neighbours(root_center_hash, true);
        let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
        for &root_hash in neigs.sorted_values().iter() {
          self.cone_coverage_centers_recur(
            depth_start, root_hash, shs_cone_radius,
            &shs_computer(cone_lon, cone_lat, cos_cone_lat),
            &minmax_array, 0, &mut bmoc_builder
          );
        }
        bmoc_builder.to_bmoc_packing()
    }
  }

  // TODO: make a generic function with cone_coverage_approx_recur or at least cone_coverage_fullin_recur
  fn cone_coverage_centers_recur<F>(
    &self,
    depth: u8,             // cell depth
    hash: u64,             // cell hash
    shs_cone_radius: f64,  // squared_half_segment corresponding to the cone radius
    shs_computer: &F,      // function to compute the squared_half_segment between a point and the cone center
    shs_minmax: &[MinMax], // min and max shs at various delta depth
    recur_depth: u8,       // conveniency delta depth index
    bmoc_builder: &mut BMOCBuilderUnsafe  // builder in which cells are appended (flag always set to true!)
  )
    where F: Fn((f64, f64)) -> f64
  {
    let center = get(depth).center(hash);
    let shs = shs_computer(center);
    let MinMax{min, max} = shs_minmax[recur_depth as usize];
    if shs <= min {
      bmoc_builder.push(depth, hash, true);
    } else if shs <= max {
      if depth == self.depth {
        if shs <= shs_cone_radius {
          bmoc_builder.push(depth, hash, true);
        } // else do nothing, recur end here
      } else { // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
        // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
        let hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.cone_coverage_centers_recur(depth, hash, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_centers_recur(depth, hash | 1_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_centers_recur(depth, hash | 2_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_centers_recur(depth, hash | 3_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
      }
    }
  }

  /// Returns a hierarchical view of the list of cells fully covered by the given cone.
  /// All BMOC flags are set to 1.
  ///
  /// # Input
  /// - `cone_lon` the longitude of the center of the cone, in radians
  /// - `cone_lat` the latitude of the center of the cone, in radians
  /// - `cone_radius` the radius of the cone, in radians
  ///
  /// # Output
  /// - the list of cells fully covered by the given cone, in a BMOC (hierarchical view with
  ///   all flags set to 1).
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 4_u8;
  /// let nested3 = get(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  ///
  /// let actual_res = nested3.cone_coverage_fullin(lon, lat, radius);
  /// let expected_res: [u64; 2] = [2081, 2082];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///      assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_coverage_fullin(&self, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    // Special case: the full sky is covered
    if cone_radius >= PI {
      return self.allsky_bmoc_builder().to_bmoc();
    }
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    let shs_cone_radius = to_squared_half_segment(cone_radius);
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius) {
      let distances = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for h in 0..12 {
        self.cone_coverage_fullin_recur(
          0, h, shs_cone_radius,
          &shs_computer(cone_lon, cone_lat, cos_cone_lat),
          &minmax_array, 0, &mut bmoc_builder
        );
      }
      return bmoc_builder.to_bmoc_packing();
    }
    // Normal case
    let depth_start = best_starting_depth(cone_radius);
    if depth_start >= self.depth {
      let neigs: Vec<u64> = self.neighbours(self.hash(cone_lon, cone_lat), true)
        .values_vec().into_iter()
        .filter(|neigh| {
          let mut fully_in = true;
          for (lon, lat) in self.vertices(*neigh) {
            fully_in &= squared_half_segment(lon - cone_lon, lat - cone_lat, lon.cos(), cos_cone_lat) <= shs_cone_radius;
          }
          fully_in
        }).collect();
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, true);
      }
      bmoc_builder.to_bmoc()
    } else {
      let distances = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
      let root_layer = get(depth_start);
      let root_center_hash = root_layer.hash(cone_lon, cone_lat);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for &root_hash in neigs.sorted_values().iter() {
        self.cone_coverage_fullin_recur(
          depth_start, root_hash, shs_cone_radius,
          &shs_computer(cone_lon, cone_lat, cos_cone_lat),
          &minmax_array, 0, &mut bmoc_builder
        );
      }
      bmoc_builder.to_bmoc_packing()
    }
  }

  // TODO: make a generic function with cone_coverage_approx_recur or at least cone_coverage_centers_recur
  fn cone_coverage_fullin_recur<F>(
    &self,
    depth: u8,             // cell depth
    hash: u64,             // cell hash
    shs_cone_radius: f64,  // squared_half_segment corresponding to the cone radius
    shs_computer: &F,      // function to compute the squared_half_segment between a point and the cone center
    shs_minmax: &[MinMax], // min and max shs at various delta depth
    recur_depth: u8,       // conveniency delta depth index
    bmoc_builder: &mut BMOCBuilderUnsafe  // builder in which cells are appended (flag always set to true!)
  )
    where F: Fn((f64, f64)) -> f64
  {
    let center = get(depth).center(hash);
    let shs = shs_computer(center);
    let MinMax{min, max} = shs_minmax[recur_depth as usize];
    if shs <= min {
      bmoc_builder.push(depth, hash, true);
    } else if shs <= max {
      if depth == self.depth {
        // Center must be in the cone
        if shs <= shs_cone_radius {
          // And the 4 vertices must also be in the cone
          for vertex_coo in self.vertices(hash) {
            // If vertex out of the cone
            if shs_computer(vertex_coo) > shs_cone_radius {
              return; // do nothing, recur end here
            }
          }
          bmoc_builder.push(depth, hash, true);
        } // else do nothing, recur end here
      } else { // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
        // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
        let hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.cone_coverage_fullin_recur(depth, hash, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_fullin_recur(depth, hash | 1_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_fullin_recur(depth, hash | 2_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
        self.cone_coverage_fullin_recur(depth, hash | 3_u64, shs_cone_radius, shs_computer, shs_minmax, recur_depth, bmoc_builder);
      }
    }
  }

  /// Returns a hierarchical view of the list of cells fully covered by the given cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
  /// This is approximated: cell tagged as partially overlapped could be fully overlapped
  /// (but we guarantee that cells tagged as fully overlapped are really fully overlapped).
  ///
  /// # Input
  /// - `cone_lon` longitude of the center of the ring, in radians
  /// - `cone_lat` latitude of the center of the ring, in radians
  /// - `cone_radius_int` internal radius of the ring, in radians
  /// - `cone_radius_ext` external radius of the ring, in radians
  ///
  /// # Output
  /// - the list of cells overlapped by the given ring, in a BMOC (hierarchical view also telling
  ///   if a cell is fully or partially covered).
  ///
  /// # Panics
  /// * if `cone_radius_ext > PI`
  /// * if `cone_radius_ext <= cone_radius_int`
  ///
  /// # Usecase
  /// * annulus region provided by the delay between Fermi and INTEGRAL for GRBs
  ///   (asked for Gravitational Waves applications).
  ///
  ///
  /// # Example
  /// ```rust
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 4_u8;
  /// let nested = get(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius_int = 5.64323_f64.to_radians();
  /// let radius_ext = 10.0_f64.to_radians();
  ///
  /// let actual_res = nested.ring_coverage_approx(lon, lat, radius_int, radius_ext);
  /// let expected_res: [u64; 40] = [2050, 2051, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061,
  ///   2062, 2063, 2080, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2094, 2176,
  ///   2177, 2178, 2817, 2820, 2821, 2822, 2823, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2880];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///      assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn ring_coverage_approx(&self, cone_lon: f64, cone_lat: f64, cone_radius_int: f64, cone_radius_ext: f64) -> BMOC {
    self.ring_coverage_approx_internal(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
      .to_bmoc_packing()

  }

  pub fn ring_coverage_approx_custom(&self, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius_int: f64, cone_radius_ext: f64) -> BMOC {
    if delta_depth == 0 {
      self.ring_coverage_approx(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
    } else {
      // TODO: change the algo not to put all cell in the MOC before pruning it:
      // - make a second recur function returning the number of sub-cell overlapped by the cone
      get(self.depth + delta_depth)
        .ring_coverage_approx_internal(cone_lon, cone_lat, cone_radius_int, cone_radius_ext)
        .to_lower_depth_bmoc_packing(self.depth)
    }
  }

  fn ring_coverage_approx_internal(&self, cone_lon: f64, cone_lat: f64, cone_radius_int: f64, cone_radius_ext: f64) -> BMOCBuilderUnsafe {
    assert!(cone_radius_ext <= PI);
    assert!(cone_radius_int < cone_radius_ext, "{} >= {} ", cone_radius_int, cone_radius_ext);
    // Common variable
    let cos_cone_lat = cone_lat.cos();
    let shs_cone_radius_int = to_squared_half_segment(cone_radius_int);
    // Special case of very large radius: test the 12 base cells
    if !has_best_starting_depth(cone_radius_ext) {
      let distances_int = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius_int);
      let minmax_int_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius_int, &distances_int);
      let distances_ext = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius_ext);
      let minmax_ext_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius_ext, &distances_ext);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius_ext));
      for h in 0..12 {
        self.ring_coverage_approx_recur(
          0, h, &shs_computer(cone_lon, cone_lat, cos_cone_lat), shs_cone_radius_int,
          &minmax_int_array, &minmax_ext_array, 0, &mut bmoc_builder
        );
      }
      return bmoc_builder;
    }
    // Normal case
    let depth_start = best_starting_depth(cone_radius_ext);
    if depth_start >= self.depth {
      let shs_max = to_squared_half_segment(cone_radius_ext
        + largest_center_to_vertex_distance_with_radius(depth_start, cone_lon, cone_lat, cone_radius_ext));
      let root_layer = get(depth_start);
      let mut neigs: Vec<u64> = root_layer.neighbours(root_layer.hash(cone_lon, cone_lat), true)
        .values_vec().iter()
        .filter(|neigh| {
          for (lon, lat) in self.vertices(**neigh) {
            if squared_half_segment(lon - cone_lon, lat - cone_lat, lon.cos(), cos_cone_lat) > shs_cone_radius_int {
              return true
            }
          }
          false
        }) // remove cell fully in the internal raidus
        .map(h_to_h_and_shs(cone_lon, cone_lat, cos_cone_lat, root_layer))
        .filter(shs_lower_than(shs_max))
        .map(self.h_and_shs_to_lower_h(depth_start))
        .collect();
      neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
      neigs.dedup();         // remove duplicates (vector must be sorted first)
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, false);
      }
      bmoc_builder
    } else {
      let distances_int = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius_int);
      let minmax_int_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius_int, &distances_int);
      let distances_ext = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius_ext);
      let minmax_ext_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius_ext, &distances_ext);
      let root_layer = get(depth_start);
      let root_center_hash = root_layer.hash(cone_lon, cone_lat);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius_ext));
      for &root_hash in neigs.sorted_values().iter() {
        self.ring_coverage_approx_recur(
          depth_start, root_hash, &shs_computer(cone_lon, cone_lat, cos_cone_lat),shs_cone_radius_int,
          &minmax_int_array, &minmax_ext_array, 0, &mut bmoc_builder
        );
      }
      bmoc_builder
    }
  }

  fn ring_coverage_approx_recur<F>(
    &self,
    depth: u8,                 // cell depth
    hash: u64,                 // cell hash
    shs_computer: &F,          // function to compute the squared_half_segment between a point and the cone center
    shs_int_cone_radius: f64,  // squared_half_segment corresponding to the internal radius
    // shs_ext_cone_radius: f64,  // squared_half_segment corresponding to the external radius
    shs_int_minmax: &[MinMax], // min and max shs at various delta depth for the internal radius
    shs_ext_minmax: &[MinMax], // min and max shs at various delta depth for the external radius
    recur_depth: u8,           // conveniency delta depth index
    bmoc_builder: &mut BMOCBuilderUnsafe  // builder in which cells are appended (flag always set to true!)
  )
    where F: Fn((f64, f64)) -> f64
  {
    let center = get(depth).center(hash);
    let shs = shs_computer(center);
    let MinMax{min: min_int, max: max_int} = shs_int_minmax[recur_depth as usize];
    let MinMax{min: min_ext, max: max_ext} = shs_ext_minmax[recur_depth as usize];
    // We recall here that 'min_ext' and 'min_int' are lower bounds
    // while 'max_ext' and 'max_ext' are upper bounds
    if min_ext >= shs && shs >= max_int { // is fully in the ring (for sure)
      bmoc_builder.push(depth, hash, true);
    } else if max_ext >= shs && shs >= min_int { // may overlap (or not) or be fully in the ring
      if depth == self.depth {
        // If all 4 vertices are in the small radius, reject the cell,
        // i.e. if at least one vertex is not in the small radius, accept the cell
        for vertex_coo in self.vertices(hash) {
          if shs_computer(vertex_coo) > shs_int_cone_radius {
            bmoc_builder.push(depth, hash, false); // could be true, but we use a fast test only
            break;
          }
        }
      } else { //because 'min_int' and 'max_ext' are approximations, we may split cells fully
        // included in the MOC. We could check the 4 corers here, but we assume it is more costly
        // than post-merging (to be verified!)
        let hash = hash << 2;
        let depth = depth + 1;
        let recur_depth = recur_depth + 1;
        self.ring_coverage_approx_recur(depth, hash, shs_computer, shs_int_cone_radius, shs_int_minmax, shs_ext_minmax, recur_depth, bmoc_builder);
        self.ring_coverage_approx_recur(depth, hash | 1_u64, shs_computer, shs_int_cone_radius, shs_int_minmax, shs_ext_minmax, recur_depth, bmoc_builder);
        self.ring_coverage_approx_recur(depth, hash | 2_u64, shs_computer, shs_int_cone_radius, shs_int_minmax, shs_ext_minmax, recur_depth, bmoc_builder);
        self.ring_coverage_approx_recur(depth, hash | 3_u64, shs_computer, shs_int_cone_radius, shs_int_minmax, shs_ext_minmax, recur_depth, bmoc_builder);
      }
    }  // else fully out of the ring (for sure)
  }

  /// Returns a hierarchical view of the list of cells overlapped by the given cone.
  /// The BMOC also tells if the cell if fully or partially overlapped by the cone.
  /// This is approximated: cell tagged as partially overlapped could be fully overlapped
  /// (but we guarantee that cells tagged as fully overlapped are really fully overlapped).
  ///
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  ///
  /// let actual_res = nested3.cone_coverage_approx(lon, lat, radius);
  /// let expected_res: [u64; 10] = [512, 514, 515, 520, 521, 522, 544, 705, 708, 709];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get(depth);
  ///
  /// let lon = 13.158329_f64.to_radians();
  /// let lat = -72.80028_f64.to_radians();
  /// let radius = 5.64323_f64.to_radians();
  /// 
  /// let actual_res = nested3.cone_coverage_approx_custom(2, lon, lat, radius);
  /// let expected_res: [u64; 8] = [514, 515, 520, 521, 522, 705, 708, 709];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn cone_coverage_approx_custom(&self, delta_depth: u8, cone_lon: f64, cone_lat: f64, cone_radius: f64) -> BMOC {
    if delta_depth == 0 {
      self.cone_coverage_approx(cone_lon, cone_lat, cone_radius)
    } else {
      // TODO: change the algo not to put all cell in the MOC before pruning it
      get(self.depth + delta_depth)
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
      let distances = largest_center_to_vertex_distances_with_radius(
        0, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for h in 0..12 {
        self.cone_coverage_approx_recur(0, h,
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                                &minmax_array, 0, &mut bmoc_builder);
      }
      return bmoc_builder;
    }
    // Normal case
    let depth_start = best_starting_depth(cone_radius);
    if depth_start >= self.depth {
      let shs_max = to_squared_half_segment(cone_radius
        + largest_center_to_vertex_distance_with_radius(depth_start, cone_lon, cone_lat, cone_radius));
      let root_layer = get(depth_start);
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
      bmoc_builder
    } else {
      let distances = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, cone_lon, cone_lat, cone_radius);
      let minmax_array: Box<[MinMax]> = to_shs_min_max_array(cone_radius, &distances);
      let root_layer = get(depth_start);
      let root_center_hash = root_layer.hash(cone_lon, cone_lat);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(cone_radius));
      for &root_hash in neigs.sorted_values().iter() {
        self.cone_coverage_approx_recur(depth_start, root_hash,
                                &shs_computer(cone_lon, cone_lat, cos_cone_lat),
                                &minmax_array, 0, &mut bmoc_builder);
      }
      bmoc_builder
    }
  }

  fn cone_coverage_approx_recur<F>(
    &self,
    depth: u8,             // cell depth
    hash: u64,             // cell hash
    shs_computer: &F,      // function to compute the squared_half_segment between a point and the cone center
    shs_minmax: &[MinMax], // min and max shs at various delta depth
    recur_depth: u8,       // conveniency delta depth index
    bmoc_builder: &mut BMOCBuilderUnsafe // builder in which cells are appended
  )
    where F: Fn((f64, f64)) -> f64
  {
    let center = get(depth).center(hash);
    let shs = shs_computer(center);
    let MinMax{min, max} = shs_minmax[recur_depth as usize];
    if shs <= min { // Fast inclusion test
      bmoc_builder.push(depth, hash, true);
    } else if shs <= max { // Fast rejection test
      if depth == self.depth {
        bmoc_builder.push(depth, hash, false);
      } else { // because 'min' is a lower bound (approx), we may split cells fully included in the MOC
        // we could check the 4 corers here, we assume it is more costly than post-merging (to be verified!)
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
  #[allow(dead_code)]
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
    // 	NEW UPPER BOUND: supposedly more robust (and faster to compute)
    // But requires delta_depth to be given in parameter!
    // At lower resolution depth, max 9 cells (partially) overlapped:
    // => grid nside max = 3 * (2^DeltaDepth)
    // Worst case, for each nside max row (or col) at depth max:
    // - 2 (both sides) x (1 cell overllapping externaly + 1 cell overlapping internally + 1 no-fusioned internall cell)
    // - x2 to be conservative
    // 12 << delta_depth (delta_depth = diff between best starting depth and MOC depth max)

    // OLD UPPER BOUND (KEEP DURING THE TRANSITION): fails in rare circumstances,
    //   e.g. depth = 14, radius = 0.001, alpha = 0.002 ,delta = -1.3;
    // cell_area = 4 * pi / ncell = 4 * pi / (3 * 4 * nside^2) = pi / (3 * nside^2) =  pi * r^2
    // cell_radius = r = 1 / (sqrt(3) * nside)
    // As a very simple and naive rule, we take 4x the number of cells needed to cover
    // the cone external annulus
    // Annulus area = 4 pi ((R + r)^2 - R^2) = 4 pi (r^2 + 2rR)
    // N cells = 4 pi (r^2 + 2rR) / 4 pi r^2 = 1 + 2 R/r = 1 + 2 * sqrt(3) * nside * R
    const TWICE_SQRT_3: f64 = 2.0_f64 * 1.732_050_807_568_877_2_f64; // sqrt(3)
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get(depth);
  ///
  /// let lon = 36.80105218_f64.to_radians();
  /// let lat = 56.78028536_f64.to_radians();
  /// let a = 14.93_f64.to_radians();
  /// let b = 4.93_f64.to_radians();
  /// let pa = 75.0_f64.to_radians();
  /// 
  /// let actual_res = nested3.elliptical_cone_coverage(lon, lat, a, b, pa);
  /// let expected_res: [u64; 16] = [27, 30, 39, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 56, 57];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
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
      get(self.depth + delta_depth)
        .elliptical_cone_coverage_internal(lon, lat, a, b, pa)
        .to_lower_depth_bmoc_packing(self.depth)
    }
  }

  pub fn elliptical_cone_coverage_internal(&self, lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOCBuilderUnsafe {
    if a >= FRAC_PI_2 {
      panic!("Unable to handle ellipses with a semi-major axis > PI/2");
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
    let root_layer = get(depth_start);
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
        .map(|h| h >> ((depth_start - self.depth) << 1)) // h_to_lower_depth
        .collect();
      neigs.sort_unstable(); // sort the array (unstable is ok since we remove duplicates)
      neigs.dedup();         // remove duplicates (vector must be sorted first)
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, neigs.len());
      for neig in neigs {
        bmoc_builder.push(self.depth,neig, false);
      }
      bmoc_builder
    } else {
      let distances: Box<[f64]> = largest_center_to_vertex_distances_with_radius(
        depth_start, self.depth + 1, lon, lat, a);
      let neigs = root_layer.neighbours(root_center_hash, true);
      let mut bmoc_builder = BMOCBuilderUnsafe::new(self.depth, self.n_moc_cell_in_cone_upper_bound(a));
      for &root_hash in neigs.sorted_values().iter() {
        self.elliptical_cone_coverage_recur(depth_start, root_hash, &sph_ellipse, &distances, 0, &mut bmoc_builder);
      }
      bmoc_builder
    }
  }
  
  fn elliptical_cone_coverage_recur(&self, depth: u8, hash: u64, 
                                       ellipse: &EllipticalCone, distances: &[f64],
                                       recur_depth: u8, bmoc_builder: &mut BMOCBuilderUnsafe) {
    let (lon, lat) = get(depth).center(hash);
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

   /// Returns a hierarchical view of the list of cells overlapped by the given box.
   /// The BMOC also tells if the cell if fully or partially overlapped by the box.
   /// 
   /// # Input
   /// - `lon` the longitude of the center of the box, in radians
   /// - `lat` the latitude of the center of the box, in radians
   /// - `a` the semi-major axis of the box (half the box width), in radians
   /// - `b` the semi-minor axis of the box (half the box height), in radians
   /// - `pa` the position angle (i.e. the angle between the north and the semi-major axis, east-of-north), in radians
   /// 
   /// # Output
   /// - the list of cells overlapped by the given box, in a BMOC 
   ///   (hierarchical view also telling if a cell is fully or partially covered).
   ///
   /// # Panics
   /// - if `a` not in `]0, pi/2]`
   /// - if `b` not in `]0, a]`
   /// - if `pa` not in `[0, pi[`
   /// 
  pub fn box_coverage(&self,  lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> BMOC {
     let center = Coo3D::from_sph_coo(lon, lat);
     let vertices = Layer::box2polygon(center.lon(), center.lat(), a, b, pa);
     self.custom_polygon_coverage(&vertices, &ContainsSouthPoleMethod::ControlPointIn(center), true)
  }

  // # Panics
  // * if `lon` not in `[0, 2\pi[`
  // * if `lat` not in `[-\pi/2, \pi/2]`
  // * if `a` not in `]0, \pi/2]` or \`]0, \pi/2]`
  // * if `b` not in `]0, a]`
  // * if `pa` not in `[0, \pi[`
  fn box2polygon(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> Vec<(f64, f64)> {
    assert!((0.0..TWICE_PI).contains(&lon), "Expected: lon in [0, 2pi[. Actual: {}", lon);
    assert!((-HALF_PI..=HALF_PI).contains(&lat), "Expected: lat in [-pi/2, pi/2]. Actual: {}", lat);
    assert!(0.0 < a && a <= HALF_PI, "Expected: a in ]0, pi/2]. Actual: {}", a);
    assert!(0.0 < b && b <= a, "Expected: b in ]0, a]. Actual: {}", b);
    assert!((0.0..PI).contains(&pa), "Expected: pa in [0, pi[. Actual: {}", pa);
    // Compute spherical coordinates
    let frame_rotation = RefToLocalRotMatrix::from_center(lon, lat);
    // By application of the Thales theorem, the new point has the property:
    //   sin(new_lat) / sin(dlat) = (cos(dlon) * cos(new_lat)) / cos(dlat)
    // With (imagine looking a the triangle in the (x, z) plane:
    // * old_x = cos(lat)
    // * new_x = cos(dlon) * cos(new_lat)
    // * old_z = sin(dlat)
    // * new_z = sin(new_lat)
    // Leading to:
    //   tan(new_lat) = cos(dlon) * tan(dlat)
    let lon = a;
    let (sin_lon, cos_lon) = lon.sin_cos();
    let lat = (cos_lon * b.tan()).atan();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_pa, cos_pa) = pa.sin_cos();
    // Rotation by the position angle
    // - upper right (before rotation by PA)
    let (x1, y1, z1) = (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat);
    // - apply rotation (sin and cos are revere since theta = pi/2 - pa)
    let (y2, z2) = (
      y1 * sin_pa - z1 * cos_pa,
      y1 * cos_pa + z1 * sin_pa,
    );
    let mut vertices = Vec::with_capacity(4);
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower right (before rotation by PA) = (y1, -z1)
    let (y2, z2) = (
      y1 * sin_pa + z1 * cos_pa,
      y1 * cos_pa - z1 * sin_pa,
    );
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - lower left (before rotation by PA) = (-y1, -z1)
    let (y2, z2) = (
      -y1 * sin_pa + z1 * cos_pa,
      -y1 * cos_pa - z1 * sin_pa,
    );
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    // - upper left (before rotation by PA) = (-y1, z1)
    let (y2, z2) = (
      -y1 * sin_pa - z1 * cos_pa,
      -y1 * cos_pa + z1 * sin_pa,
    );
    vertices.push(frame_rotation.to_global_coo(x1, y2, z2));
    vertices
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
  ///
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
  /// use cdshealpix::nested::{get, Layer};
  ///
  /// let depth = 3_u8;
  /// let nested3 = get(depth);
  /// 
  /// let actual_res = nested3.polygon_coverage(&[(0.0, 0.0), (0.0, 0.5), (0.25, 0.25)], false);
  /// let expected_res: [u64; 8] = [304, 305, 306, 307, 308, 310, 313, 316];
  /// assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
  /// for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
  ///     assert_eq!(h1, *h2);
  /// }
  /// ```
  pub fn polygon_coverage(&self, vertices: &[(f64, f64)], exact_solution: bool) -> BMOC {
    self.custom_polygon_coverage(vertices, &ContainsSouthPoleMethod::Default, exact_solution)
  }

  /// Same as [polygon_coverage](struct.Layer.html#method.polygon_coverage) but with a custom
  /// Polygon builder (in case the default behaviour is not satisfactory).
  pub fn custom_polygon_coverage(
    &self,
    vertices: &[(f64, f64)],
    south_pole_method: &ContainsSouthPoleMethod,
    exact_solution: bool) -> BMOC
  {
    let poly = Polygon::new_custom(
      vertices.iter().map(|(lon, lat)| LonLat { lon: *lon, lat: *lat} )
        .collect::<Vec<LonLat>>().into_boxed_slice(),
      south_pole_method
    );
    let bounding_cone: Cone = Cone::bounding_cone(poly.vertices());
    let mut depth_start = 0;
    let neigs: Vec<u64> = if let ContainsSouthPoleMethod::Default = south_pole_method {
      if !has_best_starting_depth(bounding_cone.radius()) {
        // poly.must_contain(bounding_cone.center()); // Opposite of bounding cone out of?
        (0..12).collect()
      } else {
        depth_start = best_starting_depth(bounding_cone.radius()).min(self.depth);
        let root_layer = get(depth_start);
        let LonLat{lon, lat} = bounding_cone.center().lonlat();
        let center_hash = root_layer.hash(lon, lat);
        let mut neigs: Vec<u64> = root_layer.neighbours(center_hash, true).values_vec();
        neigs.sort_unstable();
        neigs
      }
    } else {
      (0..12).collect()
    };
    // Compute and sort the list of cells containing at least one polygon vertex
    let mut sorted_poly_vertices_hash = self.hashs_vec(poly.vertices());
    // Special treatment for the exact solution
    if exact_solution {
      let vertices= poly.vertices();
      let mut left= &vertices[vertices.len() - 1];
      for right in vertices {
        let special_lonlats = arc_special_points(left, right, 1.0e-14, 20);
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
      if (n_vertices_in_poly > 0 &&  n_vertices_in_poly < 4) || has_intersection(poly, poly_vertices) {
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
      } else if n_vertices_in_poly == 4 {
        moc_builder.push(depth, hash, true);
      }
    }
  }
  
  fn hashs_vec<T: LonLatT>(&self, poly_vertices: &[T]) -> Vec<u64> {
    poly_vertices.iter().map(|coo| self.hash(coo.lon(), coo.lat()))
      .collect::<Vec<u64>>()
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
  ((k << 2) + ( ((i as i8) + ((k - 1) >> 7)) & 3_i8)) as u8
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
       poly.is_intersecting_great_circle_arc(&vertices[2], &vertices[1]) // N vs E
    || poly.is_intersecting_great_circle_arc(&vertices[0], &vertices[1]) // S vs E
    || poly.is_intersecting_great_circle_arc(&vertices[3], &vertices[2]) // W vs N
    || poly.is_intersecting_great_circle_arc(&vertices[3], &vertices[0]) // W vs S
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
  0x0555555555555555_u64 >> (60 - (depth << 1))
}

/// mask ...101010
/// ```rust
/// use cdshealpix::nested::{y_mask, x_mask};
/// assert_eq!(y_mask(3), 0b00101010);
/// assert_eq!(y_mask(3), x_mask(3) << 1);
/// ```
#[inline]
pub const fn y_mask(depth: u8) -> u64{
  0x0AAAAAAAAAAAAAAA_u64 >> (60 - (depth << 1))
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
  // 0x0FFFFFFFFFFFFFFF_u64 >> (60 - (depth << 1))
  (1_u64 << (depth << 1)) - 1_u64
}

#[inline]
fn to_shs_min_max_array(cone_radius: f64, distances: &[f64]) -> Box<[MinMax]> {
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

#[cfg(test)]
mod tests {
  use super::*;
  
  #[test]
  fn testok_hash_d0() {
    let layer = get(0);
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
    let layer = get(3);
    let hash = layer.hash(333.5982493968911_f64.to_radians(), -25.919634217871433_f64.to_radians());
    assert_eq!(735_u64, hash);
  }

  #[test]
  fn testok_hash_2() {
    let layer = get(0);
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
    let layer = get(0);
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
    let layer = get(3);
    let ra = 180.0_f64;
    let dec = -45.85_f64;
    let hash = layer.hash(ra.to_radians(), dec.to_radians());
    assert_eq!(682_u64, hash);
  }
  
  #[test]
  fn testok_neighbour_d0() {
    let layer = get(0);
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
    let layer = get(0);
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
    
    let layer = get(depth);
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
  fn testok_external_edge_struct_v1() {
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
    let expected_res: [u64; 18] = [26, 27, 30, 31, 47, 53, 55, 61, 63, 69, 71, 77, 79, 90, 91, 94, 95, 143];
    // let expected_res: [u64; 20] = [63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 167, 173, 175, 239, 250, 251, 254, 255];
    //println!("{:?}", &actual_res);
    for (h1, h2) in actual_res.iter().zip(expected_res.iter()) {
      assert_eq!(h1, h2);
    }
    assert_eq!(expected_res.len(), actual_res.len());
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
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }
  
  #[test]
  fn testok_cone_approx_custom_bmoc_dbg() {
    let actual_res = cone_coverage_approx_custom(2, 1,20_f64.to_radians(), 0.0_f64.to_radians(), 50.0_f64.to_radians());
    let expected_res: [u64; 50] = [0, 1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 49, 52, 53, 64, 65, 66, 67,
      68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 82, 88, 89, 90, 91, 94, 131, 134, 135, 136,
      137, 138, 139, 140, 141, 142, 143, 181, 183, 189];
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_conecenter_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius = 5.64323_f64;
    let actual_res = cone_coverage_centers(depth, lon.to_radians(),lat.to_radians(), radius.to_radians());
    let expected_res: [u64; 7] = [2058, 2059, 2080, 2081, 2082, 2083, 2088];
    // println!("draw red circle({} {} {}deg)", lon, lat, radius);
    // to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_conefullin_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius = 5.64323_f64;
    let actual_res = cone_coverage_fullin(depth, lon.to_radians(),lat.to_radians(), radius.to_radians());
    let expected_res: [u64; 2] = [2081, 2082];
    println!("draw red circle({} {} {}deg)", lon, lat, radius);
    to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_ring_bmoc_dbg() {
    let depth = 4_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius_int = 5.64323_f64;
    let radius_ext = 10.0_f64;
    let actual_res = ring_coverage_approx(depth, lon.to_radians(),lat.to_radians(), radius_int.to_radians(), radius_ext.to_radians());
    let expected_res: [u64; 40] = [2050, 2051, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062,
      2063, 2080, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090, 2091, 2092, 2094, 2176, 2177,
      2178, 2817, 2820, 2821, 2822, 2823, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2880];
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_int);
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_ext);
    // to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
    for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_ring_bmoc() {
    let depth = 5_u8;
    let lon = 13.158329_f64;
    let lat = -72.80028_f64;
    let radius_int = 5.64323_f64;
    let radius_ext = 10.0_f64;
    let actual_res = ring_coverage_approx_custom(depth, 2, lon.to_radians(),lat.to_radians(), radius_int.to_radians(), radius_ext.to_radians());
    let expected_res: [u64; 99] = [8202, 8203, 8206, 8207, 8218, 8224, 8225, 8226, 8227, 8228, 8229,
      8230, 8231, 8232, 8233, 8234, 8236, 8237, 8239, 8240, 8241, 8242, 8243, 8246, 8248, 8249,
      8250, 8251, 8252, 8254, 8320, 8333, 8335, 8336, 8337, 8338, 8339, 8340, 8342, 8344, 8345,
      8346, 8347, 8348, 8355, 8356, 8357, 8358, 8359, 8360, 8361, 8362, 8363, 8364, 8365, 8366,
      8367, 8368, 8369, 8370, 8376, 8704, 8705, 8706, 8707, 8708, 11280, 11281, 11283, 11284,
      11285, 11286, 11287, 11292, 11293, 11328, 11329, 11330, 11331, 11332, 11333, 11334, 11335,
      11336, 11337, 11340, 11341, 11344, 11345, 11346, 11347, 11348, 11349, 11350, 11351, 11352,
      11353, 11520, 11521
    ];
    // For visual verification (cargo test testok_ring_bmoc -- --nocapture)
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_int);
    // println!("draw red circle({} {} {}deg)", lon, lat, radius_ext);
    // to_aladin_moc(&actual_res);
    assert_eq!(actual_res.flat_iter().deep_size(), expected_res.len());
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
    let lon = 0.0_f64;
    let lat = 0.0_f64;
    let a = 30.0_f64;
    let b = 5.0_f64;
    let pa = 30.0_f64;
    println!("draw red ellipse({} {} {}deg {}deg {}deg)", lon, lat, a, b, 90.0 - pa);
    let lon = lon.to_radians();
    let lat = lat.to_radians();
    let a = a.to_radians();
    let b = b.to_radians();
    let pa = pa.to_radians();
    let actual_res = elliptical_cone_coverage(4, lon, lat, a, b, pa);
    // let actual_res = elliptical_cone_coverage_custom(6, 2, lon, lat, a, b, pa);

    /*let expected_res: [u64; 112] = [130, 136, 138, 1035, 1038, 1039, 1050, 1051, 1056, 1057, 1058,
      1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 
      1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1120,
      1122, 1123, 1126, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1144, 1146, 1147, 1150,
      1153, 1156, 1157, 1159, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1177 ,1180, 1181,
      1183, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225 ,1226, 1227 ,1228 ,1229,
      1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244,
      1245, 1246, 1247, 1252, 1253, 1264, 1265, 1268, 2933, 2935, 2941];*/

    to_aladin_moc(&actual_res);

    let expected_res: [u64; 86] = [136, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064,
      1065, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1078, 1079, 1080, 1081,
      1082, 1083, 1084, 1085, 1086, 1087, 1122, 1123, 1128, 1129, 1130, 1131, 1132, 1133, 1134,
      1135, 1146, 1147, 1156, 1157, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1180, 1181,
      1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1227, 1228, 1229, 1230, 1231,
      1232, 1233, 1234, 1235, 1236, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247,
      2935];

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
    assert_eq!(actual_res.flat_iter().count(), 1716);
    /*for (h1, h2) in actual_res.flat_iter().zip(expected_res.iter()) {
      assert_eq!(h1, *h2);
    }
    assert_eq!(expected_res.len(), actual_res.flat_iter().count());*/
  }
  
  fn to_aladin_moc(bmoc: &BMOC) {
    print!("draw moc {}/", bmoc.get_depth_max());
    for cell in bmoc.flat_iter() {
      print!("{}, ", cell);
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
    // println!("draw moc 3/ {:?}", actual_res_approx.to_flat_array());
    assert_eq!(expected_res_approx.len(), actual_res_approx.deep_size());
    for (h1, h2) in actual_res_approx.flat_iter().zip(expected_res_approx.iter()) {
      assert_eq!(h1, *h2);
    }
    
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    // println!("draw moc 3/ {:?}", actual_res_exact.to_flat_array());
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
    let depth = 3; // 10
    
    let mut vertices = [(174.75937396073138, -49.16744206799886), 
      (185.24062603926856, -49.16744206799887), 
      (184.63292896369916, -42.32049830486584),
      (175.3670710363009, -42.32049830486584)];
    let expected_res_exact: [u64; 4] = [596, 597, 680, 682];

    to_radians(&mut vertices);

    let actual_res_exact = polygon_coverage(depth, &vertices, true);

    /*println!("@@@@@ FLAT VIEW");
    for cell in actual_res_exact.flat_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
    
    assert!(actual_res_exact.deep_size() > 0);
    
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_3() {
    // In Aladin: draw polygon(359.05552155,+04.00678949, 003.48359255,-04.16855897, 000.85415475,-04.30100561)
    let depth = 7; // 7
    
    let mut vertices = [(359.05552155, 4.00678949),
    (3.48359255, -4.16855897),
    (0.85415475, -4.30100561)];
    let expected_res_exact: [u64; 86] =[
      69479, 69484, 69485, 69486, 69487, 69488, 69489, 69490, 69491, 69492, 69494, 69496, 69497, 
      69498, 69499, 69500, 69501, 69502, 69503, 69572, 69573, 69574, 69575, 69581, 69584, 69585,
      69586, 69587, 69588, 69589, 69590, 69591, 69592, 69593, 69594, 69595, 69596, 69597, 69598,
      69599, 69617, 69619, 69620, 69621, 69622, 69623, 69628, 69629, 69631, 72322, 72328, 72330,
      72352, 72353, 72354, 72355, 72360, 72361, 72362, 72363, 72364, 72366, 75093, 77824, 77825,
      77826, 77827, 77828, 77830, 77831, 77833, 77835, 77836, 77837, 77838, 77839, 77860, 77861,
      77863, 77869, 77872, 77874, 77880, 77882, 77883, 77969];
    
    to_radians(&mut vertices);

    let actual_res_exact = polygon_coverage(depth, &vertices, true);

    to_aladin_moc(&actual_res_exact);

    println!("@@@@@ FLAT VIEW");
    for cell in actual_res_exact.flat_iter() {
    println!("@@@@@ cell a: {:?}", cell);
    }
    
    assert!(actual_res_exact.deep_size() > 0);
    
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_4() {
    // In Aladin: draw polygon(353.8156714, -56.33202193, 6.1843286, -56.33202193, 5.27558041, -49.49378172, 354.72441959, -49.49378172)
    let depth = 3;
    let mut vertices = [(353.8156714, -56.33202193),
                        (6.1843286, -56.33202193),
                        (5.27558041, -49.49378172),
                        (354.72441959, -49.49378172)];
    let expected_res_exact: [u64; 4] =[546, 552, 721, 724];
    to_radians(&mut vertices);
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_polygone_exact_5() {
    // In Aladin: draw polygon(,269.0489904481256, 50.93107840310317,256.0142359254377,
    // 13.990475156963315,243.27946433651516, -7.489608522084967,224.32882610807906,
    // -24.613170006428422,181.0522052001408, -32.62716788567398,196.13865549198817,
    // -48.194574021576855,204.11410185161034, -61.679032617555755,204.22079738688976,
    // -75.90259740652621,452.0599233625658, -81.01038799438487,291.81039030880413,
    // -59.94991219586425,298.63694405367926, -35.41386619268488,307.97031962597856,
    // -11.966387215315455,330.9614415037472, 19.1371582387119,310.5212786376428,
    // 20.26388923574062,296.2930605411064, 25.849021337885926,283.4739050648204,
    // 34.80711164805169,269.0489904481256, 50.93107840310317)
    let depth = 6;
    let vertices = [
      (4.695790732486566, 0.8889150097255261),
      (4.4682913488764395, 0.24417985540748033),
      (4.2460276551603116, -0.1307183283958091),
      (3.915276622719796, -0.42958085596528983),
      (3.159957098738856, -0.5694515052059677),
      (3.4232653287700523, -0.8411539982726408),
      (3.5624631270616547, -1.0765021986213243),
      (3.5643253154494583, -1.3247502355595913),
      (7.889934078990009, -1.4138979988201017),
      (5.093052102418385, -1.0463233540993349),
      (5.212197941830804, -0.6180885659230597),
      (5.375096075892637, -0.2088528564758103),
      (5.776366851387001, 0.33400642074068165),
      (5.4196187097295985, 0.35367158642311125),
      (5.171289457253199, 0.45115053076437905),
      (4.947552986866946, 0.6074987013677716)
    ];
    // let expected_res_exact: [u64; 4] =[546, 552, 721, 724];
    // to_radians(&mut vertices);
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    to_aladin_moc(&actual_res_exact);
    /*assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }*/
    /*to_degrees(&mut vertices);
    print!("\ndraw polygon(");
    for (lon, lat) in vertices.iter() {
      print!(",{}, {}", lon, lat);
    }
    println!(")");*/
  }

  #[test]
  fn testok_polygone_exact_6() {
    // In Aladin: draw polygon(268.84102386, -45.73283624, 299.40278164, 38.47909742, 66.0951825, 38.76894431, 96.66953469, -45.43470273)
    let depth = 0;
    let mut vertices = [
      (268.84102386, -45.73283624),
      (299.40278164, 38.47909742),
      (66.0951825, 38.76894431),
      (96.66953469, -45.43470273),
    ];
    let expected_res_exact: [u64; 9] =[0, 3, 4, 5, 7, 8, 9, 10, 11];
    to_radians(&mut vertices);
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
    /*to_degrees(&mut vertices);
    print!("\ndraw polygon(");
    for (lon, lat) in vertices.iter() {
      print!(",{}, {}", lon, lat);
    }
    println!(")");*/
  }

  #[test]
  fn testok_polygone_exact_fxp() {
    let depth = 10;
    let mut vertices = [ 
      (000.9160848095,  01.0736381331),
      (001.5961114529, -00.7062969568),
      (359.6079412529,  01.1296198985),
      (358.7836886856, -01.3663552354),
      (001.8720201899,  00.4097184220),
      (358.4159783831,  00.2376811155),
      (359.8319515193, -01.2824324848),
      (358.7798765255, 00.9896544935),
      (001.5440798843, 00.8056786162)
    ];
    // 10/4806091 10/4806094
    to_radians(&mut vertices);
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    // to_aladin_moc(&actual_res_exact);
    for h in actual_res_exact.flat_iter() {
      if h == 4806091 || h == 4806094 {
        assert!(false, "Contains 10/4806091 or 10/4806094")
      }
    }
  }
  
  
  /*
  Just used to test visually, everything was OK
  #[test]
  fn testok_polygon_aladinlitev3(){
    let depth = 2;
    let vertices = [
      (2.762765958009174, -1.1558461939389928), (1.8511822776529407, -1.0933557485463217),
      (1.2223938385258124, -0.9257482282529843), (0.7433555354994513, -0.7588156605217202),
      (0.3057973131669505, -0.6330807489082455), (-0.13178353851370578, -0.5615526603403407),
      (-0.5856691165836956, -0.5400327446156468), (-0.4015158856736461, -0.255385904573062),
      (-0.36447295312493533, 0.055767062790012194), (-0.46118312909272624, 0.3439014598749883),
      (-0.7077120187325645, 0.5509825177593931), (-1.0770630524112423, 0.6032198101231719),
      (-1.4301394011313524, 0.46923928877261123), (-1.4878998343615297, 0.8530569092752573),
      (-1.9282441740007021, 1.1389659504677638), (-2.7514745722416727, 1.0909408969626568),
      (-3.043766621461535, 0.7764645583461347), (-3.0294959096529364, 0.4089232125719977),
      (-2.878305764522211, 0.05125802264111997), (3.0734739115502143, 0.057021963147195695),
      (2.7821939201916965, -0.05879878325187485), (2.5647762800185414, -0.27363879379223127),
      (2.4394176374317205, -0.5545388864809319), (2.444623710270669, -0.8678251267471937)
    ];
    let mut s = String::new();
    for coo in vertices.iter().map(|(a, b)| [*a, *b]).flatten() {
      s.push(',');
      s.push_str(&format!("{}", (coo as f64).to_degrees()));
    }
    println!("draw polygon({})", &s.as_str()[1..]);
    // In Aladin: 
    // draw polygon(2.762765958009174, -1.1558461939389928, 1.8511822776529407, -1.0933557485463217, 1.2223938385258124, -0.9257482282529843, 0.7433555354994513, -0.7588156605217202,  0.3057973131669505, -0.6330807489082455, -0.13178353851370578, -0.5615526603403407,  -0.5856691165836956, -0.5400327446156468, -0.4015158856736461, -0.255385904573062,  -0.36447295312493533, 0.055767062790012194, -0.46118312909272624, 0.3439014598749883,  -0.7077120187325645, 0.5509825177593931, -1.0770630524112423, 0.6032198101231719,  -1.4301394011313524, 0.46923928877261123, -1.4878998343615297, 0.8530569092752573,  -1.9282441740007021, 1.1389659504677638, -2.7514745722416727, 1.0909408969626568,  -3.043766621461535, 0.7764645583461347, -3.0294959096529364, 0.4089232125719977,  -2.878305764522211, 0.05125802264111997, 3.0734739115502143, 0.057021963147195695, 2.7821939201916965, -0.05879878325187485, 2.5647762800185414, -0.27363879379223127, 2.4394176374317205, -0.5545388864809319, 2.444623710270669, -0.8678251267471937)
    // to_radians(&mut vertices);
    let actual_res_exact = polygon_coverage(depth, &vertices, true);
    to_aladin_moc(&actual_res_exact);
    let hash = hash(depth, 2.762765958009174, -1.1558461939389928);
    println!("hash: {}", hash);
  }*/
  
  
  #[test]
  fn testok_polygone_exact_markus() {
    // In Aladin: 
    // draw polygon(...)
    let depth = 12;
    let coos = [1.8513602782977292_f64, -0.18592998414051692, 1.8512258422891248,
      -0.18579783541887313, 1.8511638051632155, -0.185675857268537, 1.8511017577699869,
      -0.1856453564601204, 1.8510397110839012, -0.1856148549652769, 1.8509776667949591,
      -0.18557418856024166, 1.8509052799363221, -0.18554368392668075, 1.8508432374924337,
      -0.18550301603535885, 1.8507708544432935, -0.1854623454440125, 1.8507088117613895,
      -0.18543184028926155, 1.8506467721028483, -0.18539117022619508, 1.8505847358196907,
      -0.18534033525596078, 1.8505330419248691, -0.18528950207341022, 1.8504710184894149,
      -0.18519800896602515, 1.850419329785655, -0.18513701052017761, 1.8503883238565562,
      -0.18507601709824864, 1.8503573186288944, -0.18501502350745627, 1.8503263170225652,
      -0.1849438655291929, 1.8502849728262747, -0.1848828686900894, 1.8502643132238532,
      -0.184811713286104, 1.8502539912630134, -0.18475072501341952, 1.8502333294756885,
      -0.18468973371793343, 1.850233354215554, -0.18460841998467828, 1.8502333758556728,
      -0.18453727047097696, 1.850243746325784, -0.18443563131867075, 1.8502437739504332,
      -0.18434415338072863, 1.850243801563952, -0.18425267544669602, 1.850233483954607,
      -0.18418152293981005, 1.8502335055540653, -0.18411037344041434, 1.850223185237512,
      -0.18404938513277874, 1.8501818522075997, -0.18397822335057393, 1.8501405202635473,
      -0.18390706126594134, 1.8500991860698157, -0.18384606309162968, 1.8500475110573917,
      -0.1837952254787585, 1.8499958370161058, -0.1837443873907134, 1.8499441675707649,
      -0.18368338461546654, 1.84990284099673, -0.18361222078765688, 1.8498718533556298,
      -0.18354106034805112, 1.8498305249747435, -0.18348006020220078, 1.8498098766287663,
      -0.1834089031712837, 1.8497995779550107, -0.18330725725741462, 1.8497892640149673,
      -0.1832462681693668, 1.8497892913795262, -0.18317511869547, 1.8498099966304176,
      -0.18309381264003968, 1.8498410342777623, -0.18302267446611675, 1.8498927334119273,
      -0.18298203607383598, 1.8499651015965932, -0.1829515684700408, 1.850037472401915,
      -0.18291093572064188, 1.8500995055507055, -0.18287029878505212, 1.8501408712721863,
      -0.18280932641072814, 1.8501718997230883, -0.1827483506134699, 1.850182255036083,
      -0.18268736846694722, 1.8501305967663695, -0.18261620321629635, 1.8500582543810962,
      -0.18258568769748984, 1.8499962452637886, -0.18256533891575163, 1.8499238933743258,
      -0.18256531429383693, 1.8498721909927536, -0.18262628139511278, 1.8498308160127352,
      -0.1827075802042478, 1.8497894475832763, -0.18276855028745648, 1.8497377337273002,
      -0.18284984456644895, 1.8496860346994894, -0.1828904815280938, 1.8496239939937291,
      -0.1829412780938848, 1.8495412828905007, -0.18298190121966892, 1.8494792487450675,
      -0.1830123677623018, 1.8494172231098334, -0.1830225051988025, 1.8493448651057844,
      -0.1830224730958927, 1.8492828585330476, -0.18299195220809547, 1.8492828925303078,
      -0.18292080274479733, 1.8492829216622935, -0.18285981749099187, 1.8492829556396726,
      -0.1827886680286644, 1.8492829944578515, -0.18270735435791027, 1.8492830284122925,
      -0.18263620489631682, 1.8492933935110967, -0.1825752244003239, 1.849324430143685,
      -0.1825142533027064, 1.8493347945288308, -0.182453272729678, 1.849365839361594,
      -0.18237197298389143, 1.849396873991587, -0.1823110014833154, 1.8494382433000762,
      -0.18225003429912384, 1.8495106086441788, -0.18220940835074173, 1.849572633374496,
      -0.18217894145457372, 1.849655331928133, -0.18213831785758774, 1.8497276901408222,
      -0.18210785331326518, 1.849800051409277, -0.1820672236239383, 1.8498724078572921,
      -0.18203675720966933, 1.8499344285723411, -0.1820062863060658, 1.8500171252356468,
      -0.1819554931554334, 1.850079147551529, -0.18191485643812733, 1.8501825076068366,
      -0.18187423161174024, 1.8502445277525188, -0.18183359306249244, 1.8503168814810491,
      -0.18179295669421697, 1.8503788941242136, -0.18177264507944257, 1.8504512486801108,
      -0.18172184276349443, 1.8505029333185445, -0.1816710345197491, 1.8505442852757363,
      -0.1816100591568629, 1.850585633907972, -0.1815592476998868, 1.8506476565420458,
      -0.1814779478344996, 1.8506786715997467, -0.18141696921824987, 1.850709688119604,
      -0.18134582621623202, 1.8507303681750746, -0.18128484516779256, 1.8507510498462858,
      -0.18121369982844085, 1.8507407311398347, -0.18114254826722087, 1.8507304105793165,
      -0.18108156089920815, 1.8506994263298109, -0.1810104050299821, 1.8506477725552963,
      -0.18095957295895188, 1.8505857869062527, -0.18090873805888288, 1.850554803069368,
      -0.18084774560122668, 1.8505858179298913, -0.1807766032700415, 1.8506168296497307,
      -0.18071562497989227, 1.850647840673332, -0.18065464651625232, 1.8507201778206104,
      -0.180603840703177, 1.8507821733848713, -0.18059368882855678, 1.8508338422657715,
      -0.18055304172508826, 1.8508648499351492, -0.1804920620550578, 1.850864864843558,
      -0.1804107483219105, 1.8508648797465577, -0.1803294345853517, 1.850864892782244,
      -0.18025828506295855, 1.8508752354417124, -0.18019730159328218, 1.8508959092420871,
      -0.1801363199016534, 1.8509269155834893, -0.18006517565476338, 1.8509992407064624,
      -0.18002453044035827, 1.8510612321725082, -0.1799940470354819, 1.8511335526348274,
      -0.17997372852525145, 1.851195539039006, -0.1799737362930964, 1.851278187577929,
      -0.1799737455840484, 1.851422822521087, -0.17997375891123377, 1.8514951406828215,
      -0.17996359995564287, 1.8515674580202341, -0.17996360428695535, 1.8516707685022664,
      -0.179963608856252, 1.851743085839695, -0.17996361092195487, 1.8518050721289219,
      -0.17996361195004668, 1.8518670584181505, -0.17995344807278466, 1.8519393756213955,
      -0.17995344760633805, 1.8520116928246395, -0.17995344620699838, 1.8520840104304412,
      -0.1799636080947026, 1.85217699158969, -0.17999409638466657, 1.8522596433499237,
      -0.18003474808418177, 1.8523216338197017, -0.18007540027750576, 1.8523836232924833,
      -0.1800957233463537, 1.8524559467486013, -0.18013637309113859, 1.8525282700518246,
      -0.18016685768320762, 1.8525902599855844, -0.1801770142472329, 1.852662584575065,
      -0.1802074971054239, 1.8527349067467758, -0.18021765059379727, 1.8528175607341215,
      -0.18022780140668926, 1.8529105465946833, -0.18023794908706292, 1.8529311978668475,
      -0.18017695989513224, 1.8529001920025943, -0.1801159803832753, 1.8528485225687163,
      -0.1800448401341806, 1.8527865285329002, -0.1800041937674299, 1.8527245370097851,
      -0.17997371093336756, 1.8526625447015321, -0.17993306319323382, 1.8525902237418577,
      -0.17990258032858317, 1.852528235689744, -0.17988225954154966, 1.8525902170349364,
      -0.1798517592273769, 1.8526522020591731, -0.17985175088889574, 1.852745179595501,
      -0.17985173709613428, 1.8528071646197006, -0.17985172704426694, 1.8528691533611787,
      -0.17987204474727547, 1.8529414794564334, -0.17992285245321354, 1.8530034737194898,
      -0.179963497108399, 1.8530757955426747, -0.17998381042047257, 1.853137789136592,
      -0.1800142893692613, 1.8532307820893295, -0.18006508861459566, 1.8532927774692713,
      -0.1800955658471284, 1.8533444476726812, -0.180146373561122, 1.8533857958125843,
      -0.1802276762574201, 1.8534271393909703, -0.18028865021083082, 1.8534684809191593,
      -0.18033945964033768, 1.8534994974343275, -0.18041060021371805, 1.8535305147629872,
      -0.1804817406124282, 1.8535718618956585, -0.18054271348844708, 1.8535925454408997,
      -0.1806036924563041, 1.8535615675880976, -0.180664687201364, 1.8534995824661735,
      -0.180695198243738, 1.8534375908057812, -0.18070538017125248, 1.853365263656628,
      -0.1807053999696568, 1.8532929365074355, -0.18070541883480354, 1.8532309443469994,
      -0.1807155984771237, 1.8531689592273146, -0.180756270078698, 1.853168976197942,
      -0.18082741958169457, 1.8531999916204283, -0.18089856179674516, 1.8532516711771596,
      -0.18095953455695846, 1.8533240071244552, -0.18097984465025554, 1.8533963464617944,
      -0.18101031802316253, 1.8534583475915447, -0.18102046497809307, 1.8535203520344468,
      -0.1810407754602831, 1.8535926932809386, -0.18107124629742702, 1.853551379553357,
      -0.1811322441480593, 1.853468719459801, -0.18115259680598478, 1.8533860559168442,
      -0.1811627840320468, 1.853313725401153, -0.18117296738149677, 1.8532413971792996,
      -0.1811933140107519, 1.8531897401786994, -0.18123398329392268, 1.853158754286251,
      -0.18129497580183623, 1.8531587711626412, -0.18136612529137458, 1.8531794557183354,
      -0.18143726997806436, 1.853210472282553, -0.1814982479071053, 1.8532724917514882,
      -0.18155921797476143, 1.8533345126147671, -0.18162018735461194, 1.8533758634520165,
      -0.18167099743739729, 1.853458540674906, -0.18168113878211153, 1.8535308838873166,
      -0.1816912819845117, 1.8536135584443616, -0.18169125683075096, 1.8536755677428784,
      -0.18170140137622234, 1.8537375632869104, -0.18167088839128523, 1.8537995509138552,
      -0.1816200462987067, 1.8538615448392612, -0.18158953194282895, 1.8539028623273504,
      -0.1815386960249918, 1.8539545090822382, -0.18147769175202094, 1.853995824861937,
      -0.18142685514865847, 1.8540578073165581, -0.18137601020034316, 1.85410945917895,
      -0.18133533292176704, 1.85418177740647, -0.18129464667012765, 1.8542541168289186,
      -0.18130478054090335, 1.8542851446671542, -0.18136575251690873, 1.854274842432865,
      -0.18143690644015867, 1.8542438675549624, -0.18149790493589874, 1.854212891989415,
      -0.1815589032600803, 1.8541819157359403, -0.18161990141274395, 1.854140604494005,
      -0.18168090359572794, 1.854088957915552, -0.18174190958079806, 1.854026967573646,
      -0.18178259066192382, 1.8539649723831393, -0.1818131068473544, 1.853902976498352,
      -0.1818436223473159, 1.8538409873050585, -0.1818944655828114, 1.853789335355931,
      -0.18195546880818889, 1.8537583521286545, -0.1820164646202317, 1.8537583733736718,
      -0.1820774498829631, 1.8537583981667958, -0.18214859935568553, 1.8537274204772618,
      -0.18222992341714173, 1.8536860961411439, -0.18228075804079932, 1.8536344321829328,
      -0.18232143141876903, 1.8535724222652281, -0.18233161484398933, 1.8535000734963125,
      -0.18233163639387198, 1.8534483839488172, -0.1822909943717122, 1.8534173601355401,
      -0.18223001777030146, 1.8534173427097982, -0.1821690325050502, 1.8534069901763455,
      -0.18210805008898873, 1.8533759678871828, -0.182047073256815, 1.8532209657591807,
      -0.1821385907538852, 1.8532209657591807, -0.1821385907538852, 1.8531589595755087,
      -0.18215893381474696, 1.853117630451771, -0.18220976424895258, 1.8530659629067465,
      -0.18225043238787383, 1.8530246322191117, -0.18230126213656264, 1.8529833049442839,
      -0.18237242000328946, 1.8529316368446482, -0.1824232511157904, 1.8528799677781902,
      -0.18247408175220095, 1.8528282977445811, -0.182524911912516, 1.8527456170839731,
      -0.18256558233245482, 1.852693944920749, -0.1826164112546596, 1.8526215951840501,
      -0.18263674988737813, 1.8525492436239446, -0.18264692337438293, 1.8524768883584832,
      -0.18262660329194685, 1.8523838605427814, -0.18259612000798175, 1.8522908361603814,
      -0.18259612781505855, 1.8522184838629452, -0.18259613281959508, 1.852104788287783,
      -0.18261646722026986, 1.8520531092781092, -0.18265712602218745, 1.8520427750984183,
      -0.18271811162760823, 1.8521151290255806, -0.1827181088916202, 1.8521771466774282,
      -0.18271810580290904, 1.8522495006045756, -0.1827181013319043, 1.8523218562395611,
      -0.18273842435046514, 1.8523632070863645, -0.18279940611373796, 1.8524355642360182,
      -0.18281972766392637, 1.8524975830524693, -0.18281972102856434, 1.8525906112771218,
      -0.18281970978849002, 1.8526526286180132, -0.18280953722531323, 1.8527146457259827,
      -0.18279936397582477, 1.8527869955569074, -0.182768859929351, 1.8528490119468843,
      -0.1827586851930074, 1.8529110241835405, -0.18272818134741206, 1.8529730357223029,
      -0.1826976768158609, 1.8530350465633145, -0.18266717159837986, 1.8530970567067129,
      -0.18263666569499454, 1.8531590661526396, -0.18260615910573078, 1.8531900594433808,
      -0.1825451666039551, 1.8532003803817598, -0.18248417888771545, 1.8531693557982265,
      -0.18241303670180406, 1.853221024111172, -0.18237236760955008, 1.8532933734227546,
      -0.18237234966205038, 1.8533553899086845, -0.18238249774591547, 1.8534174124410632,
      -0.1824129735649136, 1.8534484373148825, -0.1824739501651511, 1.853448458080047,
      -0.18254509964029003, 1.8534484758839718, -0.18260608490476157, 1.8535001684529953,
      -0.18264672692445688, 1.8535621948350693, -0.1826772011390782, 1.8536241988278674,
      -0.18263652519325926, 1.853655187053907, -0.18257553012162162, 1.8536655028126596,
      -0.18251454155155494, 1.8536758183442161, -0.18245355296267193, 1.8537171406126367,
      -0.18239255420394243, 1.8537998077329572, -0.1823417052451024, 1.8538514706437519,
      -0.18230103034319703, 1.8538721153030797, -0.18222987351480946, 1.853872092765795,
      -0.18216888825415814, 1.8538720664800044, -0.18209773878303329, 1.8539237226908758,
      -0.1820468990060636, 1.8539857287020418, -0.18203671169602123, 1.8540580694387256,
      -0.1820265196690583, 1.854140749099136, -0.18202648673473754, 1.8542130938019286,
      -0.1820264569165996, 1.854337108668429, -0.18201623941775788, 1.8544094484765554,
      -0.1820060428557937, 1.8544507543028441, -0.181934874485717, 1.8544610551126304,
      -0.18186372024756248, 1.85449202468799, -0.18179255634714783, 1.8545229984280474,
      -0.18173155648472142, 1.8545746400423422, -0.18167054650397377, 1.854626285662072,
      -0.18161970025692645, 1.8546882749522515, -0.18158917676023903, 1.8547295847169856,
      -0.1815383347544417, 1.8547398865501494, -0.1814773442116756, 1.854750182777349,
      -0.18140618943984393, 1.8547708123902436, -0.18133502930381967, 1.8547811134037422,
      -0.18127403868276232, 1.854812081019679, -0.18121303722084287, 1.8548843980213114,
      -0.18118250611807843, 1.8549567315037854, -0.18118246671147994, 1.8549980941266584,
      -0.18123326482089308, 1.8550394575272688, -0.18128406262345753, 1.8551118226373557,
      -0.1813348422576304, 1.855142872491283, -0.18141613789159422, 1.8551635831794189,
      -0.18148727523052058, 1.8551739541751902, -0.18154825439023753, 1.855173997427532,
      -0.18161940384446773, 1.8551740345118601, -0.1816803890900812, 1.8551844122249557,
      -0.18175153243476497, 1.8551844680618457, -0.18184301030035063, 1.855194846250852,
      -0.18191415362325747, 1.8552052247243542, -0.18198529692593465, 1.85521559721333,
      -0.18204627600163992, 1.8552259699446727, -0.18210725505748945, 1.8552570133751127,
      -0.18216822162789617, 1.8552983928572575, -0.18222918170647745, 1.8553604311072789,
      -0.18226980021820402, 1.8554328057910092, -0.1823104114751146, 1.8555051679248793,
      -0.18233069338463606, 1.8555878661791823, -0.18235096750615315, 1.8556498867316906,
      -0.18236109019332603, 1.8557325932425572, -0.1823915263827604, 1.855794629261471,
      -0.18242197587644818, 1.8558773224810987, -0.18243208151878457, 1.855949672603792,
      -0.18243202927598628, 1.855980726089584, -0.18249299182360623, 1.8559394208942452,
      -0.18254384304099072, 1.8558877788897787, -0.1825947013807476, 1.8558361359239577,
      -0.18264555924593065, 1.8557741484255865, -0.1826862596482374, 1.8557018381076122,
      -0.1827472948592403, 1.8556398484486467, -0.1827879937788979, 1.8555675213744636,
      -0.18282869884964942, 1.8554951932178876, -0.1828694029882902, 1.8554332140249616,
      -0.18293042803810927, 1.8553919136685315, -0.18300160364397636, 1.855381623102711,
      -0.18307275958212196, 1.8553816759809545, -0.18315407323467503, 1.8553920595113906,
      -0.18322521618783746, 1.855392099311757, -0.18328620142857147, 1.8553818082597895,
      -0.18335735737122666, 1.8553508416788267, -0.18342852618988437, 1.8553302055429384,
      -0.18348952425176293, 1.8553095689502415, -0.18355052223910293, 1.8552682623384746,
      -0.18362169694954694, 1.8552372864379694, -0.1836827009405365, 1.8552166480855876,
      -0.18374369859095774, 1.8552166986250527, -0.18382501226104841, 1.8552167428619664,
      -0.18389616172419765, 1.855206442141301, -0.18395715315164085, 1.8552064862673985,
      -0.18402830261842207, 1.8552168693296875, -0.1840994459148106, 1.8552375853711334,
      -0.18416041877281975, 1.855237629947007, -0.18423156824474607, 1.855237680907982,
      -0.18431288192966777, 1.8552377318872193, -0.18439419561762005, 1.8552274305273724,
      -0.18445518709732775, 1.8552171352810847, -0.18452634277208088, 1.8552068397748254,
      -0.18459749843080292, 1.855206884025601, -0.1846686479184675, 1.8551965818312721,
      -0.1847296393511224, 1.8551759642309624, -0.18483129373905474, 1.8551760081158577,
      -0.1849024432368502, 1.8551657051572494, -0.18496343462191608, 1.8551347639410862,
      -0.18509558760038614, 1.8551244601711432, -0.18515657891799897, 1.8550728040399498,
      -0.18523792247843815, 1.8550418228732382, -0.1853090896669299, 1.8549901581162065,
      -0.1853802682587192, 1.854948821818579, -0.18543111226157105, 1.854876454267602,
      -0.18547180849782238, 1.8548247857687108, -0.18554298557650512, 1.8547420794541367,
      -0.18559384975419665, 1.8546593662655284, -0.18563454849317843, 1.8545662944694563,
      -0.1856447583426744, 1.8545042480799074, -0.18565495212383804, 1.8544318595148424,
      -0.18566514996322853, 1.854359475413111, -0.1856855110893831, 1.8542767670714466,
      -0.18574653263011986, 1.8542250834943699, -0.18580753997754113, 1.8541837366695544,
      -0.18585837835974803, 1.8541527401739855, -0.18592954067736409, 1.8541217385690594,
      -0.18599053860317263, 1.8540803893603461, -0.18604137623222627, 1.8540287089531733,
      -0.1861227102396788, 1.8539873578459334, -0.18617354718929408, 1.8539460218020507,
      -0.18626504074818293, 1.853915015919087, -0.1863260375513494, 1.8538736698208658,
      -0.18639720213192207, 1.8538219791307038, -0.18646837003308542, 1.8537806343526355,
      -0.18654969816950612, 1.8537392847482848, -0.1866208617760621, 1.8536772461860471,
      -0.18669203180853772, 1.8536048585853184, -0.18675304014573832, 1.8535324629682506,
      -0.18679371908332434, 1.8534600662469602, -0.18683439708737642, 1.8533876568005647,
      -0.18683441721708252, 1.8533566071093006, -0.1867734401494385, 1.8533565871984605,
      -0.18670229051176154, 1.8533565672938108, -0.18663114087943428, 1.8533668824835334,
      -0.18652949579605713, 1.8533772059318696, -0.18645834339686834, 1.8533771828928605,
      -0.18637702955217636, 1.853387505920691, -0.18630587714490543, 1.8533978286788215,
      -0.18623472472359145, 1.85339781117922, -0.18617373935273773, 1.8533977849385794,
      -0.18608226130288735, 1.8533977645364603, -0.18601111171378606, 1.8533977354016076,
      -0.18590946945138145, 1.8533977150149308, -0.18583831987299143, 1.8533976946346133,
      -0.18576717029885748, 1.8533976742606526, -0.18569602072888974, 1.8533562921093778,
      -0.18564521076173426, 1.853294235280328, -0.185624898469215, 1.8532321815180122,
      -0.1856147497127729, 1.8531701304682395, -0.18561476449242317, 1.8531184331532222,
      -0.1856655974032131, 1.8530563905581923, -0.18570626782034258, 1.85298400067569,
      -0.18572661117418918, 1.8529116062702775, -0.18572662514312013, 1.852828869806907,
      -0.18572663996183605, 1.8527254475954322, -0.18571649254160064, 1.8526427083397266,
      -0.18569617616057874, 1.8525806522761148, -0.18566569173581843, 1.8524875724573289,
      -0.18564537437166742, 1.8524255210534533, -0.18564538090231686, 1.8523427858482646,
      -0.18564538854050341, 1.8522600498960486, -0.18563523073199883, 1.852177314258372,
      -0.1856250717014105, 1.852084235855293, -0.18558441918638036, 1.8520118425508667,
      -0.18555392885326194, 1.8519497911080103, -0.18548278054767153, 1.851908424134568,
      -0.1854014672210403, 1.8518773997294453, -0.1853404820287538, 1.8518463760311155,
      -0.18527949666729704, 1.8518050118457905, -0.1852286752531498, 1.8517739901487216,
      -0.18513719683245522, 1.8517326282580837, -0.1850660464477075, 1.8516705858039069,
      -0.1850050592936848, 1.851598203195129, -0.18496439940759357, 1.851525821036894,
      -0.18493390280685654, 1.8514637786409567, -0.18492373446428506, 1.8514017364801523,
      -0.18491356543464932, 1.8513396935551185, -0.18491355993727504, 1.851277648396113,
      -0.18493388219148216, 1.851225939686831, -0.18498469761055747, 1.8511845707188266,
      -0.18503551382315714, 1.8511431968466974, -0.18511682239303032, 1.8511328476108682,
      -0.18517780637161896, 1.8511431803673817, -0.18523879304451857, 1.851163854589929,
      -0.18529978100878403, 1.8511948706319659, -0.185360770150041, 1.8512362276491314,
      -0.18543192451942578, 1.8512775880031034, -0.1854827501407481, 1.85132929184779,
      -0.18552341223756463, 1.851360310101894, -0.18559456469755148, 1.8513913300837215,
      -0.18565555276513232, 1.8514326928101248, -0.18571654147556824, 1.851443030010248,
      -0.18577752762071573, 1.8514326828964742, -0.18583851218390054, 1.8513602782977292,
      -0.18592998414051692];
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(coos.len() >> 1);
    for i in (0..coos.len()).step_by(2) {
      vertices.push((coos[i], coos[i + 1]))
    }
    // let expected_res_exact: [u64; 9] =...
    // Problem was an error in Newtom-Raphson method
    let _actual_res_exact = polygon_coverage(depth, &vertices, true);
    // So far we are happy if it does not panic :)
    // to_aladin_moc(&actual_res_exact);
  }

  #[test]
  fn testok_zone_1() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
        0.0_f64.to_radians(), 0.0_f64.to_radians(),
        10.0_f64.to_radians(), 10.0_f64.to_radians(),
      );
    let expected_res_exact: [u64; 53] =[4351, 4515, 4516, 4517, 4518, 4519, 4520, 4521, 4522, 4523,
      4524, 4525, 4526, 4527, 4528, 4530, 4531, 4536, 4537, 4538, 4539, 4540, 4542, 4543, 4693,
      4864, 4865, 4866, 4867, 4868, 4869, 4870, 4871, 4873, 4876, 4877, 4878, 4879, 4880, 4881,
      4882, 4883, 4884, 4885, 4886, 4887, 4888, 4889, 4890, 4891, 4892, 4901, 4912];
    let actual_res_exact = zone_coverage(depth, lon_min, lat_min, lon_max, lat_max);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_zone_2() {
    let depth = 2;
    let (lon_min, lat_min, lon_max, lat_max) = (
      284.2_f64.to_radians(), -25.0_f64.to_radians(),
      26.7_f64.to_radians(), 85.0_f64.to_radians(),
    );
    let expected_res_exact: [u64; 53] =[2, 8, 9, 10, 11, 14, 15, 48, 49, 50, 51, 52, 53, 54, 55, 56,
      57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
      113, 116, 117, 118, 119, 125, 139, 142, 183, 187, 188, 189, 190, 191];
    let actual_res_exact = zone_coverage(depth, lon_min, lat_min, lon_max, lat_max);
    to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_zone_3() {
    let depth = 4;
    let (lon_min, lat_min, lon_max, lat_max) = (
      12.7_f64.to_radians(), 64.5_f64.to_radians(),
      350.2_f64.to_radians(), 65.0_f64.to_radians(),
    );
    let expected_res_exact: [u64; 65] =[127, 207, 211, 212, 213, 214, 216, 217, 218, 227, 228, 229,
      230, 232, 233, 383, 447, 463, 467, 468, 469, 470, 472, 473, 474, 483, 484, 485, 486, 488, 489,
      490, 639, 703, 719, 723, 724, 725, 726, 728, 729, 730, 739, 740, 741, 742, 744, 745, 746, 959,
      975, 979, 980, 981, 982, 984, 985, 986, 995, 996, 997, 998 ,1000, 1001, 1002];
    let actual_res_exact = zone_coverage(depth, lon_min, lat_min, lon_max, lat_max);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }

  #[test]
  fn testok_zone_4() {
    let depth = 5;
    let (lon_min, lat_min, lon_max, lat_max) = (
      75.7_f64.to_radians(), -75.0_f64.to_radians(),
      12.2_f64.to_radians(), -74.5_f64.to_radians(),
    );
    let expected_res_exact: [u64; 71] =[8257, 8258, 8259, 8260, 8321, 8322, 8323, 8328, 9245, 9246,
      9247, 9261, 9262, 9263, 9265, 9266, 9267, 9268, 9272, 9281, 9282, 9283, 9284, 9288, 9345,
      9346, 9347, 9348, 9352, 10269, 10270, 10271, 10285, 10286, 10287, 10289, 10290, 10291, 10292,
      10296, 10305, 10306, 10307, 10308, 10312, 10369, 10370, 10371, 10372, 10376, 11293, 11294,
      11295, 11309, 11310, 11311, 11313, 11314, 11315, 11316, 11320, 11329, 11330, 11331, 11332,
      11336, 11393, 11394, 11395, 11396, 11400];
    let actual_res_exact = zone_coverage(depth, lon_min, lat_min, lon_max, lat_max);
    // to_aladin_moc(&actual_res_exact);
    assert!(actual_res_exact.deep_size() > 0);
    assert_eq!(expected_res_exact.len(), actual_res_exact.deep_size());
    for (h1, h2) in actual_res_exact.flat_iter().zip(expected_res_exact.iter()) {
      assert_eq!(h1, *h2);
    }
  }


  fn to_radians(lonlats: &mut [(f64, f64)]) {
    for (lon, lat) in lonlats.iter_mut() {
      *lon = lon.to_radians();
      *lat = lat.to_radians();
    }
  }

  #[allow(dead_code)]
  fn to_degrees(lonlats: &mut [(f64, f64)]) {
    for (lon, lat) in lonlats.iter_mut() {
      *lon = lon.to_degrees();
      *lat = lat.to_degrees();
    }
  }

  #[test]
  fn testok_bmoc_not() {
    let actual_res = cone_coverage_approx_custom(3, 4, 36.80105218_f64.to_radians(), 56.78028536_f64.to_radians(), 14.93_f64.to_radians());
    /*println!("@@@@@ HIERARCH VIEW");
    for cell in actual_res.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }
    println!("@@@@@ HIERARCH VIEW, NOT");*/
    let complement: BMOC = actual_res.not();
    /*for cell in complement.into_iter() {
      println!("@@@@@ cell a: {:?}", cell);
    }*/
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
    let layer_0 = get(0);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    let layer_0 = get(0);
    assert_eq!(1, layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_prec_2() {
    let lon_deg = 359.99999999999994_f64;
    let lat_deg = 41.81031502783791_f64;
    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_1 = get(1);
    layer_1.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_1.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);*/
    let layer_1 = get(1);
    assert_eq!(13, layer_1.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  #[test]
  fn test_prec_3() {
    let lon_deg = 359.99999999999994_f64; //359.99999999999994_f64;
    let lat_deg = 41.81031489577861_f64; //41.81031489577857_f64;

    /*let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get(6);
    layer_0.shift_rotate_scale(&mut xy);
    println!("x: {}, y: {}", xy.0, xy.1);
    let mut ij = discretize(xy);
    println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1/*, &mut ij, xy, lon, lat*/);
    println!("d0h_bits: {}", d0h_bits);*/
    println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));*/
    
    let layer_6 = get(6);
    assert_eq!(13653, layer_6.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }

  /*
  #[test]
  fn test_prec_4() {
    let lon_deg = 292.49999999999994_f64; //359.99999999999994_f64;
    let lat_deg = 41.810314895778546_f64; //41.81031489577857_f64;

    let mut xy = proj(lon_deg.to_radians(), lat_deg.to_radians());
    // println!("proj_x: {:.17}, proj_y: {:.17}", xy.0, xy.1);
    xy.0 = ensures_x_is_positive(xy.0);
    let layer_0 = get(0);
    layer_0.shift_rotate_scale(&mut xy);
    // println!("x: {}, y: {}", xy.0, xy.1);
    let ij = discretize(xy);
    // println!("i: {}, j: {}", ij.0, ij.1);
    let ij_d0c = layer_0.base_cell_coos(&ij);
    // println!("i0: {}, j0: {}", ij_d0c.0, ij_d0c.1);/*
    let d0h_bits = layer_0.depth0_bits(ij_d0c.0, ij_d0c.1, ij, xy/*, lon, lat*/);
    // println!("d0h_bits: {}", d0h_bits);*/
    // println!("hash: {}", layer_0.hash(lon_deg.to_radians(), lat_deg.to_radians()));

    let layer_2 = get(2);
    assert_eq!(56, layer_2.hash(lon_deg.to_radians(), lat_deg.to_radians()));
  }
   */

  #[test]
  fn test_ring() {
    for depth in 0..10 {
      let layer = get(depth);
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
    /*println!("{:?}", res);
    // Result with previous version of hash_dxdy_v1 !
    assert_eq!(res, [
      (20, 0.0), 
      (38, 0.1661686383097217), 
      (33, 0.2024027885319438), 
      (20, 0.6314285731583344)
    ]);*/
    // Result with previous version of hash_dxdy_v2 !
    assert_eq!(res, [
      (20, 0.0),
      (38, 0.1661686383097213),
      (33, 0.20240278853194385),
      (20, 0.6314285731583348)
    ]);
  }
  
  #[test]
  fn test_bilinear_interpolation_2() {
    let lon_deg = 83.633478_f64;
    let lat_deg = 22.015110_f64;
    let res = bilinear_interpolation(18, lon_deg.to_radians(), lat_deg.to_radians());
    /* println!("{:?}", res);
    // Result with previous version of hash_dxdy_v1 !
    assert_eq!(res, [
      (405766747916, 0.5757471135241182), 
      (405766747917, 0.3604806280107034), 
      (405766747918, 0.039217694696856834), 
      (405766747919, 0.024554563768321474)
    ]);
    */
    assert_eq!(res, [
      (405766747916, 0.5757471135599139),
      (405766747917, 0.3604806280331154),
      (405766747918, 0.039217694661061175),
      (405766747919, 0.024554563745909478)
    ]);
  }

  #[test]
  fn test_bilinear_interpolation_3() {
    let lon_rad = [0.17453293_f64, 0.43633231_f64, 0.0_f64];
    let lat_rad = [0.08726646_f64, 0.17453293_f64, 0.78539816_f64];
    let depth = 5;
    for (lon, lat) in lon_rad.iter().zip(lat_rad.iter()) {
      let res = bilinear_interpolation(depth, *lon, *lat);
      println!("{:?}", res);
    }
  }


  #[test]
  fn test_gen_file() -> std::io::Result<()> {
    use std::io::prelude::*;
    
    let depth = 8;
    let layer = get(depth);
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

  /*#[test]
  fn test_xpm1_and_q () {
    let a = Layer::d0h_lh_in_d0c(std::f64::NAN, 0.0);
    eprintln!("{:?}", &a);
    let a = Layer::d0h_lh_in_d0c(std::f64::INFINITY, 0.0);
    eprintln!("{:?}", &a);
  }*/

}
