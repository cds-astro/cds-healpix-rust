use super::*;

pub const fn n_hash(nside: u32) -> u64 {
  let nside = nside as u64;
  12 * nside * nside
}

/// Returns the number of isolatitude in the whole sphere at the given `nside`, 
/// i.e. the number of small circles parallel to the equator containing HEALPix cell centers.
/// 
/// # Input
/// - `nside`: the sub-division of a base cell along both the South-East and the South-West axis
/// 
/// # Output
/// -  $4 * nside - 1$: $2 * nside - 1$ ring for the NPC cell, $2 * nside - 1$ for the SPC cell
/// $+1$ for the equator.
/// 
/// # Remark
/// Contrary to the NESTED scheme, `nside` here is not necessarily a power of 2.
/// The smallest possible value is `1`.
/// 
/// # Example
/// ```rust
/// use cdshealpix::ring::{n_isolatitude_rings};
/// 
/// assert_eq!(n_isolatitude_rings(1),  3);
/// assert_eq!(n_isolatitude_rings(2),  7);
/// assert_eq!(n_isolatitude_rings(3), 11);
/// assert_eq!(n_isolatitude_rings(4), 15);
/// assert_eq!(n_isolatitude_rings(5), 19);
/// assert_eq!(n_isolatitude_rings(6), 23);
/// ```
pub const fn n_isolatitude_rings(nside: u32) -> u32 {
  // Not yet stable in const fn: debug_assert!(nside > 0);
  (nside << 2) - 1
}

/// Index of the first cell which is fully in the Equatorial Region, 
/// i.e. number of cells in the 4 polar cap triangles, 
/// i.e. four time the $nside^{\mathrm{th}}$ [triangular number](https://en.wikipedia.org/wiki/Triangular_number).
/// ```math
/// 4 * \sum_{i=1}^{nside} i = 4 * \frac{nside (nside + 1)}{2} = 2 nside (nside + 1)
/// ```
/// 
/// # Warning
/// To obtain the index of the first cell on the NPC/EQR transition latitude, call the
/// `first_hash_on_npc_eqr_transition` method.
/// 
/// # Example
/// ```rust
/// use cdshealpix::ring::first_hash_in_eqr;
/// 
/// assert_eq!(first_hash_in_eqr(1),  4);
/// assert_eq!(first_hash_in_eqr(2), 12);
/// assert_eq!(first_hash_in_eqr(3), 24);
/// assert_eq!(first_hash_in_eqr(4), 40);
/// assert_eq!(first_hash_in_eqr(5), 60);
/// assert_eq!(first_hash_in_eqr(6), 84);
/// assert_eq!(first_hash_in_eqr(7), 112);
/// assert_eq!(first_hash_in_eqr(8), 144);
/// ```
#[inline]
pub const fn first_hash_in_eqr(nside: u32) -> u64 {
  // Not yet stable in const fn: debug_assert!(nside > 0);
  triangular_number_x4(nside as u64)
}

/// Index of the first cell on the North polar cap / Equatorial Region transition latitude. 
/// I.e. number of cells in the 4 polar cap triangles of side = nside - 1, 
/// I.e. four time the $(nside - 1)^{\mathrm{th}}$ [triangular number](https://en.wikipedia.org/wiki/Triangular_number).
/// ```math
/// 4 * \sum_{i=1}^{nside - 1} i = 4 * \frac{nside (nside - 1)}{2} = 2 nside (nside - 1)
/// ```
/// 
/// # Warning
/// To obtain the index of the first cell fully in the EQR, call the
/// `first_hash_in_eqr` method.
/// 
/// # Example
/// ```rust
/// use cdshealpix::ring::first_hash_on_npc_eqr_transition;
/// 
/// assert_eq!(first_hash_on_npc_eqr_transition(1),  0);
/// assert_eq!(first_hash_on_npc_eqr_transition(2),  4);
/// assert_eq!(first_hash_on_npc_eqr_transition(3), 12);
/// assert_eq!(first_hash_on_npc_eqr_transition(4), 24);
/// assert_eq!(first_hash_on_npc_eqr_transition(5), 40);
/// assert_eq!(first_hash_on_npc_eqr_transition(6), 60);
/// assert_eq!(first_hash_on_npc_eqr_transition(7), 84);
/// assert_eq!(first_hash_on_npc_eqr_transition(8), 112);
/// assert_eq!(first_hash_on_npc_eqr_transition(9), 144);
/// ```
#[inline]
pub const fn first_hash_on_npc_eqr_transition(nside: u32) -> u64 {
  // Not yet stable in const fn: debug_assert!(nside > 0);
  triangular_number_x4((nside - 1) as u64)
}

/// Index of the first cell on the Equatorial Region / South polar cap transition latitude.
/// # Warning
/// To obtain the index of the first cell fully in the SPC, call the
/// `first_hash_in_spc` method.
/// 
/// # Example
/// ```rust
/// use cdshealpix::ring::first_hash_on_eqr_spc_transition;
/// 
/// assert_eq!(first_hash_on_eqr_spc_transition(1),   8);
/// assert_eq!(first_hash_on_eqr_spc_transition(2),  36);
/// assert_eq!(first_hash_on_eqr_spc_transition(4), 152);
/// ```
#[inline]
pub const fn first_hash_on_eqr_spc_transition(nside: u32) -> u64 {
  let n = nside as u64;
  // 8n^2 + (2*n^2 - 2n) = 2n(5n - 1)
  (n * (5 * n - 1)) << 1
}

/// Index of the first cell fully in the South polar cap.
/// # Warning
/// To obtain the index of the first cell on the EQR/SPC transition latitude, call the
/// `first_hash_on_eqr_spc_transition` method.
/// 
/// # Example
/// ```rust
/// use cdshealpix::ring::{first_hash_in_spc};
/// 
/// assert_eq!(first_hash_in_spc(1),  12);
/// assert_eq!(first_hash_in_spc(2),  44);
/// assert_eq!(first_hash_in_spc(4), 168);
/// ```
#[inline]
pub const fn first_hash_in_spc(nside: u32) -> u64 {
  let n = nside as u64;
  // 8n^2 + (2*n^2 + 2n) = 2n(5n + 1)
  (n * (5 * n + 1)) << 1
}

/// Four time the [triangular number](https://en.wikipedia.org/wiki/Triangular_number), i.e.
/// ```math
/// 4 * \sum_{i=1}^{n} i = 4 * \frac{n (n + 1)}{2} = 2 n (n + 1)
/// ```
#[inline]
pub(crate) const fn triangular_number_x4(n: u64) -> u64 {
  (n * (n + 1)) << 1
}


/*pub fn hash_v2(nside: u32, lon: f64, lat: f64) -> u64 {
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
  
  if lat.sign == F64_SIGN_BIT_MASK {
    n_hash(nside_u32) - val
  } else {
    val 
  }
  apply_offset_and_signs(&mut xy, x.offset, lon.sign, lat.sign);
}*/

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
pub fn hash(nside: u32, lon: f64, lat: f64) -> u64 {
  let nside = nside as u64;
  let half_nside = 0.5 * (nside as f64);
  // Project and ensure that x is positive (in case input lon < 0)
  let mut xy = proj(lon, lat);
  xy.0 = ensures_x_is_positive(xy.0);
  // Set the origin to (x=0, y=-3) in the projection plane, scale both axis by a factor  1/2 * nside
  let mut dl = half_nside * xy.0;          debug_assert!(0.0 <= dl && dl <  4.0 * nside as f64);
  let mut dh = half_nside * (xy.1 + 3.0);  debug_assert!(0.0 <= dh && dh <= 2.5 * nside as f64);
  let mut i_ring = (dh as u64) << 1;
  let mut i_in_ring = dl as u64;
  dl -= i_in_ring as f64;      debug_assert!(0.0 <= dl && dl < 1.0);
  dh -= (i_ring >> 1) as f64;  debug_assert!(0.0 <= dh && dh < 1.0);
  deal_with_1x1_box(dl, dh, &mut i_ring, &mut i_in_ring);
  // Those tests are already performed in the proj, so if we have to improve performances we may
  // merge this code with the projection code (loosing readability).
  if i_ring >= 5 * nside { // North pole, rare case
    return i_in_ring / nside;
  }
  i_ring = 5 * nside - 1 - i_ring;
  if i_ring < nside { // North polar cap
    let off = nside - 1 - i_ring;
    i_in_ring -= (off >> 1) + (off & 1) + off * (i_in_ring / nside);
    triangular_number_x4(i_ring) + i_in_ring
  } else if i_ring >= 3 * nside { // South polar cap
    let off = i_ring + 1 - 3 * nside;
    i_in_ring -= (off >> 1) + (off & 1) + off * (i_in_ring / nside);
    n_hash(nside as u32) - triangular_number_x4(n_isolatitude_rings(nside as u32) as u64 - i_ring) + i_in_ring
  } else { // Equatorial region
    first_hash_in_eqr(nside as u32)
      + (i_ring - nside) * (nside << 2) 
      + if i_in_ring == nside << 2 { 0 } else { i_in_ring }
  }
}

///    NORTH   
///  h ^           Deals with a box of size 1x1, with: 
/// W  |_._    E   - dl in `[0, 1[` along the x-axis
/// E  |\2/|   A   - dh in `[0, 1[` along the y-axis
/// S  .3X1.   S   The box is divided in 4 quarters (q in [0, 3[) 
/// T  |/0\|   T   - q = 0 => nothing to do
///    --.----> l  - q = 1 => increase `i_ring` of  +1 and increase `i_in_ring` of +1
///                - q = 2 => increase `i_ring` of  +2
///    SOUTH       - q = 3 => increase `i_ring` of  +1
/// Warning, the inequalities are important to decide if an edge is in a cell or its neighbour!!
fn deal_with_1x1_box(dl: f64, dh: f64, i_ring: &mut u64, i_in_ring: &mut u64) {
  debug_assert!(0.0 <= dl && dl < 1.0);
  debug_assert!(0.0 <= dh && dh < 1.0);
  let in1 = (dl <= dh) as u64;       /* 1/0 */  debug_assert!(in1 == 0 || in1 == 1);
  let in2 = (dl >= 1.0 - dh) as u64; /* 0\1 */  debug_assert!(in2 == 0 || in2 == 1);
  *i_ring += in1 + in2;
  *i_in_ring += in2 >> in1; // <=> in2 & (1 - in1)
}

#[inline]
fn is_hash(nside: u32, hash: u64) -> bool {
  hash < n_hash(nside)
}

#[inline]
fn check_hash(nside: u32, hash: u64) { 
  assert!(is_hash(nside, hash), "Wrong hash value: too large.");
}

/// Center of the given cell in the Euclidean projection space.
/// # Output
/// - `(x, y)` coordinates such that $x \in [0, 8[$ and $y \in [-2, 2]$.
/// 
/// # TODO
/// TO BE TESTED!!
pub fn  center_of_projected_cell(nside: u32, hash: u64) -> (f64, f64) {
  check_hash(nside, hash);
  if hash < first_hash_on_npc_eqr_transition(nside) { // North polar cap
    let i_ring = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
    let n_in_ring = i_ring + 1;
    let i_in_ring = hash - triangular_number_x4(i_ring);
    let y = 1.0 + (nside as u64 - 1 - i_ring) as f64 / (nside as f64);
    let q = i_in_ring / n_in_ring;
    let x = ((i_in_ring + q * (nside as u64 - n_in_ring)) << 1) - i_ring;
    (x as f64 / nside as f64, y)
  } else if hash >= first_hash_in_spc(nside) { // South polar cap
    let hash = n_hash(nside) - 1 - hash; // start counting in reverse order from south polar cap
    let i_ring = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
    let n_in_ring = i_ring + 1;
    let i_in_ring = ((n_in_ring << 2) - 1) - (hash - triangular_number_x4(i_ring));
    let y = 1.0 + (nside as u64 - 1 - i_ring) as f64 / (nside as f64);
    let q = i_in_ring / n_in_ring;
    let x = ((i_in_ring + q * (nside as u64 - n_in_ring)) << 1) - i_ring;
    (x as f64 / nside as f64, -y)
  } else { // Equatorial region
    let i_ring = (hash - first_hash_on_npc_eqr_transition(nside)) / (nside as u64);
    let y = (nside as i64 - i_ring as i64) as f64 / (nside as f64);
    let x = (((hash - i_ring) << 1) + ((i_ring + 1) & 1)) as f64 / (nside as f64);
    (x, y)
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
pub fn center(nside: u32, hash: u64) -> (f64, f64) {
  let (x, y) = center_of_projected_cell(nside, hash);
  super::unproj(x, y)
}























#[test]
fn test_hash() {
  let nside = 1;
  // Easy
  assert_eq!(hash(nside,  40.0_f64.to_radians(),  45.0_f64.to_radians()), 0);
  assert_eq!(hash(nside, 130.0_f64.to_radians(),  45.0_f64.to_radians()), 1);
  assert_eq!(hash(nside, 220.0_f64.to_radians(),  45.0_f64.to_radians()), 2);
  assert_eq!(hash(nside, 310.0_f64.to_radians(),  45.0_f64.to_radians()), 3);
  assert_eq!(hash(nside,   0.0_f64.to_radians(),   0.0_f64.to_radians()), 4);
  assert_eq!(hash(nside,  90.0_f64.to_radians(),   0.0_f64.to_radians()), 5);
  assert_eq!(hash(nside, 180.0_f64.to_radians(),   0.0_f64.to_radians()), 6);
  assert_eq!(hash(nside, 270.0_f64.to_radians(),   0.0_f64.to_radians()), 7);
  assert_eq!(hash(nside,  40.0_f64.to_radians(), -45.0_f64.to_radians()), 8);
  assert_eq!(hash(nside, 130.0_f64.to_radians(), -45.0_f64.to_radians()), 9);
  assert_eq!(hash(nside, 220.0_f64.to_radians(), -45.0_f64.to_radians()), 10);
  assert_eq!(hash(nside, 310.0_f64.to_radians(), -45.0_f64.to_radians()), 11);
  // Particular points
  // - north pole
  assert_eq!(hash(nside,   40.0_f64.to_radians(),  90.0_f64.to_radians()), 0);
  assert_eq!(hash(nside,  130.0_f64.to_radians(),  90.0_f64.to_radians()), 1);
  assert_eq!(hash(nside,  220.0_f64.to_radians(),  90.0_f64.to_radians()), 2);
  assert_eq!(hash(nside,  310.0_f64.to_radians(),  90.0_f64.to_radians()), 3);
  // - south pole
  assert_eq!(hash(nside,   40.0_f64.to_radians(),  -90.0_f64.to_radians()), 8);
  assert_eq!(hash(nside,  130.0_f64.to_radians(),  -90.0_f64.to_radians()), 9);
  assert_eq!(hash(nside,  220.0_f64.to_radians(),  -90.0_f64.to_radians()), 10);
  assert_eq!(hash(nside,  310.0_f64.to_radians(),  -90.0_f64.to_radians()), 11);
  // - around corner
  assert_eq!(hash(nside, 44.0_f64.to_radians(),  0.0_f64.to_radians()), 4);
  assert_eq!(hash(nside, 46.0_f64.to_radians(),  0.0_f64.to_radians()), 5);
  assert_eq!(hash(nside, 45.0_f64.to_radians(),  1.0_f64.to_radians()), 0);
  assert_eq!(hash(nside, 45.0_f64.to_radians(),  -1.0_f64.to_radians()), 8);
}

#[test]
fn test_hash_2() {
  for depth in 0..10 {
    let nside = nside(depth);
    let layer = nested::get_or_create(depth);
    for icell in 0..layer.n_hash() {
      let (lon, lat) = layer.center(icell);
      assert_eq!(hash(nside, lon, lat), layer.to_ring(icell));
    }
  }   
}


#[test]
fn test_deal_with_1x1_box() {
  let mut i_ring = 0;
  let mut i_in_ring = 0;
  
  deal_with_1x1_box(0.0, 0.0, &mut i_ring, &mut i_in_ring);
  assert_eq!(i_ring, 1);
  assert_eq!(i_in_ring, 0);
  
  i_ring = 0; i_in_ring = 0;
  deal_with_1x1_box(0.5, 0.0, &mut i_ring, &mut i_in_ring);
  assert_eq!(i_ring, 0);
  assert_eq!(i_in_ring, 0);
  
  i_ring = 0; i_in_ring = 0;
  deal_with_1x1_box(0.9, 0.5, &mut i_ring, &mut i_in_ring);
  assert_eq!(i_ring, 1);
  assert_eq!(i_in_ring, 1);
  
  i_ring = 0; i_in_ring = 0;
  deal_with_1x1_box(0.5, 0.9, &mut i_ring, &mut i_in_ring);
  assert_eq!(i_ring, 2);
  assert_eq!(i_in_ring, 0);

  i_ring = 0; i_in_ring = 0;
  deal_with_1x1_box(0.5, 0.5, &mut i_ring, &mut i_in_ring);
  assert_eq!(i_ring, 2);
  assert_eq!(i_in_ring, 0);
}

/*
to_nested() {
  check firs that the nside is a power of 2!!
}
*/

// put also here the Cone query returning ranges!!