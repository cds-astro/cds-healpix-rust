use super::*;

pub const fn nside_max() -> u64 {
  1 << 29
}

/// Returns the number of distinct hash values (i.e. the number if cells the sphere is divided in)
/// at the gicen NSIDE.
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

// TODO: remove the previous method and the '_u32' below
/// Four time the [triangular number](https://en.wikipedia.org/wiki/Triangular_number), i.e.
/// ```math
/// 4 * \sum_{i=1}^{n} i = 4 * \frac{n (n + 1)}{2} = 2 n (n + 1)
/// ```
#[inline]
pub(crate) const fn triangular_number_x4_u32(n: u32) -> u64 {
  let n = n as u64;
  (n * (n + 1)) << 1
}

/// Returns the cell number (hash value) associated with the given position on the unit sphere
/// # Inputs
/// - `nside`: the NSIDE of the RING scheme
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
/// use cdshealpix::ring::hash;
/// 
/// let nside = 1;
/// // Easy
/// assert_eq!(hash(nside,  40.0_f64.to_radians(),  45.0_f64.to_radians()),  0);
/// assert_eq!(hash(nside, 130.0_f64.to_radians(),  45.0_f64.to_radians()),  1);
/// assert_eq!(hash(nside, 220.0_f64.to_radians(),  45.0_f64.to_radians()),  2);
/// assert_eq!(hash(nside, 310.0_f64.to_radians(),  45.0_f64.to_radians()),  3);
/// assert_eq!(hash(nside,   0.0_f64.to_radians(),   0.0_f64.to_radians()),  4);
/// assert_eq!(hash(nside,  90.0_f64.to_radians(),   0.0_f64.to_radians()),  5);
/// assert_eq!(hash(nside, 180.0_f64.to_radians(),   0.0_f64.to_radians()),  6);
/// assert_eq!(hash(nside, 270.0_f64.to_radians(),   0.0_f64.to_radians()),  7);
/// assert_eq!(hash(nside,  40.0_f64.to_radians(), -45.0_f64.to_radians()),  8);
/// assert_eq!(hash(nside, 130.0_f64.to_radians(), -45.0_f64.to_radians()),  9);
/// assert_eq!(hash(nside, 220.0_f64.to_radians(), -45.0_f64.to_radians()), 10);
/// assert_eq!(hash(nside, 310.0_f64.to_radians(), -45.0_f64.to_radians()), 11);
/// // Particular points
/// // - north pole
/// assert_eq!(hash(nside,  40.0_f64.to_radians(), 90.0_f64.to_radians()), 0);
/// assert_eq!(hash(nside, 130.0_f64.to_radians(), 90.0_f64.to_radians()), 1);
/// assert_eq!(hash(nside, 220.0_f64.to_radians(), 90.0_f64.to_radians()), 2);
/// assert_eq!(hash(nside, 310.0_f64.to_radians(), 90.0_f64.to_radians()), 3);
/// // - south pole
/// assert_eq!(hash(nside,  40.0_f64.to_radians(), -90.0_f64.to_radians()),  8);
/// assert_eq!(hash(nside, 130.0_f64.to_radians(), -90.0_f64.to_radians()),  9);
/// assert_eq!(hash(nside, 220.0_f64.to_radians(), -90.0_f64.to_radians()), 10);
/// assert_eq!(hash(nside, 310.0_f64.to_radians(), -90.0_f64.to_radians()), 11);
/// // - around the (45, 0) corner
/// assert_eq!(hash(nside, 44.0_f64.to_radians(),  0.0_f64.to_radians()), 4);
/// assert_eq!(hash(nside, 46.0_f64.to_radians(),  0.0_f64.to_radians()), 5);
/// assert_eq!(hash(nside, 45.0_f64.to_radians(),  1.0_f64.to_radians()), 0);
/// assert_eq!(hash(nside, 45.0_f64.to_radians(), -1.0_f64.to_radians()), 8);
/// ```
pub fn hash(nside: u32, lon: f64, lat: f64) -> u64 {
  hash_with_dldh(nside, lon, lat).0
}

/// Returns the cell number (hash value) associated with the given position on the unit sphere, 
/// together with the offset `(dx, dy)` on the Euclidean plane of the projected position with
/// respect to the origin of the cell (South vertex).
/// # Inputs
/// - `nside`: the NSIDE of the RING scheme
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
/// TODO
pub fn hash_with_dxdy(nside: u32, lon: f64, lat: f64) -> (u64, f64, f64) {
  let (h, dl, dh) = hash_with_dldh(nside, lon, lat);
  let (dx, dy) = dldh_to_dxdy(dl, dh);
  (h, dx, dy)
}

fn hash_with_dldh(nside: u32, lon: f64, lat: f64) -> (u64, f64, f64) {
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
    return (i_in_ring / nside, 1.0, 1.0);
  }
  i_ring = 5 * nside - 1 - i_ring;
  let hash = if i_ring < nside { // North polar cap
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
  };
  (hash, dl, dh)
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

/// Build the hash value from the ring index (starting at the North pole, increasing southward)
/// and the index inside the ring (starting at the West, increasing eastward).
/// WARNING: the index in the ring in polar caps is provided like the index in ring in the euqatorial
/// region, i.e. like if the projection rectangle was full (no missing triangles).
// TODO: REMOVE DUPLICATE CODE IN RING HASH() AND NESTED TO_RING() 
fn build_hash(nside: u32, i_ring: u32, mut i_in_ring: u32) -> u64 {
  debug_assert!(0 <= i_ring && i_ring < n_isolatitude_rings(nside));
  debug_assert!(0 <= i_in_ring && i_in_ring <= n_isolatitude_rings(nside));
  if i_ring < nside { // North polar cap
    let off = nside - 1 - i_ring;
    i_in_ring -= (off >> 1) + (off & 1) + off * (i_in_ring / nside);
    triangular_number_x4_u32(i_ring) + i_in_ring as u64
  } else if i_ring >= 3 * nside { // South polar cap
    let off = i_ring + 1 - 3 * nside;
    i_in_ring -= (off >> 1) + (off & 1) + off * (i_in_ring / nside);
    n_hash(nside) - triangular_number_x4_u32(n_isolatitude_rings(nside) - i_ring) + i_in_ring as u64
  } else { // Equatorial region
    first_hash_in_eqr(nside)
      + (i_ring - nside) as u64 * (nside << 2) as u64
      + if i_in_ring == nside << 2 { 0 } else { i_in_ring } as u64
  }
}

/*
// Return i_ring and i_in_ring
fn decode_hash(nside: u32, hash: u64) -> (u32, u32) {

}*/

/*
/// Build the hash value from the ring index (starting at the North pole, increasing southward)
/// and the index inside the ring (starting at the West, increasing eastward)
/// ** counted in the first quarter, i.e. in [0, pi/2]**.
/// WARNING: the index in the ring in polar caps is provided like the index in ring in the euqatorial
/// region, i.e. like if the projection rectangle was full (no missing triangles).
/// - `q` simply equals i_in_ring / nside
/// - `q` = 0 if $lon \in [0, pi/2[$
/// - `q` = 1 if $lon \in [pi/2, pi[$
/// - `q` = 2 if $lon \in [pi, 3pi/2[$
/// - `q` = 3 if $lon \in [3pi/2, 2pi[$
/// Simply equals 
fn build_hash(nside: u32, i_ring: u32, mut i_in_ring: u32, q: u32) -> u64 {
  debug_assert!(0 <= i_ring && i_ring < n_isolatitude_rings(nside));
  debug_assert!(0 <= i_ring && i_ring < nside);
  if i_ring < nside { // North polar cap
    let off = nside - 1 - i_ring;
    i_in_ring -= (off >> 1) + (off & 1) + (off - nside) * q;
    triangular_number_x4_u32(i_ring) as u64 + i_in_ring as u64
  } else if i_ring >= 3 * nside { // South polar cap
    let off = i_ring + 1 - 3 * nside;
    i_in_ring -= (off >> 1) + (off & 1) + (off - nside) * q;
    n_hash(nside) - triangular_number_x4_u32(n_isolatitude_rings(nside) - i_ring) + i_in_ring as u64
  } else { // Equatorial region
    i_in_ring += q * nside;
      first_hash_in_eqr(nside)
      + (i_ring - nside) as u64 * (nside << 2) as u64
      + if i_in_ring == nside << 2 { 0 } else { i_in_ring } as u64
  }
}*/

fn dldh_to_dxdy(dl: f64, dh: f64) -> (f64, f64) {
  /*let dl = dl - 0.5;
  let dh = dh - 0.5;
  let dx = dh + dl;
  let dy = dh - dl;
  q = 0 => dx + 1, dy + 1
  q = 1 => dx    , dy + 1
  q = 2 => dx, dy
  q = 3 => dx + 1, dy
  */
  let dx = dh + dl - 1.0;
  let dy = dh - dl;
  (dx + ((dx < 0.0) as u8) as f64, dy + ((dy < 0.0) as u8) as f64)
}


/* THE IDEE WOULD BE TO AVOID HAVING TO BLOCK OF IF CONDITION (ONE IN PROJ AND ONE AFTER PROJ)
fn hash_v2(nside: u32, lon: f64, lat: f64) -> u64 {
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

/*
THE IDEE WOULD BE TO USE THE 45deg ROTATION INSTEAD OF dela_with_1x1_box to compare perfs
fn hash_v2(nside: u32, lon: f64, lat: f64) -> u64
*/

#[inline]
fn is_hash(nside: u32, hash: u64) -> bool {
  hash < n_hash(nside)
}

#[inline]
fn check_hash(nside: u32, hash: u64) { 
  assert!(is_hash(nside, hash), "Wrong hash value: too large.");
}

/// Center of the given cell in the Euclidean projection space.
/// # Input
/// - `nside`: the NSIDE of the RING scheme
/// - `hash`: the hash value of the cell we look for the unprojected center
/// 
/// # Output
/// - `(x, y)` coordinates such that $x \in [0, 8[$ and $y \in [-2, 2]$.
pub fn  center_of_projected_cell(nside: u32, hash: u64) -> (f64, f64) {
  check_hash(nside, hash);
  if hash < first_hash_on_npc_eqr_transition(nside) { // North polar cap
    // Ring index from the northmost ring, increasing southward
    let i_ring = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
    // Number of cells in the ring for each base cell (Remark: n_in_ring = nside in th EQR)
    let n_in_ring = i_ring + 1;
    // Index in the ring
    let i_in_ring = hash - triangular_number_x4(i_ring);
    // Base cell containing the hash (0, 1, 2 or 3)
    let q = i_in_ring / n_in_ring;
    // Position of the center of the first cell on the ring cell, in [0, nside] <=> x in [0, 1], (i.e. inside a base cell) 
    let off_in_d0h = nside as u64 - i_ring;
    // Ring index in a base cell
    let i_in_d0h_ring = i_in_ring - q * n_in_ring;
    // x inside a base cell, between 0 and 2 * nside
    let x = ((i_in_d0h_ring) << 1) as f64 + off_in_d0h as f64;
    let y = 1.0 + (nside as u64 - 1 - i_ring) as f64 / (nside as f64);
    ((q << 1) as f64 + x as f64 / nside as f64, y)
  } else if hash >= first_hash_in_spc(nside) { // South polar cap
    let hash = n_hash(nside) - 1 - hash; // start counting in reverse order from south polar cap
    let i_ring = (((1 + (hash << 1)) as f64).sqrt() as u64 - 1) >> 1;
    let n_in_ring = i_ring + 1;
    let i_in_ring = ((n_in_ring << 2) - 1) - (hash - triangular_number_x4(i_ring));
    let q = i_in_ring / n_in_ring;
    let off_in_d0h = nside as u64 - i_ring;
    let i_in_d0h_ring = i_in_ring - q * n_in_ring;
    let x = ((i_in_d0h_ring) << 1) as f64 + off_in_d0h as f64;
    let y = 1.0 + (nside as u64 - 1 - i_ring) as f64 / (nside as f64);
    ((q << 1) as f64 + x as f64 / nside as f64, -y)
  } else { // Equatorial region
    let nsidex4 = (nside << 2) as u64;
    let i_ring = (hash - first_hash_on_npc_eqr_transition(nside)) / nsidex4;
    let i_in_ring = (hash - first_hash_on_npc_eqr_transition(nside)) - i_ring * nsidex4;
    let x = ((i_in_ring << 1) + ((i_ring + 1) & 1)) as f64 / (nside as f64);
    let y = (nside as i64 - i_ring as i64) as f64 / (nside as f64);
    (x, y)
  }
}

/// Compute the position on the unit sphere of the center (in the Euclidean projection plane)
/// of the cell associated to the given hash value.
/// 
/// # Input
/// - `nside`: the NSIDE of the RING scheme
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
/// TODO
///
pub fn center(nside: u32, hash: u64) -> (f64, f64) {
  let (x, y) = center_of_projected_cell(nside, hash);
  super::unproj(x, y)
}

/// Compute the position on the unit sphere of the position '(dx, dy)' from the south vertex of 
/// the HEALPix cell associated to the given hash value.
/// The x-axis is the South-East axis while the y-axis is the south-west axis.
/// 
/// # Input
/// - `nside`: the NSIDE of the RING scheme
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
/// TODO
///
pub fn sph_coo(nside: u32, hash: u64, dx: f64, dy: f64) -> (f64, f64) {
  assert!(0.0 <= dx && dx < 1.0);
  assert!(0.0 <= dy && dy < 1.0);
  let (mut x, mut y) = center_of_projected_cell(nside, hash);
  x += (dx - dy) / (nside as f64);
  y += (dx + dy - 1.0) / (nside as f64);
  super::unproj(ensures_x_is_positive(x), y)
}

/// Computes the positions on the unit sphere of the 4 vertices of the given cell.
/// If you want to access the position for a given direction, use method 
/// [vertices_map](#method.vertices_map).
///   
/// # Input
/// - `nside`: the NSIDE of the RING scheme
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
/// TODO
///
#[inline]
pub fn vertices(nside: u32, hash: u64) -> [(f64, f64); 4] {
  let one_over_nside = 1.0 / (nside as f64);
  let (x, y) = center_of_projected_cell(nside, hash);
  [
    super::unproj(x, y - one_over_nside), // S
    super::unproj(x + one_over_nside, y), // E
    super::unproj(x, y + one_over_nside), // N
    super::unproj(ensures_x_is_positive(x - one_over_nside), y)  // W
  ]
}

// TODO: implement path_along_cell_edge
// TODO: implement neighbours
// TODO: implement cone_coverage (will return RANGES) => need the exact cone impl in NESTED
// TODO: implement interpolation

// For polygon and elliptical_cone, I don't know yet how to do
// (we can't return MOC/BMOC since we have no guarantee that NSIDE is a power of 2)!!

#[cfg(test)]
mod tests {
  
  use crate::ring::*;
  
  #[test]
  fn test_hash() {
    // TODO: remove this (already in the doc) and test peciliar points!
    let nside = 1;
    // Easy
    assert_eq!(hash(nside,  40.0_f64.to_radians(),  45.0_f64.to_radians()),  0);
    assert_eq!(hash(nside, 130.0_f64.to_radians(),  45.0_f64.to_radians()),  1);
    assert_eq!(hash(nside, 220.0_f64.to_radians(),  45.0_f64.to_radians()),  2);
    assert_eq!(hash(nside, 310.0_f64.to_radians(),  45.0_f64.to_radians()),  3);
    assert_eq!(hash(nside,   0.0_f64.to_radians(),   0.0_f64.to_radians()),  4);
    assert_eq!(hash(nside,  90.0_f64.to_radians(),   0.0_f64.to_radians()),  5);
    assert_eq!(hash(nside, 180.0_f64.to_radians(),   0.0_f64.to_radians()),  6);
    assert_eq!(hash(nside, 270.0_f64.to_radians(),   0.0_f64.to_radians()),  7);
    assert_eq!(hash(nside,  40.0_f64.to_radians(), -45.0_f64.to_radians()),  8);
    assert_eq!(hash(nside, 130.0_f64.to_radians(), -45.0_f64.to_radians()),  9);
    assert_eq!(hash(nside, 220.0_f64.to_radians(), -45.0_f64.to_radians()), 10);
    assert_eq!(hash(nside, 310.0_f64.to_radians(), -45.0_f64.to_radians()), 11);
    // Particular points
    // - north pole
    assert_eq!(hash(nside,  40.0_f64.to_radians(), 90.0_f64.to_radians()), 0);
    assert_eq!(hash(nside, 130.0_f64.to_radians(), 90.0_f64.to_radians()), 1);
    assert_eq!(hash(nside, 220.0_f64.to_radians(), 90.0_f64.to_radians()), 2);
    assert_eq!(hash(nside, 310.0_f64.to_radians(), 90.0_f64.to_radians()), 3);
    // - south pole
    assert_eq!(hash(nside,  40.0_f64.to_radians(), -90.0_f64.to_radians()),  8);
    assert_eq!(hash(nside, 130.0_f64.to_radians(), -90.0_f64.to_radians()),  9);
    assert_eq!(hash(nside, 220.0_f64.to_radians(), -90.0_f64.to_radians()), 10);
    assert_eq!(hash(nside, 310.0_f64.to_radians(), -90.0_f64.to_radians()), 11);
    // - around the (45, 0) corner
    assert_eq!(hash(nside, 44.0_f64.to_radians(),  0.0_f64.to_radians()), 4);
    assert_eq!(hash(nside, 46.0_f64.to_radians(),  0.0_f64.to_radians()), 5);
    assert_eq!(hash(nside, 45.0_f64.to_radians(),  1.0_f64.to_radians()), 0);
    assert_eq!(hash(nside, 45.0_f64.to_radians(), -1.0_f64.to_radians()), 8);
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
  fn test_hash_3() {
    let nside = 4_u32;
    let lon = 4.71238898;
    let lat = -1.15965846;
    let h = hash(nside, lon, lat);
    println!("h: {}", &h);
  }
  
  #[test]
  fn test_deal_with_1x1_box() {
    let mut i_ring = 0;
    let mut i_in_ring = 0;

    deal_with_1x1_box(0.0, 0.0, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 1);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.5, 0.0, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 0);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.9, 0.5, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 1);
    assert_eq!(i_in_ring, 1);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.5, 0.9, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 2);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.5, 0.5, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 2);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.25, 0.75, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 2);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.75, 0.75, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 2);
    assert_eq!(i_in_ring, 0);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.75, 0.25, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 1);
    assert_eq!(i_in_ring, 1);

    i_ring = 0;
    i_in_ring = 0;
    deal_with_1x1_box(0.25, 0.25, &mut i_ring, &mut i_in_ring);
    assert_eq!(i_ring, 1);
    assert_eq!(i_in_ring, 0);
  }
  
  #[test]
  fn test_dldh_to_dxdy() {
    assert_eq!(dldh_to_dxdy(0.0, 0.0), (0.0, 0.0));
    assert_eq!(dldh_to_dxdy(0.5, 0.0), (0.5, 0.5));
    assert_eq!(dldh_to_dxdy(0.0, 0.5), (0.5, 0.5));
    assert_eq!(dldh_to_dxdy(0.5, 1.0), (0.5, 0.5));
    assert_eq!(dldh_to_dxdy(1.0, 0.5), (0.5, 0.5));

    assert_eq!(dldh_to_dxdy(0.25, 0.25), (0.5, 0.0));
    assert_eq!(dldh_to_dxdy(0.75, 0.25), (0.0, 0.5));
    assert_eq!(dldh_to_dxdy(0.75, 0.75), (0.5, 0.0));
    assert_eq!(dldh_to_dxdy(0.25, 0.75), (0.0, 0.5));
  }
  
  #[test]
  fn test_center() {
    let nside = 2;
    let ipix = 33;
    // let center = 
    println!("{:?}", center_of_projected_cell(nside, ipix));
    let (lon, lat) = center(nside, ipix);
    println!("(lon: {}, lat: {})", lon.to_degrees(), lat.to_degrees());
    /* dx 0.5 dy 0.5*/
  }
  
  #[test]
  fn test_center_2() {
    let nside = 2;
    assert_eq!(center_of_projected_cell(nside,  2), (5.0,  1.5));
    assert_eq!(center_of_projected_cell(nside, 46), (5.0, -1.5));

  }

  #[test]
  fn test_center_3() {
    let nside = 4;
    let l2 = nested::get_or_create(2);
    // NPC
    let ipix = 7;
    assert_eq!(
      center_of_projected_cell(nside, ipix),
      l2.center_of_projected_cell(l2.from_ring(ipix))
    );
    // SPC
    let ipix = 183;
    assert_eq!(
      center_of_projected_cell(nside, ipix),
      l2.center_of_projected_cell(l2.from_ring(ipix))
    );
    
    // let center = 
    //let (lon, lat) = center(nside, ipix);
    //println!("(lon: {}, lat: {})", lon.to_degrees(), lat.to_degrees());

    // println!("hash: {}", hash(nside, lon, lat));
  }
  
}
