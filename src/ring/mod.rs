use super::*;

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



/*
pub fn hash(nside: u32, lon: f64, lat: f64) -> u64 {
  let mut xy = proj(lon, lat);
  xy.0 = ensures_x_is_positive(xy.0);
  
  regarde si partie entier x * nside est 0 ou 1
  regarde si partie entier y * nside est 0 ou 1
  // Those tests are already performed in the proj, so if we have to improve performances we may
  // merge this code with the projection code (loosing readability).
  if y > 1.0 {
    
  } else if y < -1.0 {
    
  } else {
    let i_ring = ;
    let i_in_ring = ;
  }
    
  // projection
  // rotation, shift, scale
  // base cell computation
  // i, j ... ??
  0
}
*/

/*
to_nested() {
  check firs that the nside is a power of 2!!
}
*/

// put also here the Cone query returning ranges!!