use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_8, PI};

use super::sph_geom::coo3d::*;
use super::{Customf64, ONE_OVER_TRANSITION_Z, TRANSITION_Z};
// use super::sph_geom::{Polygon};

/// Returns the coordinates (lon, lat) of the "special points" the given great-circle arc contains
/// (if any). Most of the time the result will be empty or will contains a single value.
/// A list is returned for the (rare?) cases in which an great-circle arc overlap both the
/// equatorial region and a polar cap(s), and contains a special point in both.
/// # Inputs
/// - `p1` first point of the great-circle arc
/// - `p2` second point of the great_-circle arc
/// - `z_eps_max` the target precision on z (used to stop the Newton-Raphson method),
///   a reasonable choice may be `|p2.z - p1.z|/ 1000`.
///   Internally, the value can't be higher than `|z2 - z1| / 50`
///   nor lower than `1e-15`.
/// - `n_iter_max` upper limit on the number of iteration to be used in the Newton-Raphson method
///   a reasonable choice may be 20.
pub fn arc_special_points<'a>(
  mut p1: &'a Coo3D,
  mut p2: &'a Coo3D,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Box<[LonLat]> {
  // Ensure p1.z() < p2.z()
  if p1.z() > p2.z() {
    std::mem::swap(&mut p1, &mut p2);
  }
  if TRANSITION_Z <= p1.z() || p2.z() <= -TRANSITION_Z {
    // NPC only or SPC only
    match arc_special_point_in_pc(p1, p2, z_eps_max, n_iter_max) {
      Some(lonlat) => Box::new([lonlat; 1]),
      None => Vec::new().into_boxed_slice(),
    }
  } else if -TRANSITION_Z <= p1.z() && p2.z() <= TRANSITION_Z {
    // EQR only
    match arc_special_point_in_eqr(p1, p2, z_eps_max, n_iter_max) {
      Some(lonlat) => Box::new([lonlat; 1]),
      None => Vec::new().into_boxed_slice(),
    }
  } else if p1.z() < -TRANSITION_Z {
    // SPC, EQR (and maybe NPC)
    let mut res: Vec<LonLat> = Vec::with_capacity(3);
    let v_eqr_south = Coo3D::from(intersect_with_transition_lat_spc(p1, p2).unwrap());
    if let Some(lonlat) = arc_special_point_in_pc(p1, &v_eqr_south, z_eps_max, n_iter_max) {
      res.push(lonlat);
    }
    if p2.z() <= TRANSITION_Z {
      if let Some(lonlat) = arc_special_point_in_eqr(&v_eqr_south, p2, z_eps_max, n_iter_max) {
        res.push(lonlat);
      }
    } else {
      let v_eqr_north = Coo3D::from(intersect_with_transition_lat_npc(p1, p2).unwrap());
      if let Some(lonlat) =
        arc_special_point_in_eqr(&v_eqr_south, &v_eqr_north, z_eps_max, n_iter_max)
      {
        res.push(lonlat);
      }
      if let Some(lonlat) = arc_special_point_in_pc(&v_eqr_north, p2, z_eps_max, n_iter_max) {
        res.push(lonlat);
      }
    }
    res.into_boxed_slice()
  } else {
    // both EQR and NPC
    let mut res: Vec<LonLat> = Vec::with_capacity(2);
    let v_eqr_north = Coo3D::from(intersect_with_transition_lat_npc(p1, p2).unwrap());
    if let Some(lonlat) = arc_special_point_in_eqr(p1, &v_eqr_north, z_eps_max, n_iter_max) {
      res.push(lonlat);
    }
    if let Some(lonlat) = arc_special_point_in_pc(&v_eqr_north, p2, z_eps_max, n_iter_max) {
      res.push(lonlat);
    }
    res.into_boxed_slice()
  }
}

///////////////////////
// EQUATORIAL REGION //
///////////////////////

/// Returns the coordinate `z` (= sin(lat)) of the point such that the tangent line to the cone
/// on the projection plane has a slope equals to `+-1`, i.e. `d(DeltaX)/dY = +-1`.
/// The longitude can then be computed by solving the equation of the cone.
/// The sign of the slope is given by the parameter `positive_slope`.
///
/// This method is valid on the Equatorial Zone only (cylindrical equal area projection).
/// It is ok for radius > 20 mas:
/// - if radius < 20 mas = 1e-7 rad, then 1 - 2*sin^2(r/2) = 1 - 1e-14 which is close from the
///   precision of a double floating point.
///
/// The function rely on the Newton-Raphson method to solve `(DeltaX)/dY -+ 1 = 0`.
///
// TEST
// ```math
// \boxed{
//   \left\{
//     \begin{array}{lcl}
//       X & = & \alpha \times \frac{4}{\pi} \\
//       Y & = & \sin(\delta) \times \frac{3}{2}
//     \end{array}
//   \right.
// }
// \Rightarrow
// \left\{
//   \begin{array}{lcl}
//     \alpha \in [0, 2\pi] & \leadsto &  X \in [0, 8] \\
//     \sin\delta \in [-\frac{2}{3}, \frac{2}{3}] & \leadsto & Y \in [-1, 1]
//   \end{array}
// \right.
// ```
///
/// # Inputs
/// - `z` the Newton-Raphson method starting point. Reasonable choices are
///   - `sin(lat + 0.8 * radius)` when looking for the north special points (returned `z` > `z_cone_centre`)   
///   - `sin(lat - 0.8 * radius)` when looking for the south special points (returned `z` < `z_cone_centre`)
/// - `z0` the sine of the latitude of the center of the cone (= z_0)
/// - `eucl_cone_radius` the euclidean radius of the cone, i.e. $`2\sin(\frac{\theta}{2})`$
/// - `north_point`:
///   - `true`  if you look for `z` > `z_cone_centre` (i.e. equations Y = -X + b, north part of the cone)  
///   - `false` if you look for `z` < `z_cone_centre` (i.e. equations Y =  X + b, south part of the cone)
/// - `z_eps_max` the target precision on z (used to stop the Newton-Raphson method),
///   a reasonable choice may be `eucl_cone_radius / 1000`.
///   Internally, the value can't be higher than `eucl_cone_radius / 50`
///   nor lower than `1e-15`.
/// - `n_iter_max` upper limit on the number of iteration to be used in the Newton-Raphson method
///   a reasonable choice may be 20.
/// # Output
/// - `Some(z)` the sine of the latitude of the point such that the tangent line to the cone
///   on the projection plane has a slope equals to `+-1`.
/// - `None` if the computed latitude is out of the Equatorial Region
#[allow(dead_code)]
pub fn cone_special_point_in_eqr(
  mut z: f64,
  z0: f64,
  eucl_cone_radius: f64,
  north_point: bool,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Option<f64> {
  // Compute constants
  let cte = if north_point {
    // negative slope
    -ONE_OVER_TRANSITION_Z * FRAC_PI_4
  } else {
    // positive slope
    ONE_OVER_TRANSITION_Z * FRAC_PI_4
  };
  let w0 = 1.0 - z0.pow2();
  let r = 1.0 - eucl_cone_radius.pow2().half();
  //  Newton-Raphson method
  let z_eps_max = z_eps_max.min(0.2e-1 * eucl_cone_radius).max(1.0e-15);
  let mut z_eps = 1.0_f64;
  let mut n_iter = 0_u8;
  while n_iter < n_iter_max && z_eps.abs() > z_eps_max {
    z_eps = f_over_df_eqr(z, z0, w0, cte, r);
    z -= z_eps;
    n_iter += 1;
  }
  debug_assert!(z.is_finite());
  if z.abs() < TRANSITION_Z {
    debug_assert!((north_point && z >= z0) || (!north_point && z <= z0));
    Some(z)
  } else {
    None
  }
}

/// Check if the great-circle arc (defined by the smallest distance between the two
/// given points) contains a 'special point', i.e. a point such such that a tangent line
/// arc on the projection plane has a slope equals to `+-1`, i.e. `d(DeltaX)/dY = +-1`.
///
/// # Inputs
/// - `p1` first point of the great-circle arc
/// - `p2` second point of the great_-circle arc
/// - `p1_cross_p2` the (not necessarily normalized) cross-product `p1 x p2` (we pass it as a
///   parameter since it is already computed when working on polygons, the cross product can be
///   indifferently `p1 x p2` or `p2 x p1`)
/// - `z_eps_max` the target precision on z (used to stop the Newton-Raphson method),
///   a reasonable choice may be `|p2.z - p1.z|/ 1000`.
///   Internally, the value can't be higher than `|z2 - z1| / 50`
///   nor lower than `1e-15`.
/// - `n_iter_max` upper limit on the number of iteration to be used in the Newton-Raphson method
///   a reasonable choice may be 20.
/// # Output
/// - `Some(z)` the sine of the latitude of the point such that the tangent line to great-circle
///   arc (if it exists) has a slope equals to `+-1`.
/// - `None` if the arc do not contains a tangent line of slope '+-1', or if the result
///   if not in the equatorial region.
pub fn arc_special_point_in_eqr(
  p1: &Coo3D,
  p2: &Coo3D,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Option<LonLat> {
  let cone_center = cross_product(p1, p2).normalized();
  let z0 = cone_center.z();
  let z1 = p1.z();
  let z2 = p2.z();
  let north_point = z0 < 0.0;
  // Compute constants
  let cte = if north_point {
    -ONE_OVER_TRANSITION_Z * FRAC_PI_4
  } else {
    ONE_OVER_TRANSITION_Z * FRAC_PI_4
  };
  // Remark: r = 1 - 2 sin^2(pi/2 / 2) = 1 - 2 * (sqrt(2)/2)^2 = 0
  let w0 = 1.0 - z0.pow2();
  // Test if we start the method or not.
  let d1 = f_eqr(z1, z0, w0, cte, 0.0);
  let d2 = f_eqr(z2, z0, w0, cte, 0.0);
  if have_same_sign(d1, d2) {
    return None;
  }
  // Newton-Raphson method
  let z_eps_max = z_eps_max.min(0.2e-1 * (z2 - z1).abs()).max(1.0e-15);
  let mut z = (z1 + z2).half(); // mean of z1 and z2
  let mut z_eps = 1.0_f64;
  let mut n_iter = 0_u8;
  while n_iter < n_iter_max && z_eps.abs() > z_eps_max {
    z_eps = f_over_df_eqr(z, z0, w0, cte, 0.0);
    z -= z_eps;
    n_iter += 1;
  }
  // z must be in the [z1, z2] range except if Newton-Raphson method fails (divergence or to slow convergence)
  // TODO: if the condition ((z1 <= z2 && z1 <= z && z <= z2) || (z2 < z1 &&  z2 <= z && z <= z1))
  // TODO: is not met, then try a binary/dichotomic approach (slower, but more robust)
  if z.abs() < TRANSITION_Z && ((z1 <= z2 && z1 <= z && z <= z2) || (z2 < z1 && z2 <= z && z <= z1))
  {
    intersect_parallel(p1, p2, z).map(|v| v.lonlat())
  } else {
    None
  }
}

/// Computes dX / dY in the equatorial region
#[inline]
#[allow(clippy::many_single_char_names)]
fn f_eqr(z: f64, z0: f64, w0: f64, cte: f64, r: f64) -> f64 {
  let w = 1.0 - z.pow2(); // in equatortial region, -2/3 < z < 2/3
  let q = z / w; // so q is always defined
  let n = r - z * z0;
  (z0 - q * n) / (w0 * w - n.pow2()).sqrt() - cte
}

/// Computes the ratio (dX / dY) / (d^2X / dY^2) in the equatorial region
#[inline]
#[allow(clippy::many_single_char_names)]
fn f_over_df_eqr(z: f64, z0: f64, w0: f64, cte: f64, r: f64) -> f64 {
  let w = 1.0 - z.pow2();
  let q = z / w;
  let n = r - z * z0;
  let sqrt_d2_minus_n2 = (w0 * w - n.pow2()).sqrt();
  let qn = q * n;
  let dalphadz = (z0 - qn) / sqrt_d2_minus_n2;
  let f = dalphadz - cte;
  let df = (q * (z0.twice() - 3.0 * qn) - n * (1.0 / w + dalphadz.pow2())) / sqrt_d2_minus_n2;
  f / df
}

////////////////
// POLAR CAPS //
////////////////

/// Returns the coordinate `z` (=sin(lat)) of the point such that the tangent line to the cone
/// on the projection plane has a slope equals to `+-1`, i.e. `d(DeltaX)/dY = +-1`.
///
/// This method is valid on the polar caps only (Collignon projection).
/// It is ok for radius > 20 mas:
/// - if radius < 20 mas = 1e-7 rad, then 1 - 2*sin^2(r/2) = 1 - 1e-14 which is close from the
///   precision of a double floating point.
///
/// The function rely on the Newton-Raphson method to solve `(DeltaX)/dY -+ 1 = 0`.
///
/// # Inputs
/// - `z` the Newton-Raphson method starting point. Reasonable choices are
///   - `sin(lat - 0.9 * radius)` for positive slopes (returned `z` < `z_cone_centre`)
///   - `sin(lat + 0.9 * radius)` for negative slopes (returned `z` > `z_cone_centre`)
/// - `z0` the sine of the latitude of the center of the cone (= z_0)
/// - `eucl_cone_radius` the euclidean radius of the cone, i.e. `2*sin(ang_radius / 2)`
/// - `positive_slope`:
///   - `true`  if you look for `z` < `z_cone_centre` (i.e. equations Y =  X + b, south part of the cone)
///   - `false` if you look for `z` > `z_cone_centre` (i.e. equations Y = -X + b, north part of the cone)
/// - `z_eps_max` the target precision on z (used to stop the Newton-Raphson method),
///   a reasonable choice may be `cone_radius_rad / 1000`.
///   Internally, the value can't be higher than `eucl_cone_radius / 50`
///   nor lower than `1e-15`.
/// - `n_iter_max` upper limit on the number of iteration to be used in the Newton-Raphson method
///
/// # Output
/// - `Some(z)` the sine of the latitude of the point such that the tangent line to the cone
///   on the projection plane has a slope equals to `+-1`.
/// - `None` if the computed latitude is out of the North Polar Cap
///
///
/// # Warning
/// - if the longitude computed from the returned `z` is
///   - `> pi/2` in case of `east_value = true`
///   - `< 0` in case of `east_value = false`
/// - then the solution must be rejected and re-computed considering:
///   - case `north_value = true`
///     - `cone_center_lon_mod_half_pi - pi/2` and still North-East in case of `east_value = true`
///     - `cone_center_lon_mod_half_pi + pi/2` and still North-West in case of `east_value = false`
///   - case `north_value = false`
///     - `cone_center_lon_mod_half_pi - pi/2` and still South-East in case of `east_value = true`
///     - `cone_center_lon_mod_half_pi + pi/2` and still South-West in case of `east_value = false`
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn cone_special_point_in_pc(
  mut z: f64,
  cone_center_lon_mod_half_pi: f64,
  mut z0: f64,
  eucl_cone_radius: f64,
  east_value: bool,
  mut north_value: bool,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Option<f64> {
  let spc = z < 0.0; // south polar cap
  if spc {
    z = -z;
    z0 = -z0;
    north_value = !north_value;
  }
  // Compute constants
  let cte = if north_value { -FRAC_PI_8 } else { FRAC_PI_8 };
  let w0 = 1.0 - z0.pow2();
  let r = 1.0 - eucl_cone_radius.pow2().half();
  let direction = if east_value { 1.0 } else { -1.0 };
  // Newton-Raphson method
  let z_eps_max = z_eps_max.min(0.2e-1 * eucl_cone_radius).max(1.0e-15);
  let mut n_iter = 0_u8;
  let mut z_eps = 1.0_f64;
  while n_iter < n_iter_max && z_eps.abs() > z_eps_max {
    z_eps = f_over_df_npc(z, cone_center_lon_mod_half_pi, z0, w0, cte, direction, r);
    z -= z_eps;
    n_iter += 1;
  }
  // Cas cone contient le pole ou deborde sur les autres quarter A TRAITER ICI ?!
  if z.is_finite() && z > TRANSITION_Z {
    Some(if spc { -z } else { z })
  } else {
    None
  }
}

// check_multi_base_cells: must be set to true at first call!!
pub fn arc_special_point_in_pc<'a>(
  mut p1: &'a Coo3D,
  mut p2: &'a Coo3D,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Option<LonLat> {
  // Ensure p1.lon() < p2.lon()
  if p1.lon() > p2.lon() {
    std::mem::swap(&mut p1, &mut p2);
  }
  // Check if great-circle arc overlap several base cells
  debug_assert!(p1.lon() % FRAC_PI_2 >= 0.0);
  debug_assert!(p2.lon() % FRAC_PI_2 >= 0.0);
  let lon1_div_half_pi = p1.lon().div_eucl(FRAC_PI_2) as u8;
  let lon2_div_half_pi = p2.lon().div_eucl(FRAC_PI_2) as u8;
  debug_assert!(lon1_div_half_pi < 4);
  debug_assert!(lon2_div_half_pi < 4);
  debug_assert!(lon1_div_half_pi <= lon2_div_half_pi);
  // TODO: clean the code to make it more compact: i.e., method larger plane / smaller plane
  if lon1_div_half_pi != lon2_div_half_pi {
    let mut res_z1 = None;
    let mut res_z2 = None;
    // Info to understand the code
    // - the normal to the plane ra = q * pi/2, with q = 0 or 2 is n = (0, 1, 0)
    // - the normal to the plane ra = q * pi/2, with q = 1 or 3 is n = (1, 0, 0)
    if p2.lon() - p1.lon() > PI {
      // Cross lon = 0
      let p2p1_n = cross_product(p2, p1).normalized();
      if p2.lon() % FRAC_PI_2 > 0.0 {
        // First quarter, [p2.lon, ((p2.lon % PI/2) + 1) * PI/2]
        debug_assert!(lon2_div_half_pi > 0);
        let n2_y = lon2_div_half_pi & 1;
        let n2 = Coo3D::from_vec3((n2_y ^ 1) as f64, n2_y as f64, 0.0);
        let intersect2 = Coo3D::from(intersect_point_pc(p1, p2, &p2p1_n, &n2));
        debug_assert!(
          p2.lon() < intersect2.lon() || intersect2.lon() < p1.lon(),
          "p1.lon: {}, p2.lon: {}, intersect.lon: {}",
          p1.lon(),
          p2.lon(),
          intersect2.lon()
        );
        if p2.lon() < intersect2.lon() {
          res_z2 = arc_special_point_in_pc_same_quarter(p2, &intersect2, z_eps_max, n_iter_max);
        }
      }
      // Second quarter [(p1.lon % PI/2) * PI/2, p1.lon]
      debug_assert!(lon1_div_half_pi < 3);
      let n1_x = lon1_div_half_pi & 1;
      let n1 = Coo3D::from_vec3(n1_x as f64, (n1_x ^ 1) as f64, 0.0);
      let intersect1 = Coo3D::from(intersect_point_pc(p1, p2, &p2p1_n, &n1));
      debug_assert!(
        intersect1.lon() < p1.lon(),
        "p1: ({}, {}); p2: ({}, {}); intersect: ({}, {}); n; ({}, {})",
        p1.lon().to_degrees(),
        p1.lat().to_degrees(),
        p2.lon().to_degrees(),
        p2.lat().to_degrees(),
        intersect1.lon().to_degrees(),
        intersect1.lat().to_degrees(),
        n1.lon().to_degrees(),
        n1.lat().to_degrees()
      );
      res_z1 = arc_special_point_in_pc_same_quarter(&intersect1, p1, z_eps_max, n_iter_max);
    } else {
      let p1p2_n = cross_product(p1, p2).normalized();
      if p1.lon() % FRAC_PI_2 > 0.0 {
        // First quarter, [p1.lon, ((p1.lon % PI/2) + 1) * PI/2]
        debug_assert!(lon1_div_half_pi < 3);
        let n1_y = lon1_div_half_pi & 1;
        let n1 = Coo3D::from_vec3((n1_y ^ 1) as f64, n1_y as f64, 0.0);
        let intersect1 = Coo3D::from(intersect_point_pc(p1, p2, &p1p2_n, &n1));
        debug_assert!(p1.lon() < intersect1.lon());
        res_z1 = arc_special_point_in_pc_same_quarter(p1, &intersect1, z_eps_max, n_iter_max);
      }
      // Last quarter [(p2.lon % PI/2) * PI/2, p2.lon]
      debug_assert!(lon2_div_half_pi > 0);
      let n2_x = lon2_div_half_pi & 1;
      let n2 = Coo3D::from_vec3(n2_x as f64, (n2_x ^ 1) as f64, 0.0);
      let intersect2 = Coo3D::from(intersect_point_pc(p1, p2, &p1p2_n, &n2));
      debug_assert!(
        intersect2.lon() <= p2.lon(),
        "Failed: {} < {}",
        intersect2.lon(),
        p2.lon()
      );
      if intersect2.lon() < p2.lon() {
        res_z2 = arc_special_point_in_pc_same_quarter(&intersect2, p2, z_eps_max, n_iter_max);
      } else {
        res_z2 = None
      }
    }
    if res_z1.is_some() {
      res_z1
    } else {
      res_z2
    }
  } else {
    // Same quarter
    arc_special_point_in_pc_same_quarter(p1, p2, z_eps_max, n_iter_max)
  }
}

// Here we assume that the great-circle arc is in a same quarter
// (i.e. (p1.lon % pi/2) == (p2.lon % pi/2) (except if one of the two point is on a border n * pi/2)
fn arc_special_point_in_pc_same_quarter(
  p1: &Coo3D,
  p2: &Coo3D,
  z_eps_max: f64,
  n_iter_max: u8,
) -> Option<LonLat> {
  debug_assert!(
    p1.lon() < p2.lon(),
    "p1: ({}, {}); p2: ({}, {})",
    p1.lon().to_degrees(),
    p1.lat().to_degrees(),
    p2.lon().to_degrees(),
    p2.lat().to_degrees()
  );
  let mut p2_mod_half_pi = p2.lon() % FRAC_PI_2;
  if p2_mod_half_pi == 0.0 {
    p2_mod_half_pi = FRAC_PI_2;
  }
  let v1 = Coo3D::from_sph_coo(p1.lon() % FRAC_PI_2, p1.lat());
  let v2 = Coo3D::from_sph_coo(p2_mod_half_pi, p2.lat());
  let mut cone_center = cross_product(&v1, &v2).normalized();
  let lonlat = cone_center.lonlat();
  if lonlat.lon() > PI {
    cone_center = cone_center.opposite();
  }
  let cone_center_lon = cone_center.lonlat().lon();
  //debug_assert!(0 <= cone_center_lon && cone_center_lon <= );
  let mut z0 = cone_center.z();
  let mut z1 = v1.z();
  let mut z2 = v2.z();
  let mut north_value = z0 < 0.0;
  // ( here we could have but do not use the fact that we previously ensure that p1.lon() < p2.lon() )
  let east_value = ((v1.lat() > v2.lat()) ^ (v1.lon() > v2.lon())) ^ !north_value;
  // Deal with NPC / SPC
  let mut z = (z1 + z2).half(); // (0.1 * z1 + 0.9 * z2); //(z1 + z2).half(); // mean of z1 and z2
  let spc = z < 0.0; // south polar cap
  if spc {
    z = -z;
    z1 = -z1;
    z2 = -z2;
    z0 = -z0;
    north_value = !north_value;
  }
  // Compute constants
  //  - remark: r = 1 - 2 sin^2(pi/2 / 2) = 1 - 2 * (sqrt(2)/2)^2 = 0
  let cte = if north_value { -FRAC_PI_8 } else { FRAC_PI_8 };
  let w0 = 1.0 - z0.pow2();
  let direction = if east_value { 1.0 } else { -1.0 };
  // Test if we start the method or not
  let d1 = f_npc(z1, cone_center_lon, z0, w0, cte, direction, 0.0);
  let d2 = f_npc(z2, cone_center_lon, z0, w0, cte, direction, 0.0);
  if have_same_sign(d1, d2) {
    return None;
  }
  // Choose an initial value
  let dz = f_over_df_npc(z, cone_center_lon, z0, w0, cte, direction, 0.0);
  z -= dz;
  if !((z1 < z && z < z2) || (z2 < z && z < z1)) {
    z = z2 - f_over_df_npc(z2, cone_center_lon, z0, w0, cte, direction, 0.0);
    if !((z1 < z && z < z2) || (z2 < z && z < z1)) {
      z = z1 - f_over_df_npc(z1, cone_center_lon, z0, w0, cte, direction, 0.0);
    }
  }
  // Newton-Raphson method
  let z_eps_max = z_eps_max.min(0.2e-1 * (z2 - z1).abs()).max(1.0e-15);
  let mut n_iter = 0_u8;
  let mut z_eps = 1.0_f64;
  while n_iter < n_iter_max && z_eps.abs() > z_eps_max {
    z_eps = f_over_df_npc(z, cone_center_lon, z0, w0, cte, direction, 0.0);
    z -= z_eps;
    n_iter += 1;
  }
  // Return result if seems correct
  if z.is_finite() && z > TRANSITION_Z && ((z1 < z && z < z2) || (z2 < z && z < z1)) {
    if spc {
      let v = intersect_parallel(p1, p2, -z).unwrap();
      Some(v.lonlat())
    } else {
      let v = intersect_parallel(p1, p2, z).unwrap();
      Some(v.lonlat())
    }
  } else {
    None
  }
}

/// Returns the intersection point between the given arc (of given normal vector)
/// and the plane of given normal vector
/// WARNING: only valid in polar caps since we use 'z' to determine if we have to take the
/// result of (p1 x p2) x n or its complements (here we have the guarantee that
/// sign(p1.z) = sign(p2.z) must be = sign(res.z)
#[inline]
fn intersect_point_pc(p1: &Coo3D, p2: &Coo3D, p1_x_p2: &UnitVect3, n: &Coo3D) -> UnitVect3 {
  debug_assert!(p1.z().abs() >= TRANSITION_Z && p2.z().abs() >= TRANSITION_Z);
  debug_assert_eq!(p1.z() > 0.0, p2.z() > 0.0);
  let intersect = cross_product(p1_x_p2, n).normalized();
  if !have_same_sign(intersect.z(), p1.z()) {
    intersect.opposite()
  } else {
    intersect
  }
}

/// Computes dX / dY in the north polar cap
#[inline]
#[allow(clippy::many_single_char_names)]
fn f_npc(
  z: f64,
  cone_center_lon_mod_half_pi: f64,
  z0: f64,
  w0: f64,
  cte: f64,
  direction: f64,
  r: f64,
) -> f64 {
  let w = 1.0 - z;
  let w2 = 1.0 - z.pow2();
  let q = z / w2;
  let n = r - z * z0;
  let d2 = w0 * w2;
  let sqrt_d2_minus_n2 = (d2 - n.pow2()).sqrt();
  let qn = q * n;
  let arccos = (n / d2.sqrt()).acos();
  let dalphadz = (z0 - qn) / sqrt_d2_minus_n2;
  direction * w * dalphadz - 0.5 * (direction * arccos + cone_center_lon_mod_half_pi - FRAC_PI_4)
    + cte
}

/// Computes the ratio (dX / dY) / (d^2X / dY^2) in the north polar cap
#[inline]
#[allow(clippy::many_single_char_names)]
fn f_over_df_npc(
  z: f64,
  cone_center_lon_mod_half_pi: f64,
  z0: f64,
  w0: f64,
  cte: f64,
  direction: f64,
  r: f64,
) -> f64 {
  let w = 1.0 - z;
  let w2 = 1.0 - z.pow2();
  let q = z / w2;
  let n = r - z * z0;
  let d2 = w0 * w2;
  let sqrt_d2_minus_n2 = (d2 - n.pow2()).sqrt();
  let qn = q * n;
  let arccos = (n / d2.sqrt()).acos();
  let dalphadz = (z0 - qn) / sqrt_d2_minus_n2;
  let f = direction * w * dalphadz
    - 0.5 * (direction * arccos + cone_center_lon_mod_half_pi - FRAC_PI_4)
    + cte;
  let df = -ONE_OVER_TRANSITION_Z * direction * dalphadz
    + (direction * w / sqrt_d2_minus_n2)
      * (q * (z0.twice() - 3.0 * qn) - n * (1.0 / w2 + dalphadz.pow2()));
  f / df
}

/// Returns the intersection of the given great-circle arc (defined by the smallest distance
/// between the two given points) and the small circle of given z (equation $`z=cte`$).
/// Let's use the following notations:
/// - Coordinates of $`\vec{a}\times\vec{b} = (x_0, y_0, z_0)`$
/// - Coordinates of the points we are looking for $`\vec{i} = (x, y, z=cte)`$
///
/// We look for `x` and `y` solving
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     x^2 + y^2 + z^2 & = & 1 \\
///     xx_0 + yy_0 + zz_0 & = & 0 \\
///     (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 & = & 2 \mathrm(unused)
///   \end{array}
/// \right.
/// ```
/// It leads to
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     y & = & - \left(x\frac{x_0}{y_0} - \frac{zz_0}{y_0}\right) \\
///     0 & = & x^2(1+\frac{x_0^2}{y_0^2}) + x\frac{2x_0zz_0}{y_0^2} + (\frac{zz_0}{y_0})^2 + z^2 - 1
///   \end{array}
/// \right.
/// ```
/// We solve the quadratic equation
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     ax^2 + bx + c & = & 0 \\
///     \Delta & = & b^2 - 4ac \\
///     x & = & \frac{-b\pm \sqrt{\Delta}}{2a}
///   \end{array}
/// \right.
/// ```
/// If $`y_0 = 0`$, we directly derive the result:
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     x & = & -\frac{zz_0}{x_0} \\
///     y & = & \pm\sqrt{1 - x^2 - z^2} = \pm\sqrt{1 - z^2(1 + \frac{z_0^2}{x_0^2})}
///   \end{array}
/// \right.
/// ```
/// In both cases, two solutions are available.
/// We select the pair $`(x, y)`$ such that both scalar products with the two great-circle arc
/// points are higher than the two points scalar product.
pub fn intersect_parallel<T1, T2>(p1: &T1, p2: &T2, z: f64) -> Option<UnitVect3>
where
  T1: UnitVec3,
  T2: UnitVec3,
{
  debug_assert!((-1.0..=1.0).contains(&z));
  if (p1.z() < z && z < p2.z()) || (p2.z() < z && z < p1.z()) {
    let p1_dot_p2 = dot_product(p1, p2);
    let p1_x_p2 = cross_product(p1, p2);
    let p1_x_p2_norm = p1_x_p2.norm();
    let ang_dist = p1_x_p2_norm.atan2(p1_dot_p2);
    let p1_x_p2 = p1_x_p2.normalized();
    let x0 = p1_x_p2.x();
    let y0 = p1_x_p2.y();
    let z0 = p1_x_p2.z();
    if y0.abs() <= 1e-14 {
      let x = -(z * z0) / x0;
      let y1 = (1.0 - (x.pow2() + z.pow2())).sqrt();
      let y2 = -y1;
      if p1.x() * x + p1.y() * y1 + p1.z() * z >= p1_dot_p2
        && p2.x() * x + p2.y() * y1 + p2.z() * z >= p1_dot_p2
      {
        Some(UnitVect3::new_unsafe(x, y1, z))
      } else if p1.x() * x + p1.y() * y2 + p1.z() * z >= p1_dot_p2
        && p2.x() * x + p2.y() * y2 + p2.z() * z >= p1_dot_p2
      {
        Some(UnitVect3::new_unsafe(x, y2, z))
      } else {
        unreachable!();
      }
    } else if ang_dist < (1_f64 / 3600.0).to_degrees() {
      // dist < 1 arcsec: use local flat approximation and compute position on the segment
      //   x = (x2 - x1)t + x1
      //   y = (y2 - y1)t + y1
      //   z = (z2 - z1)t + z1
      // => t = (z - z1) / (z2 - z1)
      let t = (z - p1.z()) / (p2.z() - p1.z());
      let x = (p2.x() - p1.x()) * t + p1.x();
      let y = (p2.y() - p1.y()) * t + p1.y();
      Some(UnitVect3::new(x, y, z))
    } else {
      let x0_y0 = x0 / y0;
      let zz0_y0 = z * z0 / y0;
      let a = 1.0 + x0_y0.pow2();
      let b = 2.0 * x0_y0 * zz0_y0;
      let c = zz0_y0.pow2() + z.pow2() - 1.0;
      let delta = b.pow2() - 4.0 * a * c;
      let sqrt_delta = delta.sqrt();
      let x1 = (-b + sqrt_delta) / a.twice();
      let y1 = -x1 * x0_y0 - zz0_y0;
      let x2 = (-b - sqrt_delta) / a.twice();
      let y2 = -x2 * x0_y0 - zz0_y0;
      if p1.x() * x1 + p1.y() * y1 + p1.z() * z >= p1_dot_p2
        && p2.x() * x1 + p2.y() * y1 + p2.z() * z >= p1_dot_p2
      {
        Some(UnitVect3::new_unsafe(x1, y1, z))
      } else if p1.x() * x2 + p1.y() * y2 + p1.z() * z >= p1_dot_p2
        && p2.x() * x2 + p2.y() * y2 + p2.z() * z >= p1_dot_p2
      {
        Some(UnitVect3::new_unsafe(x2, y2, z))
      } else {
        unreachable!();
      }
    }
  } else {
    None
  }
}

#[inline]
fn have_same_sign(d1: f64, d2: f64) -> bool {
  d1 == 0.0 || d2 == 0.0 || ((d1 > 0.0) == (d2 > 0.0))
}

/// Returns the intersection of the given great-circle arc (defined by the smallest distance
/// between the two given points) and the small circle of equation $`z=2/3`$.
/// (Internally, we simply call [intersect_parallel](fn.intersect_parallel.html) with
/// z = 2/3).
#[inline]
fn intersect_with_transition_lat_npc<T1, T2>(p1: &T1, p2: &T2) -> Option<UnitVect3>
where
  T1: UnitVec3,
  T2: UnitVec3,
{
  intersect_parallel(p1, p2, TRANSITION_Z)
}

/// Returns the intersection of the given great-circle arc (defined by the smallest distance
/// between the two given points) and the small circle of equation $`z=-2/3`$.
/// (Internally, we simply call [intersect_parallel](fn.intersect_parallel.html) with
/// z = -2/3).
#[inline]
fn intersect_with_transition_lat_spc<T1, T2>(p1: &T1, p2: &T2) -> Option<UnitVect3>
where
  T1: UnitVec3,
  T2: UnitVec3,
{
  intersect_parallel(p1, p2, -TRANSITION_Z)
}
