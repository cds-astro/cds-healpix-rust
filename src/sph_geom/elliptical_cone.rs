
use std::f64::NAN;

use super::coo3d::*;
use super::proj::*;
use super::super::HALF_PI;
use super::super::xy_geom::ellipse::*;
use super::super::Customf64;

#[derive(Debug)]
pub struct EllipticalCone {
  center: ProjSIN, // the center of the projection is the center of the cone 
  ellipse: Ellipse,
  // Ellipse parameters
  a: f64,  // semi-major axis, in radians
  b: f64,  // semi-minor axis, in radians
  theta_sin_cos: (f64, f64), // (pi/2 - position angle).sin_cos()
}

impl EllipticalCone {
  
  /// #Inputs
  /// - `lon` longitude of the ellipse center, in radians
  /// - `lat` latitude of the ellipse center, in radians
  /// - `a` semi-major axis of the ellipse, in radians
  /// - `b` semi-minor axis of the ellipse, in radians
  /// - `pa` position angle (east of north), in radians
  pub fn new(lon: f64, lat: f64, a: f64, b: f64, pa: f64) -> EllipticalCone {
    let center = ProjSIN::new(lon, lat);
    let theta_sin_cos = (HALF_PI - pa).sin_cos();
    let ellipse = Ellipse::from_oriented(a.sin(), b.sin(), theta_sin_cos);
    EllipticalCone {
      center,
      ellipse,
      a,
      b,
      theta_sin_cos,
    }
  }

  /// Returns `true` if the given point on the unit sphere is inside the elliptical cone.
  /// # Inputs
  /// - `lon` longitude of the point we want to test, in radians
  /// - `lat` latitude of the point we want to test, in radians
  pub fn contains(&self, lon: f64, lat: f64) -> bool {
    match self.center.proj(lon, lat) {
      Some((x, y)) => self.ellipse.contains(x, y),
      None => false,
    }
  }

  /*pub fn may_overlap_cone(&self, lon: f64, lat: f64, radius: f64) -> bool {// !! TEST DISTANCE < r+R
    if self.a + radius > HALF_PI { // in this case, possible out-of-bound projection for the center of the cone
      true
    } else {
      self.overlap_cone(lon, lat, radius)
    }
  }*/
  
  /// Returns `true` if the given cone overlap the elliptical cone.
  /// # Inputs
  /// - `lon` longitude of the center of the cone, in radians
  /// - `lat` latitude of the center of the cone, in radians
  /// - `radius` cone radius, in radians
  pub fn overlap_cone(&self, lon: f64, lat: f64, radius: f64) -> bool {
    if self.a + radius <= HALF_PI {
      match self.center.proj(lon, lat) {
        Some((x, y)) => {
          let eucl_a = (self.a + radius).sin();
          let eucl_b = (self.b + radius).sin();
          Ellipse::from_oriented(eucl_a, eucl_b, self.theta_sin_cos).contains(x, y)
        },
        None => false,
      }
    } else {
      let ((x, y), ang_dist, top_hemisphere) =  self.center.forced_proj_and_distance(lon, lat);
      if radius < 0.0 ||  self.a + radius < ang_dist {
        // radius < 0.0 => we are testing for inclusion (contains_cone)
        return false;
      }
      // The projection of a cone on a plane is an ellipse
      let proj_d_min = (ang_dist - radius).sin();
      let proj_d_max = (ang_dist + radius).sin();
      let proj_a = radius.sin();                     // projected ellipse semi-major axis
      let proj_b = 0.5 * (proj_d_max - proj_d_min).abs(); // projected ellipse semi-minor axis
      // ellipse theta angle = tan(t + pi/2) = cos(t) / (-sin(t)) = x / -y
      // let angle = (x / -y).atan(); //
      let one_over_norm = 1.0 / (x.pow2() + y.pow2()).sqrt(); // distance from (0.0) the cone projected center (!= ellipse center)
      if !one_over_norm.is_finite() {
        // Special case: the center of the HEALPix cell is the center of the ellipse
        return radius >= 0.0 || radius.abs() < self.b;
      }
      // phi = angle of the (x, y) position in the euclidean plane
      let cos_phi = x * one_over_norm;
      let sin_phi = y * one_over_norm;
      let sin_cos_angle = if sin_phi >= 0.0 {
        (-cos_phi, sin_phi)
      } else {
        (cos_phi, -sin_phi)
      };
      // eprintln!("Self a: {}; b: {}; theta: {}", self.a, self.b, self.theta_sin_cos.0.atan2(self.theta_sin_cos.1).to_degrees());
      // eprintln!("Othe a: {}; b: {}; theta: {}", proj_a, proj_b, sin_cos_angle.0.atan2(sin_cos_angle.1).to_degrees());
      let proj_ell = Ellipse::from_oriented(proj_a, proj_b, sin_cos_angle);
      let proj_ell_dist = 0.5 * (proj_d_max + proj_d_min); // distance from (0, 0) to ellipse center
      let proj_ell_x = proj_ell_dist * cos_phi;
      let proj_ell_y = proj_ell_dist * sin_phi;
      //eprintln!("proj_ell_x: {}; proj_ell_y: {}", proj_ell_x, proj_ell_y);
      if top_hemisphere {
        self.ellipse.overlap(proj_ell_x, proj_ell_y, &proj_ell)
      } else {
        // The projection of the center of the cone is in the back hemisphere
        // so the area to be considered look like a crescent moon.
        // I am not sure that the following is OK, let's try...
        // So far we apply the same algo, knowing that at deeper depth cells will have less changes
        // to spuriously overalp
        self.ellipse.overlap(proj_ell_x, proj_ell_y, &proj_ell)
      }
    }
  }

  /// Returns `true` if the given cone is fully inside the elliptical cone.
  /// # Inputs
  /// - `lon` longitude of the center of the cone, in radians
  /// - `lat` latitude of the center of the cone, in radians
  /// - `radius` cone radius, in radians
  pub fn contains_cone(&self, lon: f64, lat: f64, radius: f64) -> bool {
    if radius >= self.b {
      false
    } else {
      self.overlap_cone(lon, lat, -radius)
    }
  }


  pub fn path_along_edge(lon: f64, lat: f64, a: f64, b: f64, theta: f64, half_num_points: usize) -> Box<[(f64, f64)]> {
    let center = ProjSIN::new(lon, lat);
    let mut ar = Ellipse::path_along_edge(a.sin(), b.sin(), HALF_PI - theta, half_num_points);
    for coo in ar.iter_mut() {
      *coo = match center.unproj(coo.0, coo.1) {
        Some(lonlat) => lonlat,
        None => (NAN, NAN),
      };
    }
    ar
  }

}
