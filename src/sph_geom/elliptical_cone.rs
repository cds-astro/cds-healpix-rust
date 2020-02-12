use std::f64::NAN;
use std::f64::consts::FRAC_PI_2;

use super::proj::*;
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
    let theta_sin_cos = (FRAC_PI_2 - pa).sin_cos();
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
  
  /// Returns `true` if the given cone overlap the elliptical cone.
  /// # Inputs
  /// - `lon` longitude of the center of the cone, in radians
  /// - `lat` latitude of the center of the cone, in radians
  /// - `radius` cone radius, in radians
  /// 
  /// # Info
  /// To properly address the question, one should try to see if the following equations system has
  /// at least one solution:
  /// ```math
  /// \frac{x^2}{\tan^2a} + \frac{y^2}{\tan^2b} = z^2
  /// (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2  = 4\sin^2\frac{\theta}{2}
  /// x^2 + y^2 + z^2 = 0
  /// ```
  /// with `a` and `b` which are angles and
  /// ```math
  /// x = \cos\delta\cos\alpha
  /// y = \cos\delta\sin\alpha
  /// z = \sin\delta
  /// ```
  /// This system can be replaced by (the second equation being the angular distance formula):
  /// ```math
  /// \frac{\cos^2\alpha}{\tan^2a} + \frac{\sin^2\alpha}{\tan^2b} = \tan^2\delta
  /// \sin\delta\sin\delta_0 + \cos\delta\cos\delta_0\cos(\alpha-\alpha_0) = \cos\theta
  /// ```
  /// - The first equation is the elliptical cone equation
  /// - The second is the classical cone equation
  /// - $(x_ 0, y_ 0, z_ 0)$ are obtained by rotation of the original cone center such that the 
  ///   elliptical cone is center is $(x = 0, y = 0, z = 1)$, and the x-axis is along the semi-major
  ///   axis.
  /// 
  /// # Info 2
  /// Contrary to the above Info, another (the best?) approach is to compute the coordinates of the
  /// two ellipse foci F0 and F1. The sum of the distances f0 and f1 to both foci is constant with:
  /// ```math
  /// f0 + f1 = 2a
  /// ```
  /// 
  pub fn overlap_cone(&self, lon: f64, lat: f64, radius: f64) -> bool {
    assert!(radius > 0.0);
    let ((x, y), ang_dist, top_hemisphere) =  self.center.forced_proj_and_distance(lon, lat);
    if self.a + radius < ang_dist { // Quick rejection test
        return false;
    }
    // The projection of a cone on a plane is an ellipse
    let proj_d_min = (ang_dist - radius).sin();
    let proj_d_max = (ang_dist + radius).sin();
    let proj_a = radius.sin();                          // projected ellipse semi-major axis
    let proj_b = 0.5 * (proj_d_max - proj_d_min).abs(); // projected ellipse semi-minor axis
    let one_over_norm = 1.0 / (x.pow2() + y.pow2()).sqrt(); // distance from (0, 0) to the cone projected center (!= ellipse center)
    if !one_over_norm.is_finite() {
       // Special case: the center of the HEALPix cell is the center of the ellipse
       return radius <= self.b;
    }
    // phi = angle of the (x, y) position in the euclidean plane
    let cos_phi = x * one_over_norm;
    let sin_phi = y * one_over_norm;
    let sin_cos_angle = if sin_phi >= 0.0 {
      (-cos_phi, sin_phi)
    } else {
      (cos_phi, -sin_phi)
    };
    let proj_ell = Ellipse::from_oriented(proj_a, proj_b, sin_cos_angle);
    let proj_ell_dist = 0.5 * (proj_d_max + proj_d_min); // distance from (0, 0) to ellipse center
    let proj_ell_x = proj_ell_dist * cos_phi;
    let proj_ell_y = proj_ell_dist * sin_phi;
    if top_hemisphere {
      self.ellipse.overlap(proj_ell_x, proj_ell_y, &proj_ell)
    } else {
      // The projection of the center of the cone is in the back hemisphere
      // so the area to be considered looks like a crescent moon.
      // So far we apply the same algo, knowing that at deeper depth cells will have less changes
      // to spuriously overlap
      self.ellipse.overlap(proj_ell_x, proj_ell_y, &proj_ell)
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
      match self.center.proj(lon, lat) {
        Some((x, y)) => {
          let eucl_a = (self.a - radius).sin();
          let eucl_b = (self.b - radius).sin();
          // WARNING: not sure this is a 100% reliable for large distances!
          // A 100% reliable solution would be to compute the distance to both foci:
          //   (f0 + f1)/2 <= a - r
          Ellipse::from_oriented(eucl_a, eucl_b, self.theta_sin_cos).contains(x, y)
        },
        None => false,
        }
    }
  }

  pub fn path_along_edge(lon: f64, lat: f64, a: f64, b: f64, theta: f64, half_num_points: usize) -> Box<[(f64, f64)]> {
    let center = ProjSIN::new(lon, lat);
    let mut ar = Ellipse::path_along_edge(a.sin(), b.sin(), FRAC_PI_2 - theta, half_num_points);
    for coo in ar.iter_mut() {
      *coo = match center.unproj(coo.0, coo.1) {
        Some(lonlat) => lonlat,
        None => (NAN, NAN),
      };
    }
    ar
  }

}
