
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
  pub fn contains(&self, lon: f64, lat: f64)-> bool {
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
  pub fn overlap_cone(&self, lon: f64, lat: f64, radius: f64) -> bool {
    match self.center.proj(lon, lat) {
      Some((x, y)) => {
        let eucl_a = (self.a + radius).sin();
        let eucl_b = (self.b + radius).sin();
        Ellipse::from_oriented(eucl_a, eucl_b, self.theta_sin_cos).contains(x, y)
      },
      None => false,
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
