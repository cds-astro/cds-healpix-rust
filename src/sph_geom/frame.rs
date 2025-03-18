use crate::TWICE_PI;

/// Components of the 3x3 rotation matrix transforming a vector into the
/// reference frame to a vector into the local frame (i.e. the frame in
/// which the position of the projection center is (1, 0, 0).
/// Remark:
/// * `r22 =  cos(lon)`
/// * `r21 = -sin(lon)`
/// * `r33 =  cos(lat)`
/// * `r13 =  sin(lat)`
#[derive(Debug)]
pub struct RefToLocalRotMatrix {
  r11: f64,
  r12: f64,
  r13: f64,
  r21: f64,
  r22: f64,
  r23: f64,
  r31: f64,
  r32: f64,
  r33: f64,
}

impl RefToLocalRotMatrix {
  /// Compute the reference to local rotation matrix from a projection center
  pub fn from_center(lon: f64, lat: f64) -> RefToLocalRotMatrix {
    let (sa, ca) = lon.sin_cos(); // ca, sa stands for cos(alpha), sin(alpha)
    let (sd, cd) = lat.sin_cos(); // cd, sd stands for cos(delta), sin(delta)
    RefToLocalRotMatrix {
      r11: ca * cd,
      r12: sa * cd,
      r13: sd,
      r21: -sa,
      r22: ca,
      r23: 0.0,
      r31: -ca * sd,
      r32: -sa * sd,
      r33: cd,
    }
  }

  /// Transform local to global (or reference) coordinates, by
  /// rotating the input (local) vector using the transpose of the rotation matrix.
  pub fn to_global_xyz(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    (
      self.r11 * x + self.r21 * y + self.r31 * z,
      self.r12 * x + self.r22 * y + self.r32 * z,
      self.r13 * x + self.r23 * y + self.r33 * z,
    )
  }

  /// Transform local to global (or reference) coordinates, by
  /// rotating the input (local) vector using the transpose of the rotation matrix.
  pub fn to_global_coo(&self, x: f64, y: f64, z: f64) -> (f64, f64) {
    let (x, y, z) = self.to_global_xyz(x, y, z);
    xyz_to_lonlat(x, y, z)
  }
}

fn xyz_to_lonlat(x: f64, y: f64, z: f64) -> (f64, f64) {
  // Length of the projection on the xy plane
  let r2 = x * x + y * y;
  // Latitude in [-pi/2, pi/2] (ok, since cos always positive here)
  let lat = z.atan2(r2.sqrt());
  // Compute the longitude in [-pi, pi]
  let r2 = y.atan2(x);
  // Conforms to convention: Longitude in [0, 2*PI]
  let lon = if r2 < 0.0 { TWICE_PI + r2 } else { r2 };
  (lon, lat)
}
