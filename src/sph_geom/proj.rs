
use super::super::HALF_PI;
use super::super::TWICE_PI;
use super::super::Customf64;


/// Represents a spherical projection, i.e. the projection of spherical coordinates
/// (on the unit sphere) on a two dimensional plane.
pub trait Proj: Default {
  /// Set the center of the projection
  /// # Inputs
  /// - `lon` longitude of the centre of the projection, in `[0, 2pi]` radians
  /// - `lat` latitude of the centre of the projection, in `[-pi/2, pi/2]` radians
  /// # info
  ///  If the longitude or the latitude is out of bounds, a time consuming normalization may be performed 
  fn set_center(&mut self, lon: f64, lat :f64);
  
  /// Computes the projected coordinates of the given position on the unit sphere
  /// # Inputs
  /// - `lon` longitude of the coordinate we want to project
  /// - `lat` latitude of the coordinate we want to project
  /// # Output
  /// - `(x, y)` the coordinates in the projection plane (if they exists)
  fn proj(&self, lon: f64, lat: f64) -> Option<(f64, f64)>;
  
  /// Computes the position in the unit sphere of the given coordinates on the projection plane
  /// # Inputs
  /// - `x` coordinates along the x-axis on the projection plane
  /// - `y` coordinates along the y-axis on the projection plane
  /// # Output
  /// - `(lon, lat)` the position on the unit sphere (in radians), if the input coordinates `(x, y)`
  ///    are valid
  fn unproj(&self, x: f64, y: f64) -> Option<(f64, f64)>;
}

/// Normalize the given longitude and latitude so that their values are in `[0, 2pi]` and  `[-pi/2, pi/2]`.
/// WARNING: the normalization (without a priori) is time consuming and should be avoided as much as possible
fn normalize_lonlat(lon: &mut f64, lat: &mut f64) {
  if *lon < 0.0 || TWICE_PI <= *lon || *lat < -HALF_PI || HALF_PI < *lat {
    let (cos_l, sin_l) = (*lon).sin_cos();
    let (cos_b, sin_b) = (*lat).sin_cos();
    let x = cos_b * cos_l;
    let y = cos_b * sin_l;
    let z = sin_b;
    *lon = y.atan2(x);
    if *lon < 0.0_f64 {
      *lon += TWICE_PI;
    }
    *lat = z.atan2((x.pow2() + y.pow2()).sqrt()); 
  }
}

/// Orthographic projection
#[derive(Debug)]
pub struct ProjSIN {
  center_lon: f64,
  center_lat: f64,
  // derived quantities
  cos_center_lat: f64,
  sin_center_lat: f64,
}

impl ProjSIN {
  
  /// # Inputs
  /// - `lon` longitude of the centre of the projection, in `[0, 2pi]` radians
  /// - `lat` latitude of the centre of the projection, in `[-pi/2, pi/2]` radians
  /// # info
  ///  If the longitude or the latitude is out of bounds, a time consuming normalization is performed 
  pub fn new(lon: f64, lat :f64) -> ProjSIN {
    let mut proj: ProjSIN = Default::default();
    proj.set_center(lon, lat);
    proj
  }
}

impl Default for ProjSIN {
  fn default() -> Self {
    ProjSIN {
      center_lon: 0.0,
      center_lat: 0.0,
      // derived quantities
      cos_center_lat: 1.0,
      sin_center_lat: 0.0,
    }
  }
}

impl Proj for ProjSIN {
  
  fn set_center(&mut self, lon: f64, lat :f64) {
    // I put assert tests here because the normalization is time consuming and I would like to
    // avoid it as much as possible
    debug_assert!(0.0 <= lon && lon < TWICE_PI);
    debug_assert!(-HALF_PI <= lat && lat < HALF_PI);
    self.center_lon = lon;
    self.center_lat = lat;
    normalize_lonlat(&mut self.center_lon, &mut self.center_lat);
    let (sin_lat, cos_lat) = self.center_lat.sin_cos();
    self.cos_center_lat = cos_lat;
    self.sin_center_lat = sin_lat;
  }
  
  fn proj(&self, lon: f64, lat: f64) -> Option<(f64, f64)> {
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_dlon, cos_dlon) = (lon - self.center_lon).sin_cos();
    if self.sin_center_lat * sin_lat - self.cos_center_lat * cos_lat * cos_dlon > 0.0 {
      Some((
        cos_lat * sin_dlon,
        self.cos_center_lat * sin_lat - self.sin_center_lat * cos_lat * cos_dlon
      ))
    } else {
      None
    }
  }
  
  fn unproj(&self, x: f64, y: f64) -> Option<(f64, f64)> {
    let rho2 = x.pow2() + y.pow2();
    if rho2 < 1.0 {
      let r = (1.0 - rho2).sqrt(); // = cos(arcsin(rho)) = rho / tan(arcsin(rho))
      let lat = (self.sin_center_lat * r + y * self.cos_center_lat).asin();
      let lon = (self.cos_center_lat * r - y * self.sin_center_lat);
      let lon = self.center_lon + x.atan2(lon);
      Some((
        if lon < 0.0 {
          lon + TWICE_PI
        } else if lon > TWICE_PI {
          lon - TWICE_PI
        } else {
          lon
        }, lat))
    } else {
      None
    }
  }
  
}