
use crate::{PI, TWICE_PI};
use crate::sph_geom::coo3d::{HALF_PI, vec3_of};
use crate::sph_geom::cone::{Cone, mec_3};

#[derive(Debug)]
pub struct Zone {
  /// Minimal longitude, in `[0, 2\pi[` radians
  lon_min: f64,
  /// Minimal latitude in `[-\pi/2, \pi/2]` radians
  lat_min: f64,
  /// Maximal longitude, in `[0, 2\pi[` radians
  lon_max: f64,
  /// Maximal latitude in `[-\pi/2, \pi/2]` radians
  lat_max: f64,
  /// Tells if the zone cross the rpimary meridian (in this case, lon_min > lon_max)
  cross_primary_meridian: bool,
}

impl Zone {
  /// # Remark
  /// * If `lon_min > lon_max` then we consider that the zone crosses the primary meridian.
  /// * The north pole is included only if `lon_min == 0 && lat_max == pi/2`
  /// # Panics
  /// * if `lon_min` or `lon_max` not in `[0, 2\pi[`
  /// * if `lat_min` or `lat_max` not in `[-\pi/2, \pi/2[`
  /// * `lat_min >= lat_max`
  pub fn new(lon_min: f64, lat_min: f64, lon_max: f64, lat_max: f64) -> Zone {
    assert!(0.0 <= lon_min && lon_min < TWICE_PI && 0.0 < lon_max && lon_max <= TWICE_PI);
    assert!(-HALF_PI <= lat_min && lat_min < HALF_PI && -HALF_PI < lat_max && lat_max <= HALF_PI);
    assert!(lat_min < lat_max);
    // Because of inequalities (< lat_max), we have to make an exception for the north pole
    let lat_max = if lat_max == HALF_PI { HALF_PI + 1e-15 } else { lat_max };
    Zone {
      lon_min,
      lat_min,
      lon_max,
      lat_max,
      cross_primary_meridian: lon_min > lon_max
    }
  }

  pub fn dlon(&self) -> f64 {
    if self.cross_primary_meridian {
      PI - (self.lon_min - self.lon_max)
    } else {
      self.lon_max - self.lon_min
    }
  }

  pub fn dlat(&self) -> f64 {
    self.lat_max - self.lat_min
  }

  pub fn crossed_vertically(&self, lon: f64, lat_up: f64, lat_down: f64) -> bool {
    lat_down < self.lat_min && self.lat_max <= lat_up &&
      if self.cross_primary_meridian {
        self.lon_min <= lon || lon < self.lon_max
      } else {
        self.lon_min <= lon && lon < self.lon_max
      }
  }

  pub fn crossed_horizontally(&self, lon_left: f64, lon_right: f64, lat: f64) -> bool {
    self.lat_min <= lat && lat < self.lat_max &&
      if self.cross_primary_meridian {
        if lon_right < lon_left {
          lon_left < self.lon_min && self.lon_max <= lon_right
        } else {
          false
        }
      } else {
        if lon_right < lon_left {
          lon_left < self.lon_min || self.lon_max <= lon_right
        } else {
          lon_left < self.lon_min && self.lon_max <= lon_right
        }
      }
  }

  /// Returns the center and the radius ((lon, lat), r) of the smallest
  /// cone containing the zone.
  /// If the cone is nonreflex (i.e. larger than an hemisphere), the result is `None`.
  pub fn smallest_enclosing_cone(&self) -> Option<Cone> {
    // Compute the minimum enclosing cone
    let v1 = vec3_of(self.lon_min, self.lat_min);
    let v2 = vec3_of(self.lon_max, self.lat_max);
    let v3 = vec3_of(self.lon_min, self.lat_max);
    let mec = mec_3(&v1, &v2, &v3);
    let lonlat = mec.center().lonlat();
    // Check that the center is in the zone!!
    if self.contains(lonlat.lon, lonlat.lat) {
      Some(mec)
    } else {
      None
    }
  }

  pub fn contains(&self, lon: f64, lat: f64) -> bool {
    self.lat_min <= lat && lat < self.lat_max &&
      if self.cross_primary_meridian {
        self.lon_min <= lon || lon < self.lon_max
      } else {
        self.lon_min <= lon && lon < self.lon_max
      }
  }

  pub fn vertices(&self) -> [(f64, f64); 4] {
    [
      (self.lon_min, self.lat_min),
      (self.lon_min, self.lat_max),
      (self.lon_max, self.lat_max),
      (self.lon_max, self.lat_min)
    ]
  }
}
