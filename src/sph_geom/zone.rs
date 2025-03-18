use crate::sph_geom::cone::{mec_3, Cone};
use crate::sph_geom::coo3d::{vec3_of, HALF_PI};
use crate::TWICE_PI;

#[derive(Debug)]
pub struct Zone {
  /// Minimal longitude, in `[0, 2\pi[` radians
  lon_min: f64,
  /// Minimal latitude in `[-\pi/2, \pi/2[` radians
  lat_min: f64,
  /// Maximal longitude, in `]0, 2\pi]` radians
  lon_max: f64,
  /// Maximal latitude in `]-\pi/2, \pi/2]` radians
  lat_max: f64,
  /// Tells if the zone cross the primary meridian (in this case, lon_min > lon_max)
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
    assert!((0.0..TWICE_PI).contains(&lon_min) && 0.0 < lon_max && lon_max <= TWICE_PI);
    assert!((-HALF_PI..HALF_PI).contains(&lat_min) && -HALF_PI < lat_max && lat_max <= HALF_PI);
    assert!(lat_min < lat_max);
    // Because of inequalities (< lat_max), we have to make an exception for the north pole
    let lat_max = if lat_max == HALF_PI {
      HALF_PI + 1e-15
    } else {
      lat_max
    };
    Zone {
      lon_min,
      lat_min,
      lon_max,
      lat_max,
      cross_primary_meridian: lon_min > lon_max,
    }
  }

  /*
  /// Returns the "center" of the zone, i.e. the mean longitude and the mean latitude.
  pub fn center(&self) -> (f64, f64) {
    let b = 0.5 * (self.lat_max + self.lat_min);
    let l = 0.5 * if self.cross_primary_meridian {
      let right = TWICE_PI - self.lon_min;
      let left = self.lon_max;
      if right > left {
        right - left
      } else {
        left - right
      }
    } else {
      (self.lon_max + self.lon_min)
    };
    (l, b)
  }
  */

  pub fn dlon(&self) -> f64 {
    if self.cross_primary_meridian {
      TWICE_PI - self.lon_min + self.lon_max
    } else {
      self.lon_max - self.lon_min
    }
  }

  pub fn dlat(&self) -> f64 {
    self.lat_max - self.lat_min
  }

  /// Returns `true` if the "vertical" great circle arc of given longitude `lon` and
  /// going from latitude `last_down` to latitude `lat_up` crosses totally this zone.
  /// # WARNING
  /// * the input must satisfy `last_down` < `lat_up`
  // We already consider case with a vertex inside
  pub fn crossed_vertically(&self, lon: f64, lat_down: f64, lat_up: f64) -> bool {
    let is_in_lat_range = lat_down < self.lat_min && self.lat_max <= lat_up;
    is_in_lat_range
      && if self.cross_primary_meridian {
        self.lon_min <= lon || lon < self.lon_max
      } else {
        self.lon_min <= lon && lon < self.lon_max
      }
  }

  /// Returns `true` if the small circle arc of given latitude `lat` going from
  /// longitude `lon_left` to longitude `lon_right` crosses totally this zone.
  /// # Remark
  /// * if the small circle does not cross the primary meridian, the input must satisfy `lon_left` < `lon_right`
  /// * else (if `lon_left` > `lon_right`), we consider that the small circle crosses the primary meridian.
  // We already consider case with a vertex inside
  #[allow(clippy::collapsible_else_if)]
  pub fn crossed_horizontally(&self, lon_left: f64, lon_right: f64, lat: f64) -> bool {
    let small_arc_circle_crosses_prim_meridian = lon_right < lon_left;
    let is_lat_compatible = (self.lat_min..self.lat_max).contains(&lat);
    is_lat_compatible
      && if self.cross_primary_meridian {
        if small_arc_circle_crosses_prim_meridian {
          lon_left < self.lon_min && self.lon_max <= lon_right
        } else {
          false
        }
      } else {
        if small_arc_circle_crosses_prim_meridian {
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

  /// Check whether or not the given lon range is fully inside the zone lon range.
  pub fn is_lon_range_compatible(&self, lon_w: f64, lon_e: f64) -> bool {
    let is_range_crossing_prim_meridian = lon_e < lon_w;
    self.cross_primary_meridian == is_range_crossing_prim_meridian
      && self.lon_min <= lon_w
      && lon_e < self.lon_max
  }

  /// Do not accept points on the SE-NE great circle arc of the zone
  /// nor on the NW-NE small circle arc of the zone.
  pub fn contains(&self, lon: f64, lat: f64) -> bool {
    let is_lat_compatible = (self.lat_min..self.lat_max).contains(&lat);
    is_lat_compatible
      && if self.cross_primary_meridian {
        lon < self.lon_max || self.lon_min <= lon
      } else {
        (self.lon_min..self.lon_max).contains(&lon)
      }
  }

  /// Do not accept points on the SE-NE great circle arc of the zone
  /// nor on the NW-NE small circle arc of the zone.
  /// Do not accept the given point if it is on a border.
  pub fn contains_exclusive(&self, lon: f64, lat: f64) -> bool {
    let is_lat_compatible = self.lat_min < lat && lat < self.lat_max;
    is_lat_compatible
      && if self.cross_primary_meridian {
        lon < self.lon_max || self.lon_min < lon
      } else {
        self.lon_min < lon && lon < self.lon_max
      }
  }

  /*
  /// Do accept points on both the SE-NE great circle arc, and on
  /// the NW-NE small circle arc of the zone.
  pub fn contains_inclusive(&self, lon: f64, lat: f64) -> bool {
    let is_lat_compatible = (self.lat_min..=self.lat_max).contains(&lat);
    is_lat_compatible
      && if self.cross_primary_meridian {
        lon <= self.lon_max || self.lon_min <= lon
      } else {
        (self.lon_min..=self.lon_max).contains(&lon)
      }
  }
  */

  /// Returned vertices order:
  /// * 0: SW
  /// * 1: NW
  /// * 2: NE
  /// * 3: SE
  pub fn vertices(&self) -> [(f64, f64); 4] {
    [
      (self.lon_min, self.lat_min),
      (self.lon_min, self.lat_max),
      (self.lon_max, self.lat_max),
      (self.lon_max, self.lat_min),
    ]
  }
}
