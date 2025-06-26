use crate::nested::Layer;
use log::error;

#[cfg(feature = "memmap")]
pub mod cindex;
#[cfg(feature = "sort")]
mod feature_gate; // get_hpx and get_hpx_opt don't rely on any heavy dependencies and are needed elsewhere in the code. To all external viewers, this module is inlined.

/// Returns a function computing given layer depth Healpix hash values from an input line of
/// `separator` separated fields in which the field at index `ilon` contains the longitude (in decimal degrees)
/// and the field at index `ilat` contains the latitude (also in decimal degrees).
/// WARNING: the returned hash value is `0` in case of error.
pub fn get_hpx(
  ilon: usize,
  ilat: usize,
  separator: char,
  layer: &'static Layer,
) -> impl Fn(&String) -> u64 {
  let (index_first_col, offset_to_second_col, ilon, ilat) = if ilon < ilat {
    (ilon, ilat - ilon - 1, 0, 1)
  } else {
    (ilat, ilon - ilat - 1, 1, 0)
  };

  move |line: &String| {
    let mut field_it = line.as_str().split(separator);
    let coos = [
      field_it.nth(index_first_col).and_then(|s| {
        s.parse::<f64>()
          .map_err(|e| error!("Error parsing 1st coo: '{}': {:?}", s, e))
          .ok()
      }),
      field_it.nth(offset_to_second_col).and_then(|s| {
        s.parse::<f64>()
          .map_err(|e| error!("Error parsing 2nd coo: '{}': {:?}", s, e))
          .ok()
      }),
    ];
    match (coos[ilon], coos[ilat]) {
      (Some(lon), Some(lat)) => layer.hash(lon.to_radians(), lat.to_radians()),
      _ => {
        error!(
          "Error parsing coordinates at line: {}. Hash set to 0.",
          line
        );
        0
      }
    }
  }
}

/// Returns a function computing given layer depth Healpix hash values from an input line of
/// `separator` separated fields in which the field at index `ilon` contains the longitude (in decimal degrees)
/// and the field at index `ilat` contains the latitude (also in decimal degrees).
/// WARNING: returns None in case of error.
pub fn get_hpx_opt(
  ilon: usize,
  ilat: usize,
  separator: char,
  layer: &'static Layer,
) -> impl Fn(&String) -> Option<u64> {
  let (index_first_col, offset_to_second_col, ilon, ilat) = if ilon < ilat {
    (ilon, ilat - ilon - 1, 0, 1)
  } else {
    (ilat, ilon - ilat - 1, 1, 0)
  };
  move |line: &String| {
    let mut field_it = line.as_str().split(separator);
    let coos = [
      field_it.nth(index_first_col).and_then(|s| {
        s.parse::<f64>()
          .map_err(|e| error!("Error parsing 1st coo: '{}': {:?}", s, e))
          .ok()
      }),
      field_it.nth(offset_to_second_col).and_then(|s| {
        s.parse::<f64>()
          .map_err(|e| error!("Error parsing 2nd coo: '{}': {:?}", s, e))
          .ok()
      }),
    ];
    match (coos[ilon], coos[ilat]) {
      (Some(lon), Some(lat)) => Some(layer.hash(lon.to_radians(), lat.to_radians())),
      _ => {
        error!(
          "Error parsing coordinates at line: {}. Lon: {:?}. Lat: {:?}. Hash set to 0.",
          line, coos[ilon], coos[ilat]
        );
        None
      }
    }
  }
}

#[cfg(feature = "sort")]
pub use feature_gate::*;