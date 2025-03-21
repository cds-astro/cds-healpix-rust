//! Code common to `map` and `mom` view.

use std::{ops::RangeInclusive, str::FromStr};

use clap::Subcommand;

#[derive(Debug, Subcommand)]
pub enum Mode {
  /// Creates an allsky view using the Mollweide projection.
  #[clap(name = "allsky")]
  AllSky {
    /// Number of pixels along the y-axis
    y_size: u16,
  },
  // TODO: add possibility to create a PNG corresponding to a given HEALPix cell!! (=> on-the-fly low res HiPS)
  /*/// Generate either a Mollweide or a Sinus projection centered on the mean non-empty center with
  /// automatic bounds.
  #[clap(name = "auto")]
  Auto {
    // Build the mom, adding values. Compute the weighted mean!
    /// Number of pixels along the y-axis
    y_size: u16,
  },*/
  /// Full control on the visualization (projection, center, bounds, ...)
  #[clap(name = "custom")]
  Custom {
    /// The chosen projection: 'car', 'cea', 'cyp', 'mer', 'hpx', 'ait', 'mol', 'par', 'sfl', 'sin'
    /// 'air', 'arc', 'feye', 'sin', 'stg', 'tan', 'zea'.
    proj: String,
    /// Size of the image along the x-axis, in pixels
    img_size_x: u16,
    /// Size of the image along the y-axis, in pixels
    img_size_y: u16,
    /// Longitude of the center of the projection, in degrees
    #[clap(short = 'l', long = "lon", default_value = "0.0")]
    center_lon: f64,
    /// Latitude of the center of the projection, in degrees
    #[clap(
      allow_negative_numbers = true,
      short = 'b',
      long = "lat",
      default_value = "0.0"
    )]
    center_lat: f64,
    /// Bounds, in the projection plane, matching both image edges along the the x-axis. E.g. [-0.2..0.2].
    #[clap(long = "x-bounds")]
    proj_bounds_x: Option<Bound>,
    /// Bounds, in the projection plane, matching both image edges along the the y-axis. E.g. [-0.2..0.2].
    #[clap(long = "y-bounds")]
    proj_bounds_y: Option<Bound>,
  },
}

#[derive(Clone, Debug)]
pub struct Bound(pub RangeInclusive<f64>);
impl FromStr for Bound {
  type Err = String;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    s.strip_prefix('[')
      .ok_or_else(|| format!("Invalid bound: {}. Must start with '['", s))
      .and_then(|s| {
        s.strip_suffix(']')
          .ok_or_else(|| format!("Invalid bound: {}. Must start with '['", s))
      })
      .and_then(|s| {
        s.split_once("..")
          .ok_or_else(|| format!("'..' separator not found in '{}'", s))
      })
      .and_then(|(l, r)| {
        l.parse::<f64>().map_err(|e| e.to_string()).and_then(|l| {
          r.parse::<f64>()
            .map(|r| Bound(l..=r))
            .map_err(|e| e.to_string())
        })
      })
  }
}
