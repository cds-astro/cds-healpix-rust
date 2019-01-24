#[macro_use] extern crate structopt;

use structopt::*;

use cdshealpix::nested::{get_or_create};

#[derive(Debug, StructOpt)]
#[structopt(name = "hpx")]
#[structopt(raw(setting = "clap::AppSettings::ColoredHelp"))]
/// CDS Healpix library standalone
/// 
/// Example: ./hpx hash 8 256.369854 -12.365984
enum Args {
  #[structopt(name = "hash")]
  #[structopt(raw(setting = "clap::AppSettings::AllowLeadingHyphen"))]
  #[structopt(raw(setting = "clap::AppSettings::ColoredHelp"))]
  /// Compute the cell number of the given position at the given depth
  Hash {
    /// Healpix depth
    depth: u8,
    /// longitude, in degrees
    lon_deg: f64,
    /// latitude in degrees
    lat_deg: f64,
  }
}

fn main() {
  // println!("INFO: is this bin use BMI instruction: {}", is_x86_feature_detected!("bmi2"));
  //Args::allow_hyphen_values(true);
  match Args::from_args() {
    Args::Hash{depth, lon_deg, lat_deg} => {
      println!("{}", get_or_create(depth).hash(lon_deg.to_radians(), lat_deg.to_radians()));
    },
  };
}

