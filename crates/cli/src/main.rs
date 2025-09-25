use std::{
  error::Error,
  io::{StdoutLock, Write},
};

use clap::Parser;

use hpxlib::{proj, unproj};

use hpx_cli::{
  coverage::Coverage,
  hcidx::HealpixCumulIndex,
  input::{coo::CooInput, pos::PosInput},
  map::Map,
  mom::Mom,
  nested::Nested,
  qhcidx::QueryHCIndex,
  sort::Sort,
};

// Avoid musl's default allocator due to lackluster performance
// https://nickb.dev/blog/default-musl-allocator-considered-harmful-to-performance
#[cfg(all(target_env = "musl", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Perform HEALPix related operations on the command line.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
enum Args {
  /// Computes the projected coordinates (x, y) ∊ ([0..8[, [-2..2]) of input equatorial coordinates
  #[command(subcommand)]
  Proj(PosInput),
  /// Computes the equatorial coordinates of Healpix projected coordinates (x, y) ∊ ([0..8[, [-2..2])
  #[command(subcommand)]
  Unproj(CooInput),
  /// Operations in the HEALPix NESTED scheme
  #[command(subcommand)]
  Nested(Nested),
  Sort(Sort),
  #[clap(name = "hcidx")]
  HCIdx(HealpixCumulIndex), // Mom fromskymap
  #[clap(name = "qhcidx")]
  QHCIdx(QueryHCIndex),
  #[command(subcommand)]
  Map(Map),
  #[command(subcommand)]
  Mom(Mom),
  Cov(Coverage),
  // BMOC op/view/convert(bintable,CountMom?)
}

impl Args {
  fn exec(self) -> Result<(), Box<dyn Error>> {
    match self {
      Self::Proj(e) => e.exec(
        |lon_rad: f64,
         lat_rad: f64,
         write: &mut StdoutLock<'static>,
         sep: char|
         -> Result<(), Box<dyn Error>> {
          let (x, y) = proj(lon_rad, lat_rad);
          write!(write, "{:17.15}{}{:+017.15}", x, sep, y).map_err(|e| e.into())
        },
      ),
      Self::Unproj(e) => e.exec(
        |x: f64,
         y: f64,
         write: &mut StdoutLock<'static>,
         sep: char|
         -> Result<(), Box<dyn Error>> {
          let (lon_rad, lat_rad) = unproj(x, y);
          write!(
            write,
            "{:017.13}{}{:+017.13}",
            lon_rad.to_degrees(),
            sep,
            lat_rad.to_degrees()
          )
          .map_err(|e| e.into())
        },
      ),
      Self::Nested(e) => e.exec(),
      Self::Sort(e) => e.exec(),
      Self::HCIdx(e) => e.exec(),
      Self::QHCIdx(e) => e.exec(),
      Self::Map(e) => e.exec(),
      Self::Mom(e) => e.exec(),
      Self::Cov(e) => e.exec(),
    }
  }
}

fn main() -> Result<(), Box<dyn Error>> {
  env_logger::init();
  let args = Args::parse();
  args.exec()
}
