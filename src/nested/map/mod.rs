use std::fmt::{Debug, Display};

use num_traits::{AsPrimitive, PrimInt};

pub mod astrometry;
pub mod fits;
#[cfg(feature = "skymap")]
pub mod img;
pub mod mom;
pub mod skymap;

/// `HHash` stands for HEALPix Hash.
pub trait HHash:
  // 'static mean that Idx does not contains any reference
  'static + PrimInt + AsPrimitive<usize> + Send + Sync + Debug + Display + Clone
{
  fn to_u64(&self) -> u64;
  fn from_u64(v: u64) -> Self;
  fn from_usize(v: usize) -> Self;
}

impl HHash for u32 {
  fn to_u64(&self) -> u64 {
    *self as u64
  }
  fn from_u64(v: u64) -> Self {
    v as u32
  }
  fn from_usize(v: usize) -> Self {
    v as u32
  }
}
impl HHash for u64 {
  fn to_u64(&self) -> u64 {
    *self
  }
  fn from_u64(v: u64) -> Self {
    v
  }
  fn from_usize(v: usize) -> Self {
    v as u64
  }
}

#[cfg(all(test, feature = "skymap"))]
mod tests {
  use super::{
    img::{ColorMapFunctionType, PosConversion},
    skymap::SkyMapEnum,
  };
  use mapproj::pseudocyl::mol::Mol;

  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_skymap_spec() {
    let path = "test/resources/skymap/skymap.2mass.depth6.fits";
    // let path = "test/resources/skymap/gaiadr2.skymap.order10.fits";
    let skymap = SkyMapEnum::from_fits_file(path).unwrap();
    skymap
      .to_skymap_png_file::<Mol, _>(
        (1600, 800),
        None,
        None,
        None,
        Some(PosConversion::EqMap2GalImg),
        None,
        Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
        "test/resources/skymap/skymap.2mass.depth6.png",
        false,
      )
      .unwrap();

    let path = "test/resources/skymap/skymap.gaiadr3.depth6.fits";
    // let path = "test/resources/skymap/gaiadr2.skymap.order10.fits";
    let skymap = SkyMapEnum::from_fits_file(path).unwrap();
    skymap
      .to_skymap_png_file::<Mol, _>(
        (1600, 800),
        None,
        None,
        None,
        Some(PosConversion::EqMap2GalImg),
        None,
        Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
        "test/resources/skymap/skymap.gaiadr3.depth6.png",
        false,
      )
      .unwrap();

    let path = "test/resources/skymap/hats_test.fits";
    let skymap = SkyMapEnum::from_fits_file(path).unwrap();
    skymap
      .to_skymap_png_file::<Mol, _>(
        (1600, 800),
        None,
        None,
        None,
        Some(PosConversion::EqMap2GalImg),
        None,
        Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
        "test/resources/skymap/hats_test.png",
        false,
      )
      .unwrap();
  }

  fn init_logger() {
    let log_level = log::LevelFilter::max();
    // let log_level = log::LevelFilter::Error;

    let _ = env_logger::builder()
      // Include all events in tests
      .filter_level(log_level)
      // Ensure events are captured by `cargo test`
      .is_test(true)
      // Ignore errors initializing the logger if tests race to configure it
      .try_init();
  }

  #[test]
  fn test_read_skymap() {
    init_logger();
    let path = "test/resources/skymap/gaiadr3.nside64.densmap.fits";
    let skymap = SkyMapEnum::from_fits_file(path).unwrap();
  }
}
