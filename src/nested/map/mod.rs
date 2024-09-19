pub mod astrometry;
pub mod fits;
pub mod img;
pub mod skymap;

#[cfg(test)]
mod tests {
  use super::{
    img::{ColorMapFunctionType, PosConversion},
    skymap::SkyMap,
  };
  use mapproj::pseudocyl::mol::Mol;

  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_skymap() {
    let path = "test/resources/skymap/skymap.fits";
    let skymap = SkyMap::from_fits_file(path).unwrap();
    skymap
      .to_png_file::<Mol, _>(
        (1600, 800),
        None,
        None,
        None,
        Some(PosConversion::EqMap2GalImg),
        None,
        Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
        "test/resources/skymap/skymap.png",
        false,
      )
      .unwrap();
  }
}
