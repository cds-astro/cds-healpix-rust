//! Module containing MOM implementations in sub-modules

pub mod bslice;
pub mod ported_from_mom_builder;
pub mod zvec;

#[cfg(all(test, feature = "skymap", feature = "memmap"))]
mod tests {
  use std::f64::consts::PI;

  use itertools::Itertools;

  use mapproj::pseudocyl::mol::Mol;

  use crate::nested::{
    map::{
      img::{to_mom_png_file, ColorMapFunctionType, PosConversion},
      mom::{
        impls::{bslice::FITSMom, zvec::MomVecImpl},
        Mom, WritableMom, ZUniqHashT,
      },
      skymap::SkyMapEnum,
    },
    n_hash,
  };

  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_skymap_to_mom_tofromfits_topng() {
    let img_size = (1366, 768);

    let path = "test/resources/skymap/skymap.fits";
    // Read map
    let count_map = SkyMapEnum::from_fits_file(path)
      .unwrap()
      .to_count_map_u32()
      .unwrap();
    // Chi2 merge to MOM
    let count_mom = count_map.to_chi2_mom(16.266);
    // Count MOM to Dens MOM
    let dens_mom = MomVecImpl::<u32, f64>::from_map(count_mom, |z, count| {
      let one_over_area = (n_hash(u32::depth_from_zuniq(z)) >> 2) as f64 / PI;
      count as f64 * one_over_area
    });
    // Write MOM in FITS
    dens_mom
      .to_fits_file("test/resources/skymap/skymap.mom.density.fits", "density")
      .unwrap();
    // Read MOM
    match FITSMom::from_fits_file("test/resources/skymap/skymap.mom.density.fits").unwrap() {
      FITSMom::U32F64(fits_mom) => {
        // Make PNG from MOM
        to_mom_png_file::<'_, _, Mol, _>(
          &fits_mom.get_mom(),
          img_size,
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/skymap.mom.density.png",
          false,
        )
        .unwrap();

        let mut n_comp = 0;
        for (i, (e1, e2)) in dens_mom
          .owned_entries()
          .zip_eq(fits_mom.get_mom().owned_entries())
          .enumerate()
        {
          n_comp += 1;
          assert_eq!(e1, e2);
          assert_eq!(e1.0, fits_mom.get_mom().get_z(i));
          assert_eq!(e1, fits_mom.get_mom().get_entry(i));
        }
        println!("n_comp = {}", n_comp);
        assert_eq!(n_comp, fits_mom.get_mom().len())
      }
      _ => panic!(),
    }
  }

  /* Do not run on github!
  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_skymap_to_mom_tofromfits_topng_big() {
    let img_size = (1366, 768);

    let path = "test/resources/skymap/gaiadr3.skymap.order10.fits";
    // Read map
    let count_map = SkyMapEnum::from_fits_file(path)
      .unwrap()
      .to_count_map_u32()
      .unwrap();
    /*let dens_map = count_map.to_dens_map();
    to_skymap_png_file::<'_, _, Mol, _>(
      dens_map.as_implicit_skymap_array(),
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), // LinearLog
      "test/resources/skymap/gaiadr3.skymap.density.order10.fits",
      false,
    )
    .unwrap();
    let dens_mom_1 = dens_map.to_chi2_mom();
    to_mom_png_file::<'_, _, Mol, _>(
      &dens_mom_1,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
      "test/resources/skymap/gaiadr3.mom.density.order10.v1.png",
      false,
    )
    .unwrap();*/
    // Chi2 merge to MOM
    let count_mom = count_map.to_chi2_mom();
    // Count MOM to Dens MOM
    let dens_mom = MomVecImpl::<u32, f64>::from_map(count_mom, |z, count| {
      let one_over_area = (n_hash(u32::depth_from_zuniq(z)) >> 2) as f64 / PI;
      count as f64 * one_over_area
    });
    /*to_mom_png_file::<'_, _, Mol, _>(
      &dens_mom,
      img_size,
      None,
      None,
      None,
      Some(PosConversion::EqMap2GalImg),
      None,
      Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
      "test/resources/skymap/gaiadr3.mom.density.order10.v2.png",
      false,
    )
    .unwrap();*/
    // Write MOM in FITS
    dens_mom
      .to_fits_file(
        "test/resources/skymap/gaiadr3.mom.density.order10.fits",
        "density",
      )
      .unwrap();
    // Read MOM
    match FITSMom::from_fits_file("test/resources/skymap/gaiadr3.mom.density.order10.fits").unwrap()
    {
      FITSMom::U32F64(fits_mom) => {
        // Make PNG from MOM
        to_mom_png_file::<'_, _, Mol, _>(
          &fits_mom.get_mom(),
          img_size,
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/gaiadr3.mom.density.order10.v3.png",
          false,
        )
        .unwrap();

        /*let dens_mom2: MomVecImpl<u32, f64> = fits_mom.get_mom().into();
        to_mom_png_file::<'_, _, Mol, _>(
          &dens_mom2,
          img_size,
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/gaiadr3.mom.density.order10.v4.png",
          false,
        )
        .unwrap();*/

        let mut n_comp = 0;
        for (i, (e1, e2)) in dens_mom
          .owned_entries()
          .zip_eq(fits_mom.get_mom().owned_entries())
          .enumerate()
        {
          n_comp += 1;
          assert_eq!(e1, e2);
          assert_eq!(e1.0, fits_mom.get_mom().get_z(i));
          assert_eq!(e1, fits_mom.get_mom().get_entry(i));
        }
        println!("n_comp = {}", n_comp);
        assert_eq!(n_comp, fits_mom.get_mom().len())
      }
      _ => panic!(),
    }

    /// Order 11
    let path = "test/resources/skymap/gaiadr3.skymap.order11.fits";
    // Read map
    let count_map = SkyMapEnum::from_fits_file(path)
      .unwrap()
      .to_count_map_u32()
      .unwrap();
    // Chi2 merge to MOM
    let count_mom = count_map.to_chi2_mom();
    // Count MOM to Dens MOM
    let dens_mom = MomVecImpl::<u32, f64>::from_map(count_mom, |z, count| {
      let one_over_area = (n_hash(u32::depth_from_zuniq(z)) >> 2) as f64 / PI;
      count as f64 * one_over_area
    });
    // Write MOM in FITS
    dens_mom
      .to_fits_file(
        "test/resources/skymap/gaiadr3.mom.density.order11.fits",
        "density",
      )
      .unwrap();
    dens_mom
      .to_fits_bintable_file(
        "test/resources/skymap/gaiadr3.mom.density.order11.bintable.fits",
        "density",
      )
      .unwrap();
    // Read MOM
    match FITSMom::from_fits_file("test/resources/skymap/gaiadr3.mom.density.order11.fits").unwrap()
    {
      FITSMom::U32F64(fits_mom) => {
        // Make PNG from MOM
        to_mom_png_file::<'_, _, Mol, _>(
          &fits_mom.get_mom(),
          img_size,
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/gaiadr3.mom.density.order11.png",
          false,
        )
        .unwrap();
      }
      _ => panic!(),
    }
  }*/
}
