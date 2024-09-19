use super::math::*;

/// Computed using the ZYZ Euler transformation from the FK4 B1950 frame in which, based on
/// Blaauw (1958) and Liu (2011, A&A536, A16):
/// * North Galactic Pole FK4 B1950 longitude: 12h49m= +192.25 deg (1st rotation: z-axis of 192.25 deg)
/// * North Galactic Pole FK4 B1950 latitude : +27.4 deg (2nd rotation: x'-axis by 90 - 27.4 = 62.6 deg)
/// * Position angle, at the Galactic pole, from the equatorial pole,
///   of the great circle passing through the Galactic center: +123 deg.
///   Noting NEP (North Equatorial Pole), NGP (North Galactic Pole) and GC (Galactic Center), it is
///   the angle NEP-NGP-GP. The rotation from NGP so set the origin of the frame at the GC is thus
///   (3rd rotation: z''-axis by 180 - 123 = 57 deg).
/// # Note
/// * The additional precision comes from `bc/src/galactic.bc` and have been cross-checked
///   with a code using 80bits floats from F. Ochsenbein
const FK4B19502GAL: M3x3 = M3x3(
  // Matrix from BC
  XYZt(
    -0.066_988_739_415_150_72,
    -0.872_755_765_851_992_7,
    -0.483_538_914_632_184_24,
  ),
  XYZt(
    0.492_728_466_075_323_6,
    -0.450_346_958_019_961_33,
    0.744_584_633_283_031,
  ),
  XYZt(
    -0.867_600_811_151_434_8,
    -0.188_374_601_722_920_45,
    0.460_199_784_783_851_65,
  ),
);

/// Computed using the ZYZ Euler transformation from the **truncated** (ESA convention) coordinates
/// of the North Galactic Pole (NGP) plus the position angle of the NGP with respect to the
/// equatorial pole. See:
/// * Section 7.1 of Murray (1989),
/// * Eq. (1.5.9) and (1.5.10) of Hipparcos and Tycho Vol1,
/// * Eq. (2) of Liu (2011, A&A 536, 102).
/// * Section 2 from Liu (2011, A&A536, A16)
/// * Eq. (31) in Fox document
/// That is:
/// * NGP longitude in FK5 J2000: 12h51m26.2755s = +192.85948 deg (1st rotation: z-axis of 192.85948 deg)
/// * NGP latitude  in FK5 J2000: +27.12825 deg (2nd rotation: x'-axis by 90 - 27.12825 = 62.87175 deg)
/// * Position angle, at the Galactic pole, from the equatorial pole,
///   of the great circle passing through the Galactic center: +122.93192 deg.
///   Noting NEP (North Equatorial Pole), NGP (North Galactic Pole) and GC (Galactic Center), it is
///   the angle NEP-NGP-GP. The rotation from NGP so set the origin of the frame at the GC is thus
///   (3rd rotation: z''-axis by 180 - 122.93192 = 57.06808 deg).
/// # Notes
/// * follows the ESA convention defined for Hipparcos, Tycho and Gaia
///     + no difference made between FK5 and ICRS
///     + based on truncated values computed from the conversion of the galactic center coordinates
///       (and position angle) from FK4 to FK5
/// * The additional precision comes from `bc/src/galactic.bc` and have been cross-checked
///   with SOFA (`icrs2g.c` routine).
const FK5J20002GAL_ESA: M3x3 = M3x3(
  XYZt(
    -0.054_875_560_416_215_368,
    -0.873_437_090_234_885,
    -0.483_835_015_548_713_2,
  ),
  XYZt(
    0.494_109_427_875_583_65,
    -0.444_829_629_960_011_2,
    0.746_982_244_497_218_8,
  ),
  XYZt(
    -0.867_666_149_019_004_7,
    -0.198_076_373_431_201_52,
    0.455_983_776_175_066_9,
  ),
);

// So far, no reliable and MIT high precision crate in rust:
// * rug is LGPG
// * sin of astro-float does not seems reliable enough
// Thus I decided to write a 'bc' script to compute the matrix coefficient

/// Accounts from the difference between FK5 and ICRS in Mignard (2002).
/// * Eq. (18) from Liu (2011) A&A536, A16
/// * Eq. (30) in Fox document
const ICRS2GAL: M3x3 = M3x3(
  XYZt(
    -0.054_875_657_712_619_678_1,
    -0.873_437_051_955_779_129_8,
    -0.483_835_073_616_418_380_3,
  ),
  XYZt(
    0.494_109_437_197_107_641_2,
    -0.444_829_721_222_053_763_5,
    0.746_982_183_839_845_094_133,
  ),
  XYZt(
    -0.867_666_137_557_162_561_5,
    -0.198_076_337_275_070_594_6,
    0.455_983_813_691_152_347_6,
  ),
);

/// The Galactic coordinate system is barycentric, its xy-plane is the galactic plane (it is assumed
/// the the barycenter of the solar system is in the galactic plane, which is not correct)
/// and its x-axis points toward the galactic center.
/// It is defined from Equatorial frames (such a FK4, FK5 or ICRS) by the position of the
/// North Galactic Pole (NGC) in the equatorial frame, plus the position angle formed by
/// the North Equatorial Pole (NEP), the NGC and the great circle passing through the NEP and
/// the Galactic Center (GC).
pub struct Galactic {
  /// Rotation matrix
  m: M3x3,
}

impl Galactic {
  /// New Galactic frame for the conversion with the FK4 B19500 frame.
  pub const fn new_for_fk4_b1950() -> Self {
    Galactic { m: FK4B19502GAL }
  }

  /// New Galactic frame for the conversion with the FK5 J2000 frame, or the ICRS frame
  /// following the ESA (Hipparcos, Tycho, Gaia) conventions.
  /// Its also correspond the SOFA icrs <-> gal conversion.
  pub const fn new_for_fk5_j2000_and_icrs() -> Self {
    Galactic {
      m: FK5J20002GAL_ESA,
    }
  }

  /// New Galactic frame for the conversion with the ICRS frame accounting for the Mignard (2002)
  /// offsets between FK5 and ICRS.
  pub const fn new_for_icrs_including_fk5_icrs_offsets_from_mignard2002() -> Self {
    Galactic { m: ICRS2GAL }
  }

  /// Convert the given coordinate from Equatorial (FK4 B1950, FK5 J2000 or ICRS) to Galactic.  
  pub fn coo_eq2gal(&self, pos: &Coo) -> Coo {
    Coo::from_xyz(&self.xyz_eq2gal(&pos.xyz()))
  }

  /// Convert the given coordinate from Galactic to Equatorial (FK4 B1950, FK5 J2000 or ICRS).
  pub fn coo_gal2eq(&self, pos: &Coo) -> Coo {
    Coo::from_xyz(&self.xyz_gal2eq(&pos.xyz()))
  }

  /// Convert the given Euclidean position from Equatorial (FK4 B1950, FK5 J2000 or ICRS) to
  /// Galactic.  
  pub fn xyz_eq2gal(&self, pos: &XYZ) -> XYZ {
    self.m.rotate(pos)
  }

  /// Convert the given Euclidean position from Galactic to Equatorial
  /// (FK4 B1950, FK5 J2000 or ICRS).
  pub fn xyz_gal2eq(&self, pos: &XYZ) -> XYZ {
    self.m.unrotate(pos)
  }
}

#[cfg(test)]
mod tests {
  use crate::nested::map::astrometry::{gal::Galactic, math::Coo};

  #[test]
  fn test_gal2eq() {
    let lon = 0.0;
    let lat = 0.0;
    let gal_center = Galactic::new_for_fk5_j2000_and_icrs().coo_gal2eq(&Coo::new(lon, lat));
    println!(
      "Gal Center: ({:.6}, {:.6})",
      gal_center.lon.to_degrees(),
      gal_center.lat.to_degrees()
    )
  }
}
