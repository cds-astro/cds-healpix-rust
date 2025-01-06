//! This code has been developed to be ported in GLSL in order to be used in WebGL (for AladinLite).
//! It computes HEALPix hash from Euclidean coordinates provided on single precision floats.

use std::f32::consts::PI;

use crate::nside;

pub const TRANSITION_Z: f32 = 2.0 / 3.0;
pub const TRANSITION_Z_INV: f32 = 3.0 / 2.0;

/// Returns the projection, in the 2D Euclidean plane, of the given position on the unit sphere.
/// # Inputs:
/// - `x`: in `[-1.0, 1.0]`
/// - `y`: in `[-1.0, 1.0]`
/// - `z`: in `[-1.0, 1.0]`
/// # Output
/// - `X`: coordinate along the X-axis in the projection plane, in `[-4, 4]`
/// - `Y`: coordinate along the Y-axis in the projection plane, in `[-2, 2]`
/// # Remark
/// For the HPX projection as defined in Calabreta, use:
///   - `X *= PI / 4`
///   - `Y *= PI / 4`
/// # WARNING
/// - The function assumes, without checking, that the input vector is a unit vector
///   (hence `x^2 + y^2 + z^2 = 1`) !!
pub fn proj_gpu(x: f32, y: f32, z: f32) -> (f32, f32) {
  assert!((-1.0..=1.0).contains(&x));
  assert!((-1.0..=1.0).contains(&y));
  assert!((-1.0..=1.0).contains(&z));
  debug_assert!(
    1.0 - (x * x + y * y + z * z) < 1e-5,
    "{}",
    1.0 - (x * x + y * y + z * z)
  );
  if z > TRANSITION_Z {
    // North polar cap, Collignon projection.
    let (x_pm1, offset) = xpm1_and_offset(x, y);
    let sqrt_3_one_min_z = (3.0 * one_minus_z_pos(x, y, z)).sqrt();
    (
      (x_pm1 * sqrt_3_one_min_z) + offset as f32,
      2.0 - sqrt_3_one_min_z,
    )
  } else if z < -TRANSITION_Z {
    // South polar cap, Collignon projection
    let (x_pm1, offset) = xpm1_and_offset(x, y);
    let sqrt_3_one_min_z = (3.0 * one_minus_z_neg(x, y, z)).sqrt();
    (
      (x_pm1 * sqrt_3_one_min_z) + offset as f32,
      -2.0 + sqrt_3_one_min_z,
    )
  } else {
    // Equatorial region, Cylindrical equal area projection
    (y.atan2(x) * 4.0 / PI, z * TRANSITION_Z_INV)
  }
}

/// Returns the position in on the unit sphere `(x, y, z)` of the give position in the HEALPix
/// 2D Euclidean projection plane.
/// # Inputs
/// - `X`: coordinate along the X-axis in the projection plane, in `[-4, 4]`
/// - `Y`: coordinate along the Y-axis in the projection plane, in `[-2, 2]`
/// # Output:
/// - `x`: in `[-1.0, 1.0]`
/// - `y`: in `[-1.0, 1.0]`
/// - `z`: in `[-1.0, 1.0]`
/// # Remark
/// From the HPX projection as defined in Calabreta, use:
///   - `X /= PI / 4`
///   - `Y /= PI / 4`
pub fn unproj_gpu(x: f32, y: f32) -> (f32, f32, f32) {
  use std::f32::consts::{FRAC_PI_2, FRAC_PI_4};
  const ONE_OVER_SQRT6: f32 = 0.408_248_3_f32;
  assert!((-2f32..=2f32).contains(&y)); // check y
  if y > 1.0 {
    // North Polar Cap
    assert!((-4f32..=4f32).contains(&x));
    let x = if x < 0.0 { 8.0 + x } else { x };
    debug_assert!(x > 0.0);
    let OffsetAndPM1 { offset, pm1 } = pm1_offset_decompose(x);
    let sqrt_of_three_time_one_minus_sin_of = 2.0 - y;
    let x = if sqrt_of_three_time_one_minus_sin_of > 1e-6 {
      deal_with_numerical_approx_in_edges(pm1 / sqrt_of_three_time_one_minus_sin_of)
    } else {
      pm1
    } + offset as f32;
    // It would be faster, but less accurate, to use:
    // let z = 1.0 - sqrt_of_three_time_one_minus_sin_of.pow2() / 3.0;
    // let cos_lat = sqrt(1 - z^2);
    let lat = 2.0 * f32::acos(sqrt_of_three_time_one_minus_sin_of * ONE_OVER_SQRT6) - FRAC_PI_2;
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = (x * FRAC_PI_4).sin_cos();
    (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat)
  } else if y < -1.0 {
    // South polar cap
    let x = if x < 0.0 { 8.0 + x } else { x };
    debug_assert!(x > 0.0);
    let OffsetAndPM1 { offset, pm1 } = pm1_offset_decompose(x);
    let sqrt_of_three_time_one_minus_sin_of = 2.0 + y;
    let x = if sqrt_of_three_time_one_minus_sin_of > 1e-6 {
      deal_with_numerical_approx_in_edges(pm1 / sqrt_of_three_time_one_minus_sin_of)
    } else {
      pm1
    } + offset as f32;
    // It would be faster, but less accurate, to use:
    // let z = -1.0 + sqrt_of_three_time_one_minus_sin_of.pow2() / 3.0;
    // let cos_lat = sqrt(1 - z^2);
    let lat = FRAC_PI_2 - 2.0 * f32::acos(sqrt_of_three_time_one_minus_sin_of * ONE_OVER_SQRT6);
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = (x * FRAC_PI_4).sin_cos();
    (cos_lon * cos_lat, sin_lon * cos_lat, sin_lat)
  } else {
    // Equatorial region
    let z = y * TRANSITION_Z; // = sin(lat)
    let cos_lat = if z < 1e-2 {
      // sqrt(1 - x²) = 1 - x²/2 - x⁴/8 - x⁶/16
      let tmp = 0.5 * z * z;
      1.0 - tmp - 0.5 * tmp * tmp
    } else {
      (1.0 - z * z).sqrt()
    };
    let (sin_lon, cos_lon) = (x * FRAC_PI_4).sin_cos();
    (cos_lon * cos_lat, sin_lon * cos_lat, z)
  }
}

// Decompose the given positive real value in
// --* an integer offset in [1, 3, 5, 7] (*PI/4) and
// --* a real value in [-1.0, 1.0] (*PI/4)
pub(crate) struct OffsetAndPM1 {
  offset: u8, // = 1, 3, 5 or 7
  pm1: f32,   // in [-1.0, 1.0]
}
#[inline]
pub(crate) fn pm1_offset_decompose(x: f32) -> OffsetAndPM1 {
  let floor: u8 = x as u8;
  let odd_floor: u8 = floor | 1u8;
  OffsetAndPM1 {
    offset: odd_floor & 7u8, // value modulo 8 = 1/3/5/7
    pm1: x - (odd_floor as f32),
  }
}

fn deal_with_numerical_approx_in_edges(x: f32) -> f32 {
  x.clamp(-1.0, 1.0)
}

/// Returns the cell number (hash value) associated with the given position on the unit sphere,
/// together with the offset `(dx, dy)` on the Euclidean plane of the projected position with
/// respect to the origin of the cell (South vertex).
/// # Inputs:
/// - `depth` in `[0, 14]` (so that and HEALPix cell number can be stored on an unsigned integer)
/// - `x`: in `[-1.0, 1.0]`
/// - `y`: in `[-1.0, 1.0]`
/// - `z`: in `[-1.0, 1.0]`
/// # Output
/// - the cell number (hash value) associated with the given position on the unit sphere,
///   in `[0, 12*nside^2[`
/// - `dx`: the positional offset $\in [0, 1[$ along the south-to-east axis
/// - `dy`: the positional offset $\in [0, 1[$ along the south-to-west axis
/// # WARNING
/// - The function assumes, without checking, that the input vector is a unit vector
///   (hence `x^2 + y^2 + z^2 = 1`) !!
/// - Operations being made on simple precision float, the precision is lower than `~0.2 arcsec` only!!
/// - At depth 13, the precision on `(dx, dy)` is better than `(1/512, 1/512)`, i.e. 2e-3.
#[allow(clippy::many_single_char_names)]
pub fn hash_with_dxdy(depth: u8, x: f32, y: f32, z: f32) -> (u32, f32, f32) {
  assert!(depth <= 14);
  assert!((-1.0..=1.0).contains(&x));
  assert!((-1.0..=1.0).contains(&y));
  assert!((-1.0..=1.0).contains(&z));
  // println!("norm: {}", (x *  x + y * y + z * z));
  debug_assert!(
    1.0 - (x * x + y * y + z * z) < 1e-5,
    "{}",
    1.0 - (x * x + y * y + z * z)
  );
  // A f32 mantissa contains 23 bits.
  // - it basically means that when storing (x, y) coordinates,
  //   we can go as deep as depth 24 (or maybe 25)
  let nside = nside(depth);
  let half_nside = nside as f32 * 0.5;
  let (x_pm1, q) = xpm1_and_q(x, y);
  let (d0h, x_in_d0c, y_in_d0c) = if z > TRANSITION_Z {
    // North polar cap, Collignon projection.
    // - set the origin to (PI/4, 0)
    let sqrt_3_one_min_z = (3.0 * one_minus_z_pos(x, y, z)).sqrt();
    let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
    let d0h = q;
    (d0h, x_proj, y_proj)
  } else if z < -TRANSITION_Z {
    // South polar cap, Collignon projection
    // - set the origin to (PI/4, -PI/2)
    let sqrt_3_one_min_z = (3.0 * one_minus_z_neg(x, y, z)).sqrt();
    let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, sqrt_3_one_min_z);
    let d0h = q + 8;
    (d0h, x_proj, y_proj)
  } else {
    // Equatorial region, Cylindrical equal area projection
    // - set the origin to (PI/4, 0)               if q = 2
    // - set the origin to (PI/4, -PI/2)           if q = 0
    // - set the origin to (0, -TRANSITION_LAT)    if q = 3
    // - set the origin to (PI/2, -TRANSITION_LAT) if q = 1
    // let zero_or_one = (x_cea as u8) & 1;
    let y_pm1 = z * TRANSITION_Z_INV;
    // |\2/|
    // .3X1.
    // |/0\|
    let q01 = (x_pm1 > y_pm1) as u8; /* 0/1 */
    debug_assert!(q01 == 0 || q01 == 1);
    let q12 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */
    debug_assert!(q12 == 0 || q12 == 1);
    let q03 = 1 - q12; /* 1\0 */
    //let q13 = q01 ^ q12;                              debug_assert!(q13 == 0 || q13 == 1);
    let q1 = q01 & q12; /* = 1 if q1, 0 else */
    debug_assert!(q1 == 0 || q1 == 1);
    // x: x_pm1 + 1 if q3 | x_pm1 - 1 if q1 | x_pm1 if q0 or q2
    let x_proj = x_pm1 - ((q01 + q12) as i8 - 1) as f32;
    // y: y_pm1 + 0 if q2 | y_pm1 + 1 if q1 or q3 | y_pm1 + 2 if q0
    let y_proj = y_pm1 + (q01 + q03) as f32;
    // d0h: +8 if q0 | +4 if q3 | +5 if q1
    let d0h = ((q01 + q03) << 2) + ((q + q1) & 3);
    (d0h, x_proj, y_proj)
  };
  // Coords inside the base cell
  let x = half_nside * (y_in_d0c + x_in_d0c);
  debug_assert!(
    (0.0 - 1e-5) < x && x < (nside as f32 + 1e-5),
    "x: {}, x_proj: {}; y_proj: {}",
    &x,
    &x_in_d0c,
    &y_in_d0c
  );
  let y = half_nside * (y_in_d0c - x_in_d0c);
  debug_assert!((0.0 - 1e-5) < x && y < (nside as f32 + 1e-5), "y: {}", &y);
  // Ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
  let mut i = x as u32;
  let mut j = y as u32;
  if i == nside {
    i -= 1;
  } // Deal with numerical inaccuracies, rare so branch miss-prediction negligible
  if j == nside {
    j -= 1;
  } // Deal with numerical inaccuracies, rare so branch miss-prediction negligible
  (
    ((d0h as u32) << (depth << 1)) | ij2z(i, j),
    (x - (i as f32)),
    (y - (j as f32)),
  )
}

fn xpm1_and_q(x: f32, y: f32) -> (f32, u8) {
  let x_neg = (x < 0.0) as u8;
  debug_assert!(x_neg <= 1);
  let y_neg = (y < 0.0) as u8;
  debug_assert!(y_neg <= 1);
  let q = (x_neg + y_neg) | (y_neg << 1);
  debug_assert!(y_neg <= 3);
  // The purpose is to have the same numerical precision for each base cell
  // by avoiding subtraction by 1.0 or 3.0 or 5.0 or 7.0
  let lon = y.abs().atan2(x.abs());
  debug_assert!((0.0..=PI / 2.0).contains(&lon));
  let x02 = lon * 4.0 / PI;
  debug_assert!((0.0..=2.0).contains(&x02));
  if x_neg != y_neg {
    // Could be replaced by a sign copy from (x_neg ^ y_neg) << 32
    (1.0 - x02, q)
  } else {
    (x02 - 1.0, q)
  }
}

fn xpm1_and_offset(x: f32, y: f32) -> (f32, i8) {
  let x_neg = (x < 0.0) as i8;
  debug_assert!(x_neg == 0 || x_neg == 1);
  let y_neg = (y < 0.0) as i8;
  debug_assert!(y_neg == 0 || y_neg == 1);
  // x>0, y>0 => [    0,  pi/2[ => offset =  1
  // x<0, y>0 => [pi/2 ,    pi[ => offset =  3
  // x<0, y<0 => [3pi/2,    pi[ => offset = -3
  // x>0, y<0 => [pi   , 3pi/2[ => offset = -1
  let offset = ((-y_neg) << 2) + 1 + ((x_neg ^ y_neg) << 1);
  let lon = y.abs().atan2(x.abs());
  debug_assert!((0.0..=PI / 2.0).contains(&lon));
  let x02 = lon * 4.0 / PI;
  debug_assert!((0.0..=2.0).contains(&x02));
  if x_neg != y_neg {
    // Could be replaced by a sign copy from (x_neg ^ y_neg) << 32
    (1.0 - x02, offset)
  } else {
    (x02 - 1.0, offset)
  }
}

fn one_minus_z_pos(x: f32, y: f32, z: f32) -> f32 {
  debug_assert!(z > 0.0);
  let d2: f32 = x * x + y * y; // z = sqrt(1 - d2) AND sqrt(1 - x) = 1 - x / 2 - x^2 / 8 - x^3 / 16 - 5 x^4/128 - 7 * x^5/256
                               /*if d2 < 1e-2 { // <=> dec > 84.27 deg
                                 d2 * (0.5 + d2 * (0.125 + 0.0625 * d2))
                               }*/
  if d2 < 1e-1 {
    // <=> dec > 84.27 deg
    d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
  } else {
    1.0 - z
  }
}

fn one_minus_z_neg(x: f32, y: f32, z: f32) -> f32 {
  debug_assert!(z < 0.0);
  let d2: f32 = x * x + y * y; // z = sqrt(1 - d2) AND sqrt(1 - x) = 1 - x / 2 - x^2 / 8 - x^3 / 16 - 5 x^4/128 - 7 * x^5/256
  if d2 < 1e-1 {
    // <=> dec < -84.27 deg
    // 0.5 * d2 + 0.125 * d2 * d2
    d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
  } else {
    z + 1.0
  }
}

/// Z-Order curve projection.
fn ij2z(mut i: u32, mut j: u32) -> u32 {
  i |= j << 16;
  j = (i ^ (i >> 8)) & 0x0000FF00_u32;
  i = i ^ j ^ (j << 8);
  j = (i ^ (i >> 4)) & 0x00F000F0_u32;
  i = i ^ j ^ (j << 4);
  j = (i ^ (i >> 2)) & 0x0C0C0C0C_u32;
  i = i ^ j ^ (j << 2);
  j = (i ^ (i >> 1)) & 0x22222222_u32;
  i = i ^ j ^ (j << 1);
  i
}

#[cfg(test)]
mod tests {

  use super::*;
  use crate::nested;
  use crate::proj;

  #[test]
  fn testok_hash_gpu_spe_1() {
    const DELTA_DEPTH: u8 = 1;
    // Input
    let depth = 13; // done with depth 13
                    // Computations
    let layer = nested::get(depth + DELTA_DEPTH);
    let h = 259907375;
    // for h in 0..layer.n_hash() {
    let (ra, de) = layer.center(h);

    println!("{:?}", proj(ra, de));

    let (sl, cl) = ra.sin_cos();
    let (sb, cb) = de.sin_cos();
    let (x, y, z) = (cb * cl, cb * sl, sb);
    let (nh, dx, dy) = hash_with_dxdy(depth, x as f32, y as f32, z as f32);
    assert_eq!((h >> (DELTA_DEPTH << 1)) as u32, nh);
    let prec = 1.0 / 512.0; //512.0; // HiPS images: 512 x 512
    match h & 3 {
      0 => {
        assert!(
          (dx - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.25).abs()
        ); // precision should depends on depth!!
        assert!(
          (dy - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.25).abs()
        ); // 1e-3 => images 1000x1000 while HiPS uses 512x512
      }
      1 => {
        assert!(
          (dx - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.75).abs()
        );
        assert!(
          (dy - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.25).abs()
        );
      }
      2 => {
        assert!(
          (dx - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.25).abs()
        );
        assert!(
          (dy - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.75).abs()
        );
      }
      3 => {
        assert!(
          (dx - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.75).abs()
        );
        assert!(
          (dy - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.75).abs()
        );
      }
      _ => unreachable!(),
    }
    // }
  }

  #[test]
  fn testok_hash_gpu_spe_2() {
    const DELTA_DEPTH: u8 = 1;
    // Input
    let depth = 13; // done with depth 13
                    // Computations
    let layer = nested::get(depth + DELTA_DEPTH);
    let h = 430712695;
    // for h in 0..layer.n_hash() {
    let (ra, de) = layer.center(h);

    println!("RA: {}; Dec: {}", ra.to_degrees(), de.to_degrees());
    println!("Proj: {:?}", proj(ra, de));

    let (sl, cl) = ra.sin_cos();
    let (sb, cb) = de.sin_cos();
    let (x, y, z) = (cb * cl, cb * sl, sb);
    let (nh, dx, dy) = hash_with_dxdy(depth, x as f32, y as f32, z as f32);
    // let (onh, odx, ody) = hash_with_dxdy_old(depth, x as f32, y as f32, z as f32);
    assert_eq!((h >> (DELTA_DEPTH << 1)) as u32, nh);
    let prec = 1.0 / 512.0; //512.0; // HiPS images: 512 x 512
    match h & 3 {
      0 => {
        assert!(
          (dx - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.25).abs()
        ); // precision should depends on depth!!
        assert!(
          (dy - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.25).abs()
        ); // 1e-3 => images 1000x1000 while HiPS uses 512x512
      }
      1 => {
        assert!(
          (dx - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.75).abs()
        );
        assert!(
          (dy - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.25).abs()
        );
      }
      2 => {
        assert!(
          (dx - 0.25).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.25).abs()
        );
        assert!(
          (dy - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.75).abs()
        );
      }
      3 => {
        assert!(
          (dx - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dx - 0.75).abs()
        );
        assert!(
          (dy - 0.75).abs() <= prec,
          "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
          h,
          ra.to_degrees(),
          de.to_degrees(),
          dx,
          dy,
          (dy - 0.75).abs()
        );
      }
      _ => unreachable!(),
    }
    // }
  }

  #[test]
  fn testok_hash_gpu_sys() {
    const DELTA_DEPTH: u8 = 1;
    // Input
    let depth = 4; // done with depth 13 (but too long to be tested systematically
                   // Computations
    let layer = nested::get(depth + DELTA_DEPTH);
    let mut d0h = 0;
    for h in 0..layer.n_hash() {
      if (h >> ((depth + 1) << 1)) != d0h {
        d0h += 1;
        println!("New base cell: {}", d0h);
      }
      // println!("h: {}", &h);
      let (ra, de) = layer.center(h);
      let (sl, cl) = ra.sin_cos();
      let (sb, cb) = de.sin_cos();
      let (x, y, z) = (cb * cl, cb * sl, sb);
      let (nh, dx, dy) = hash_with_dxdy(depth, x as f32, y as f32, z as f32);
      assert_eq!((h >> (DELTA_DEPTH << 1)) as u32, nh);
      let prec = 1.0 / 512.0; //512.0; // HiPS images: 512 x 512 with level max = 14
      match h & 3 {
        0 => {
          assert!(
            (dx - 0.25).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dx - 0.25).abs()
          ); // precision should depends on depth!!
          assert!(
            (dy - 0.25).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dy - 0.25).abs()
          ); // 1e-3 => images 1000x1000 while HiPS uses 512x512
        }
        1 => {
          assert!(
            (dx - 0.75).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dx - 0.75).abs()
          );
          assert!(
            (dy - 0.25).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dy - 0.25).abs()
          );
        }
        2 => {
          assert!(
            (dx - 0.25).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dx - 0.25).abs()
          );
          assert!(
            (dy - 0.75).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dy - 0.75).abs()
          );
        }
        3 => {
          assert!(
            (dx - 0.75).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dx - 0.75).abs()
          );
          assert!(
            (dy - 0.75).abs() <= prec,
            "h: {}; ra: {}; dec: {}; dx: {}; dy: {}; prec: {}",
            h,
            ra.to_degrees(),
            de.to_degrees(),
            dx,
            dy,
            (dy - 0.75).abs()
          );
        }
        _ => unreachable!(),
      }
    }
  }

  // No a real test, used to check the projection code visually
  #[test]
  fn testok_proj_gpu() {
    println!("ra,dec,x,y,color,ra_unproj,dec_unproj");
    for ra in 0..=360 {
      for dec in 0..=180 {
        let ra = (ra as f32).to_radians();
        let dec = ((dec - 90) as f32).to_radians();
        let x = ra.cos() * dec.cos();
        let y = ra.sin() * dec.cos();
        let z = dec.sin();
        let (x_proj, y_proj) = proj_gpu(x, y, z);
        let (x, y, z) = unproj_gpu(x_proj, y_proj);
        let r = (x * x + y * y).sqrt();
        let mut ra_unproj = y.atan2(x).to_degrees();
        let dec_unproj = z.atan2(r).to_degrees();
        if ra_unproj < 0.0 {
          ra_unproj += 360.0;
        }

        println!(
          "{},{},{},{},{},{},{}",
          ra.to_degrees(),
          dec.to_degrees(),
          x_proj,
          y_proj,
          ra + (dec + 2.0) * 2.0 * PI,
          ra_unproj,
          dec_unproj
        );
      }
    }
  }
}
