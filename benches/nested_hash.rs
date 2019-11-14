use criterion::{Criterion, black_box, criterion_group, criterion_main, BenchmarkId};
use rand::Rng;

use std::f64::consts::PI;

use cdshealpix::{
  HALF_PI,
  TWICE_PI,
  TRANSITION_LATITUDE, 
  ONE_OVER_TRANSITION_Z, 
  SQRT6, 
  PI_OVER_FOUR, 
  FOUR_OVER_PI,
  F64_SIGN_BIT_MASK,
  F64_BUT_SIGN_BIT_MASK
};
use cdshealpix::nested::{self, get_or_create};
use cdshealpix::nested::zordercurve::{ZOrderCurve, get_zoc};

pub fn hash_v1(depth: u8, lon: f64, lat: f64) -> u64 {
  get_or_create(depth).hash_v1(lon, lat)
}

pub fn hash_v2(depth: u8, lon: f64, lat: f64) -> u64 {
  /*let nside = (1u32 << depth);
  let nside_minus_1 = nside - 1;
  let time_half_nside = ((depth - 1) as i64) << 52; // WARNING DO NOT WORK WITH depth=0
  
  let (d0h, l_in_d0c, h_in_d0c) = d0h_lh_in_d0c(lon, lat);
  // Coords inside the base cell
  //  - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
  let i = f64::from_bits((time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64) as u32;
  let j = f64::from_bits((time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64) as u32;
  //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
  let i = if i == nside { nside_minus_1 } else { i };
  let j = if j == nside { nside_minus_1 } else { j };
  build_hash_from_parts(depth, d0h, i, j)*/
  get_or_create(depth).hash_v2(lon, lat)
}

pub fn hash_v3(depth: u8, lon: f64, lat: f64) -> u64 {
  let nside = (1u32 << depth);
  let nside_minus_1 = nside - 1;
  let time_half_nside = ((depth - 1) as i64) << 52; // WARNING DO NOT WORK WITH depth=0
  
  let (d0h, l_in_d0c, h_in_d0c) = d0h_lh_in_d0c_v2(lon, lat);
  // Coords inside the base cell
  //  - ok to cast on u32 since small negative values due to numerical inaccuracies (like -1e-15), are rounded to 0
  let i = f64::from_bits((time_half_nside + (h_in_d0c + l_in_d0c).to_bits() as i64) as u64) as u32;
  let j = f64::from_bits((time_half_nside + (h_in_d0c - l_in_d0c).to_bits() as i64) as u64) as u32;
  //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
  let i = if i == nside { nside_minus_1 } else { i };
  let j = if j == nside { nside_minus_1 } else { j };
  build_hash_from_parts(depth, d0h, i, j)
}

#[inline]
fn build_hash_from_parts(depth: u8, d0h: u8, i: u32, j: u32) -> u64 {
  build_hash(depth,(d0h as u64) << (depth << 1_u8), i, j)
}

#[inline]
fn build_hash(depth: u8, d0h_bits: u64, i: u32, j: u32) -> u64 {
  d0h_bits | get_zoc(depth).ij2h(i, j)
}

fn d0h_lh_in_d0c(lon: f64, lat: f64) -> (u8, f64, f64) {
  let (x_pm1, q) = xpm1_and_q(lon);
  if lat > TRANSITION_LATITUDE {
    // North polar cap, Collignon projection.
    // - set the origin to (PI/4, 0)
    let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 + PI_OVER_FOUR).cos();
    let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
    let d0h = q;
    (d0h, x_proj, y_proj)
  } else if lat < -TRANSITION_LATITUDE {
    // South polar cap, Collignon projection
    // - set the origin to (PI/4, -PI/2)
    let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 - PI_OVER_FOUR).cos(); // cos(-x) = cos(x)
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
    let y_pm1 = lat.sin() * ONE_OVER_TRANSITION_Z;
    // Inequalities have been carefully chosen so that S->E and S->W axis are part of the cell,
    // and not E->N and W->N
    // Branch free version
    // |\2/|
    // .3X1.
    // |/0\|
    let q01 = (x_pm1 >   y_pm1) as u8; /* 0/1 */ debug_assert!(q01 == 0 || q01 == 1);
    let q12 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */ debug_assert!(q12 == 0 || q12 == 1);
    let q1 = q01 & q12; /* = 1 if q1, 0 else */      debug_assert!( q1 == 0 ||  q1 == 1);
    let q013 = q01 + (1 - q12); // = q01 + q03; /* 0/1 + 1\0 +  */
    // x: x_pm1 + 1 if q3 | x_pm1 - 1 if q1 | x_pm1 if q0 or q2
    let x_proj = x_pm1 - ((q01 + q12) as i8 - 1) as f64;
    // y: y_pm1 + 0 if q2 | y_pm1 + 1 if q1 or q3 | y_pm1 + 2 if q0 
    let y_proj = y_pm1 + q013 as f64;
    // d0h: +8 if q0 | +4 if q3 | +5 if q1
    let d0h = (q013 << 2) + ((q + q1) & 3);
    (d0h, x_proj, y_proj)
  }
}

fn d0h_lh_in_d0c_v2(lon: f64, lat: f64) -> (u8, f64, f64) {
  let (x_pm1, q) = xpm1_and_q(lon);
  if lat > TRANSITION_LATITUDE {
    // North polar cap, Collignon projection.
    // - set the origin to (PI/4, 0)
    let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 + PI_OVER_FOUR).cos();
    let (x_proj, y_proj) = (x_pm1 * sqrt_3_one_min_z, 2.0 - sqrt_3_one_min_z);
    let d0h = q;
    (d0h, x_proj, y_proj)
  } else if lat < -TRANSITION_LATITUDE {
    // South polar cap, Collignon projection
    // - set the origin to (PI/4, -PI/2)
    let sqrt_3_one_min_z = SQRT6 * (lat / 2.0 - PI_OVER_FOUR).cos(); // cos(-x) = cos(x)
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
    let y_pm1 = lat.sin() * ONE_OVER_TRANSITION_Z;
    // Inequalities have been carefully chosen so that S->E and S->W axis are part of the cell,
    // and not E->N and W->N
    // Version with branch
    // |\3/|
    // .2X1.
    // |/0\|
    let q13 = (x_pm1 >= -y_pm1) as u8; /* 0\1 */ debug_assert!(q13 == 0 || q13 == 1);
    let q23 = (x_pm1 <=  y_pm1) as u8; /* 1/0 */ debug_assert!(q23 == 0 || q23 == 1);
    match q13 | (q23 << 1) {
      0 => ( q         , x_pm1      , y_pm1 + 2.0),
      1 => ((q + 5) & 7, x_pm1 - 1.0, y_pm1 + 1.0), // (q + 5) & 7 <=> (q + 1) | 4
      2 => ( q + 4     , x_pm1 + 1.0, y_pm1 + 1.0),
      3 => ( q + 8     , x_pm1      , y_pm1),
      _ => unreachable!(),
    }
  }
}

/// Transform the input longitude, in radians, in a value `x` in `[-1, 1[` plus a quarter in `[0, 3]`,
/// such that `lon = (x + 1) * PI / 4 + q * PI / 2`. 
#[inline]
fn xpm1_and_q(lon: f64) -> (f64, u8) {
  let lon_bits = lon.to_bits();
  let lon_abs = f64::from_bits(lon_bits & F64_BUT_SIGN_BIT_MASK);
  let lon_sign = lon_bits & F64_SIGN_BIT_MASK;
  let x = lon_abs * FOUR_OVER_PI;
  let q = (x as u8 | 1_u8) & 7_u8;    debug_assert!(0 <= q && q < 8);
  // Remark: to avoid the branch, we could have copied lon_sign on x - q, 
  //         but I so far lack of idea to deal with q efficiently.
  //         And we are not supposed to have negative longitudes in ICRS 
  //         (the most used reference system in astronomy).
  if lon_sign == 0 { // => lon >= 0
    (x - (q as f64), q >> 1)
  } else { // case lon < 0 should be rare => few risks of branch miss-prediction
    // Since q in [0, 3]: 3 - (q >> 1)) <=> 3 & !(q >> 1)
    // WARNING: BE SURE TO HANDLE THIS CORRECTLY IN THE REMAINING OF THE CODE!
    //  - Case lon =  3/4 pi = 270 deg => x = -1, q=3
    //  - Case lon = -1/2 pi = -90 deg => x =  1, q=2
    ((q as f64) - x, 3 - (q >> 1))
  }
}


fn gen_rand_lonlat(n: usize) -> Vec<(f64, f64)> {
  // let v: Vec<f64, f64> = Vec::with_capacity(n);
  let mut rng = rand::thread_rng();
  (0..n).into_iter()
    .map(|_| (rng.gen::<f64>() * TWICE_PI, rng.gen::<f64>() * PI - HALF_PI))
    .collect()
}

pub fn benchmark_hash_v1(depth: u8, positions: &Vec<(f64, f64)>) -> u64 {
  let mut sum: u64 = 0;
  for (lon, lat) in positions {
    sum |= hash_v1(depth, *lon, *lat);
  }
  sum
}

pub fn benchmark_hash_v2(depth: u8, positions: &Vec<(f64, f64)>) -> u64 {
  let mut sum: u64 = 0;
  for (lon, lat) in positions {
    sum |= hash_v2(depth, *lon, *lat);
  }
  sum
}

pub fn benchmark_hash_v3(depth: u8, positions: &Vec<(f64, f64)>) -> u64 {
  let mut sum: u64 = 0;
  for (lon, lat) in positions {
    sum |= hash_v3(depth, *lon, *lat);
  }
  sum
}

fn bench_time_pow2(c: &mut Criterion) {
  let mut group = c.benchmark_group("Time power of 2");

  for depth in (1..29) {
    let time_half_nside = ((depth - 1) as i64) << 52;
    let half_nside = (1u32 << depth) as f64 / 2.0;
    let one_over_half_nside = 1.0 / half_nside;
    let val: f64 = black_box(PI);
    group.bench_with_input(BenchmarkId::new("x2^n v1", 1), &depth, |b, pos| b.iter(||
      f64::from_bits((time_half_nside + val.to_bits() as i64) as u64)
    ));
    group.bench_with_input(BenchmarkId::new("x2^n v2", 2), &depth,  |b, pos| b.iter(|| 
      val / half_nside
    ));
    group.bench_with_input(BenchmarkId::new("x2^n v3", 3), &depth,  |b, pos| b.iter(||
      val * one_over_half_nside
    ));
  }
  group.finish();
}

fn bench_hash(c: &mut Criterion) {
  let mut group = c.benchmark_group("Hash");
  group.sample_size(10);
  
  let depth = 16;
  let positions = gen_rand_lonlat(black_box(1000000));
  group.bench_with_input(BenchmarkId::new("Hash v1", 1), &1,
                         |b, pos| b.iter(|| benchmark_hash_v1(depth, &positions)));
  group.bench_with_input(BenchmarkId::new("Hash v2", 2), &2,
                         |b, pos| b.iter(|| benchmark_hash_v2(depth, &positions)));
  group.bench_with_input(BenchmarkId::new("Hash v3", 3), &3,
                         |b, pos| b.iter(|| benchmark_hash_v3(depth, &positions)));
  group.finish();
}

criterion_group!(hash_benches, bench_hash/*, bench_time_pow2*/);
// criterion_group!(pow2_time_benches, bench_time_pow2);

criterion_main!(hash_benches);



