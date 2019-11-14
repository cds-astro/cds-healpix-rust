use criterion::{Criterion, black_box, criterion_group, criterion_main, BenchmarkId};

use cdshealpix::nested::zordercurve::{ZOrderCurve};


use cdshealpix::nested::zordercurve::LARGE_ZOC_LUT;
fn bench_morton_lupt(n: u32) -> u64 {
  let mut t = 0u64;
  for i in (0..n).step_by(149999910) {
    for j in (0..n).step_by(14999910) {
      t ^= LARGE_ZOC_LUT.ij2h(i, j);
    }
  }
  t
}


use cdshealpix::nested::zordercurve::LARGE_ZOC_XOR;
fn bench_morton_xor(n: u32) -> u64 {
  let mut t = 0u64;
  for i in (0..n).step_by(149999910) {
    for j in (0..n).step_by(14999910) {
      t ^= LARGE_ZOC_XOR.ij2h(i, j);
    }
  }
  t
}

#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
use cdshealpix::nested::zordercurve::LARGE_ZOC_BMI;
#[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
fn bench_morton_bmi(n: u32) -> u64 {
  let mut t = 0u64;
  for i in (0..n).step_by(149999910) {
    for j in (0..n).step_by(14999910) {
      t ^= LARGE_ZOC_BMI.ij2h(i, j);
    }
  }
  t
}

fn bench_zordercurve(c: &mut Criterion) {
  let mut group = c.benchmark_group("ZOrderCurve");
  group.sample_size(10);
  
  let n = black_box(0xFFFFFFFF_u32);
  group.bench_with_input(BenchmarkId::new("LUPT", 1), &1, |b, pos| b.iter(||
    bench_morton_lupt(n)
  ));
  group.bench_with_input(BenchmarkId::new("XOR", 2), &2,  |b, pos| b.iter(||
    bench_morton_xor(n)
  ));
  #[cfg(all(any(target_arch = "x86", target_arch = "x86_64"), target_feature = "bmi2"))]
  group.bench_with_input(BenchmarkId::new("BMI", 3), &3,  |b, pos| b.iter(||
    bench_morton_bmi(n)
  ));
  group.finish();
}

criterion_group!(zordercurce_benches, bench_zordercurve);

criterion_main!(zordercurce_benches);