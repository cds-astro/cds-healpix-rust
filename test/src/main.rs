
use std::time::{Duration, Instant};

use cdshealpix::nested::bmoc::*;

use cdshealpix_test::*;

fn test_sdss_not() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let sdss_not = load_sdss_not().unwrap().to_bmoc();

  let now = Instant::now();
  let not_sdss = sdss.not();
  println!("SDSS 'not' operation done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), sdss.entries.len(), not_sdss.entries.len());

  // sdss_not.assert_equals(&not_sdss);
  assert!(sdss_not.equals(&not_sdss));
}
fn test_glimpse_not() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let not_glimpse = load_not_glimpse().unwrap().to_bmoc();

  let now = Instant::now();
  let glimpse_not = glimpse.not();
  println!("GLIMPSE 'not' operation done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), glimpse.entries.len(), not_glimpse.entries.len());

  assert!(not_glimpse.equals(&glimpse_not));
}
fn test_and() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();

  let s_and_g = load_s_and_g().unwrap().to_bmoc();

  let now = Instant::now();
  let sandg = sdss.and(&glimpse);
  println!("SDSS/GLIMPSE 'and' operation done in {} ms. In / Out sizes: {} and {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sandg.entries.len());

  assert!(s_and_g.equals(&sandg));
}
fn test_or() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_or_g = load_s_or_g().unwrap().to_bmoc();

  let now = Instant::now();
  let sorg = &sdss.or(&glimpse);
  println!("SDSS/GLIMPSE 'or' operation done in {} ms. In / Out sizes: {} or {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sorg.entries.len());

  s_or_g.assert_equals(&sorg);
  assert!(s_or_g.equals(&sorg));
}
fn test_xor() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_xor_g = load_s_xor_g().unwrap().to_bmoc();

  let now = Instant::now();
  let sxorg = sdss.xor(&glimpse);
  println!("SDSS/GLIMPSE 'xor' operation done in {} ms. In / Out sizes: {} xor {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sxorg.entries.len());

  s_xor_g.assert_equals(&sxorg);
  assert!(s_xor_g.equals(&sdss.xor(&glimpse)));
}
fn test_sdss_build_from_ordered_input() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let deep_size = sdss.deep_size();

  let mut builder = BMOCBuilderFixedDepth::new(sdss.get_depth_max(), true);

  let now = Instant::now();
  for h in sdss.flat_iter() {
    builder.push(h);
  }
  let sdss2 = builder.to_bmoc().unwrap();
  println!("SDSS 'build' operation (from SDSS flat iterator, deep size: {}) done in {} ms", deep_size, now.elapsed().as_millis());


  sdss.assert_equals(&sdss2);
  assert!(sdss.equals(&sdss2));
}

fn test_compress_glimpse() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let now = Instant::now();
  let compressed = glimpse.compress_lossy();
  println!("GLIMPSE moc compression done in {} ms. Out size: {} bytes",
           now.elapsed().as_millis(), compressed.byte_size());
  let now = Instant::now();
  let b64 = compressed.to_b64();
  println!("GLIMPSE b64 done in {} ms. Out size: {} characters",
           now.elapsed().as_millis(), b64.len());
  let now = Instant::now();
  let decompressed = compressed.self_decompress();
  println!("GLIMPSE decompression done in {} ms. Out moc size: {} cells",
           now.elapsed().as_millis(), decompressed.size());

  assert_eq!(decompressed.compress_lossy().to_b64(), b64);
  
  println!("b64 MOC: {}", b64);
}

fn test_compress_sdss() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let now = Instant::now();
  let compressed = sdss.compress_lossy();
  println!("SDSS moc compression done in {} ms. Out size: {} bytes",
           now.elapsed().as_millis(), compressed.byte_size());
  let now = Instant::now();
  let b64 = compressed.to_b64();
  println!("SDSS b64 done in {} ms. Out size: {} characters",
           now.elapsed().as_millis(), b64.len());
  let now = Instant::now();
  let decompressed = compressed.self_decompress();
  println!("SDSS decompression done in {} ms. Out moc size: {} cells",
           now.elapsed().as_millis(), decompressed.size());

  assert_eq!(decompressed.compress_lossy().to_b64(), b64);

  // println!("b64 MOC: {}", b64);

}

pub fn main() {
  test_sdss_not();
  test_glimpse_not();
  test_and();
  test_or();
  test_xor();
  test_sdss_build_from_ordered_input();

  test_compress_sdss();
  test_compress_glimpse();

}