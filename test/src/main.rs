
use base64::{encode, decode, DecodeError};

use std::time::{Duration, Instant};
use std::cmp::max;

use cdshealpix::nested::bmoc::*;
use cdshealpix::nested::moc::{
  self,
  HpxCell, OwnedMOC, LazyMOCIter, MOCIterator,
  OwnedRangeMOC, RangeMOCIterator
};
use cdshealpix::nested::moc::op::{not, not_unchecked, and, and_unchecked, or, or_unchecked};
use cdshealpix::nested::moc::compressed::{compress_unchecked, uncompress};

use cdshealpix_test::*;

fn test_sdss_not() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let sdss_not = load_sdss_not().unwrap().to_bmoc();
  
  // Test classic version
  let now = Instant::now();
  let not_sdss = sdss.not();
  println!("SDSS 'not' operation done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), sdss.entries.len(), not_sdss.entries.len());
  assert!(sdss_not.equals(&not_sdss));
  
  // Test iterator version
  let now = Instant::now();
  let depth_max = sdss.get_depth_max();
  let mut builder = BMOCBuilderUnsafe::new(depth_max, sdss.size());
  for HpxCell {depth, hash} in not_unchecked(
    LazyMOCIter::new(depth_max, sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)))) {
    builder.push(depth, hash, true);
  }
  let not_sdss = builder.to_bmoc();
  println!("SDSS 'not' operation from iterators done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), sdss.entries.len(), not_sdss.entries.len());
  assert!(sdss_not.equals(&not_sdss));
}

fn test_glimpse_not() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let not_glimpse = load_not_glimpse().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let glimpse_not = glimpse.not();
  println!("GLIMPSE 'not' operation done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), glimpse.entries.len(), not_glimpse.entries.len());
  assert!(not_glimpse.equals(&glimpse_not));

  // Test iterator version
  let now = Instant::now();
  let depth_max = glimpse.get_depth_max();
  let mut builder = BMOCBuilderUnsafe::new(depth_max, glimpse.size());
  for HpxCell {depth, hash} in not_unchecked(
    LazyMOCIter::new(depth_max, glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)))) {
    builder.push(depth, hash, true);
  }
  let glimpse_not = builder.to_bmoc();
  println!("GLIMPSE 'not' operation from iterators done in {} ms. In / Out sizes: {} / {}",
           now.elapsed().as_millis(), glimpse.entries.len(), not_glimpse.entries.len());
  assert!(not_glimpse.equals(&glimpse_not));
  
}

fn test_and() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_and_g = load_s_and_g().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let sandg = sdss.and(&glimpse);
  println!("SDSS/GLIMPSE 'and' operation done in {} ms. In / Out sizes: {} and {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sandg.entries.len());
  assert!(s_and_g.equals(&sandg));
  
  // Test iterator version
  let now = Instant::now();
  let mut builder = BMOCBuilderUnsafe::new(glimpse.get_depth_max().max(sdss.get_depth_max()), sdss.size());
  for HpxCell {depth, hash} in and_unchecked(
    LazyMOCIter::new(glimpse.get_depth_max(), glimpse.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) })),
    LazyMOCIter::new(sdss.get_depth_max(), sdss.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) }))
  ) {
    builder.push(depth, hash, true);
  }
  let sandg = builder.to_bmoc();
  println!("SDSS/GLIMPSE 'and' operation from iterators done in {} ms. In / Out sizes: {} and {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sandg.entries.len());
  assert!(s_and_g.equals(&sandg));
}

fn test_or() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_or_g = load_s_or_g().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let sorg = &sdss.or(&glimpse);
  println!("SDSS/GLIMPSE 'or' operation done in {} ms. In / Out sizes: {} or {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sorg.entries.len());
  s_or_g.assert_equals(&sorg);
  assert!(s_or_g.equals(&sorg));

  // Test iterator version
  let now = Instant::now();
  let mut builder = BMOCBuilderUnsafe::new(glimpse.get_depth_max().max(sdss.get_depth_max()), sdss.size());
  /*Code used to debug
  let mut i = 0;
  for (HpxCell {depth, hash}, HpxCell {depth: depth2, hash: hash2}) in or(
    glimpse.get_depth_max(),
    glimpse.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) }),
    sdss.get_depth_max(),
    sdss.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) })
  ).zip(s_or_g.into_iter().map(|c| HpxCell::new(c.depth, c.hash))) {
    if depth != depth2 || hash != hash2 {
      println!("KO d: {}, h: {}, d2: {}, h2: {}, i: {}", depth, hash, depth2, hash2, i);
      break;
    }/* else {
      println!("OK d: {}, h: {}, d2: {}, h2: {}, i: {}", depth, hash, depth2, hash2, i);
    }*/
    i += 1;
  }*/
  for HpxCell {depth, hash} in or_unchecked(
    LazyMOCIter::new(glimpse.get_depth_max(), glimpse.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) })),
    LazyMOCIter::new(sdss.get_depth_max(), sdss.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) }))
  ) {
    builder.push(depth, hash, true);
  }
  let sorg = builder.to_bmoc();
  println!("SDSS/GLIMPSE 'or' operation from iterators done in {} ms. In / Out sizes: {} or {} => {}",
           now.elapsed().as_millis(), sdss.entries.len(), glimpse.entries.len(), sorg.entries.len());
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
 // println!("b64 MOC: {}", b64);
  
  // Iteartor version
  let now = Instant::now();
  let compressed_2: Vec<u8> = compress_unchecked(
    LazyMOCIter::new(glimpse.get_depth_max(), glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)))
  ).collect();
  println!("GLIMPSE moc compression from iterator done in {} ms. Out size: {} bytes",
           now.elapsed().as_millis(), compressed_2.len());
  let now = Instant::now();
  let b64_2 = encode(&compressed_2[..]);
  println!("GLIMPSE b64 done in {} ms. Out size: {} characters",
           now.elapsed().as_millis(), b64_2.len());
  assert_eq!(b64, b64_2);

  let now = Instant::now();
  let decompressed_2: Vec<HpxCell<u64>> = uncompress(&compressed_2[..]).collect();
  println!("GLIMPSE decompression for iterator done in {} ms. Out moc size: {} cells",
           now.elapsed().as_millis(), decompressed_2.len());
  assert_eq!(glimpse.iter().count(), decompressed_2.len());
  for (HpxCell {depth, hash}, HpxCell {depth: depth2, hash: hash2}) in
    glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).zip(uncompress::<u64>(&compressed_2[..])) {
    assert_eq!(depth, depth2);
    assert_eq!(hash,hash2);
  }
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

  // Iteartor version
  let now = Instant::now();
  let compressed_2: Vec<u8> = compress_unchecked(
    LazyMOCIter::new(sdss.get_depth_max(), sdss.into_iter().map(|c| { HpxCell::new(c.depth, c.hash) }))
  ).collect();
  println!("SDSS moc compression from iterator done in {} ms. Out size: {} bytes",
           now.elapsed().as_millis(), compressed_2.len());
  let now = Instant::now();
  let b64_2 = encode(&compressed_2[..]);
  println!("SDSS b64 done in {} ms. Out size: {} characters",
           now.elapsed().as_millis(), b64_2.len());
  assert_eq!(b64, b64_2);

  let now = Instant::now();
  let decompressed_2: Vec<HpxCell<u64>> = uncompress(&compressed_2[..]).collect();
  println!("SDSS decompression for iterator done in {} ms. Out moc size: {} cells",
           now.elapsed().as_millis(), decompressed_2.len());
  assert_eq!(sdss.iter().count(), decompressed_2.len());
  for (HpxCell {depth, hash}, HpxCell {depth: depth2, hash: hash2}) in
    sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).zip(uncompress::<u64>(&compressed_2[..])) {
    assert_eq!(depth, depth2);
    assert_eq!(hash,hash2);
  }

}

fn test_and_compressed_iterators() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_and_g = load_s_and_g().unwrap().to_bmoc();
  
  let compressed_glimpse: Vec<u8> = compress_unchecked(
    LazyMOCIter::new(glimpse.get_depth_max(), glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)))
  ).collect();

  let compressed_sdss: Vec<u8> = compress_unchecked(
    LazyMOCIter::new(sdss.get_depth_max(), sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)))
  ).collect();

  // Decompress both MOCs while performing the AND operation and compressing on the fly the result
  let now = Instant::now();
  let compressed_and: Vec<u8> = compress_unchecked(
    and_unchecked(
      LazyMOCIter::new(glimpse.get_depth_max(), uncompress::<u64>(&compressed_glimpse[..])),
      LazyMOCIter::new(sdss.get_depth_max(), uncompress::<u64>(&compressed_sdss[..]))
    )
  ).collect();
  println!("SDSS/GLIMPSE 'and' operation from compressed MOCs done in {} ms. \
            Compressed In / Compressed Out sizes: {} and {} => {}",
           now.elapsed().as_millis(), compressed_glimpse.len(), compressed_sdss.len(), compressed_and.len());
  
}

fn test_ranges() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc2 = LazyMOCIter::new(moc.get_depth_max(), moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)));
  // println!("N ranges: {} {}", moc1.len(), &moc2.to_range_iter().count());
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = moc2.to_range_iter().to_moc_iter().collect();
  println!("To range/to cells/collect of size: {} in {} ms", moc1.len(), now.elapsed().as_millis());
  assert_eq!(moc1.len(), moc3.len());
  for (HpxCell {depth, hash}, HpxCell {depth: depth2, hash: hash2}) in
    moc1.into_iter().zip(&mut moc3.into_iter()) {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_ranges_2() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc2 = OwnedRangeMOC::from_it(
    moc.get_depth_max(),
    LazyMOCIter::new(
      moc.get_depth_max(),
      moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash))
    ).to_range_iter()
  );
  // println!("N ranges: {} {}", moc1.len(), &moc2.to_range_iter().count());
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = moc2.into_range_moc_iter().to_moc_iter().collect();
  println!("To range + collect of size: {} done in {} ms", moc1.len(), now.elapsed().as_millis());
  assert_eq!(moc1.len(), moc3.len());
  for (HpxCell {depth, hash}, HpxCell {depth: depth2, hash: hash2}) in
    moc1.into_iter().zip(&mut moc3.into_iter()) {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_expand() {
  let moc = load_sdss().unwrap().to_bmoc();
  // let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc1 = OwnedMOC::from_it_unchecked(moc.get_depth_max(), moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)));
  let len = moc1.len();
  let now = Instant::now();
  let moc2: Vec<HpxCell<u64>> = moc1.expand().collect();
  println!("Expand + collect of sizes: {} => {} at depth {} done in {} ms", 
           len, moc2.len(), moc.get_depth_max(), now.elapsed().as_millis());
}


fn test_expand_v2() {
  let moc = load_sdss().unwrap().to_bmoc();
  // let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc1 = OwnedMOC::from_it_unchecked(moc.get_depth_max(), moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)));
  let moc2 = OwnedMOC::from_it_unchecked(moc.get_depth_max(), moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)));
  let len = moc1.len();
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = or(moc1.into_moc_iter(), moc2.external_border()).collect();
  println!("Expand + collect of sizes: {} => {} at depth {} done in {} ms",
           len, moc3.len(), moc.get_depth_max(), now.elapsed().as_millis());
}

fn test_expand_ranges() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1 = OwnedMOC::<u64>::from_it_unchecked(moc.get_depth_max(), moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)));
  let len = moc1.len();
  let ranges = OwnedRangeMOC::from_it(moc.get_depth_max(), moc1.into_moc_iter().to_range_iter());
  let now = Instant::now();
  let moc2: Vec<HpxCell<u64>> = ranges.expand().collect();
  println!("Expand from ranges + collect of sizes: {} => {} at depth {} done in {} ms",
           len, moc2.len(), moc.get_depth_max(), now.elapsed().as_millis());
}

pub fn main() {
  test_sdss_not();
  test_glimpse_not();
  test_and();
  test_or();
  test_xor();

  test_compress_glimpse();
  test_compress_sdss();

  test_sdss_build_from_ordered_input();

 test_and_compressed_iterators();

  test_ranges();
  test_ranges_2();

  test_expand();
  test_expand_v2();
  test_expand_ranges();
}