use base64::{decode, encode, DecodeError};

use std::cmp::max;
use std::time::{Duration, Instant};

use cdshealpix::nested::bmoc::*;
use cdshealpix::nested::moc::compressed::{compress_unchecked, uncompress};
use cdshealpix::nested::moc::op::{and, and_unchecked, not, not_unchecked, or, or_unchecked};
use cdshealpix::nested::moc::{
  self, HpxCell, LazyMOCIter, MOCIterator, OwnedMOC, OwnedRangeMOC, RangeMOCIterator,
};

use cdshealpix_test::*;

fn test_sdss_not() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let sdss_not = load_sdss_not().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let not_sdss = sdss.not();
  println!(
    "SDSS 'not' operation done in {} ms. In / Out sizes: {} / {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    not_sdss.entries.len()
  );
  assert!(sdss_not.equals(&not_sdss));

  // Test iterator version
  let now = Instant::now();
  let depth_max = sdss.get_depth_max();
  let mut builder = BMOCBuilderUnsafe::new(depth_max, sdss.size());
  for HpxCell { depth, hash } in not_unchecked(LazyMOCIter::new(
    depth_max,
    sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  )) {
    builder.push(depth, hash, true);
  }
  let not_sdss = builder.to_bmoc();
  println!(
    "SDSS 'not' operation from iterators done in {} ms. In / Out sizes: {} / {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    not_sdss.entries.len()
  );
  assert!(sdss_not.equals(&not_sdss));
}

fn test_glimpse_not() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let not_glimpse = load_not_glimpse().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let glimpse_not = glimpse.not();
  println!(
    "GLIMPSE 'not' operation done in {} ms. In / Out sizes: {} / {}",
    now.elapsed().as_millis(),
    glimpse.entries.len(),
    not_glimpse.entries.len()
  );
  assert!(not_glimpse.equals(&glimpse_not));

  // Test iterator version
  let now = Instant::now();
  let depth_max = glimpse.get_depth_max();
  let mut builder = BMOCBuilderUnsafe::new(depth_max, glimpse.size());
  for HpxCell { depth, hash } in not_unchecked(LazyMOCIter::new(
    depth_max,
    glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  )) {
    builder.push(depth, hash, true);
  }
  let glimpse_not = builder.to_bmoc();
  println!(
    "GLIMPSE 'not' operation from iterators done in {} ms. In / Out sizes: {} / {}",
    now.elapsed().as_millis(),
    glimpse.entries.len(),
    not_glimpse.entries.len()
  );
  assert!(not_glimpse.equals(&glimpse_not));
}

fn test_and() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_and_g = load_s_and_g().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let sandg = sdss.and(&glimpse);
  println!(
    "SDSS/GLIMPSE 'and' operation done in {} ms. In / Out sizes: {} and {} => {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    glimpse.entries.len(),
    sandg.entries.len()
  );
  assert!(s_and_g.equals(&sandg));

  // Test iterator version
  let now = Instant::now();
  let mut builder = BMOCBuilderUnsafe::new(
    glimpse.get_depth_max().max(sdss.get_depth_max()),
    sdss.size(),
  );
  for HpxCell { depth, hash } in and_unchecked(
    LazyMOCIter::new(
      glimpse.get_depth_max(),
      glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
    ),
    LazyMOCIter::new(
      sdss.get_depth_max(),
      sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
    ),
  ) {
    builder.push(depth, hash, true);
  }
  let sandg = builder.to_bmoc();
  println!(
    "SDSS/GLIMPSE 'and' operation from iterators done in {} ms. In / Out sizes: {} and {} => {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    glimpse.entries.len(),
    sandg.entries.len()
  );
  assert!(s_and_g.equals(&sandg));
}

fn test_or() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_or_g = load_s_or_g().unwrap().to_bmoc();

  // Test classic version
  let now = Instant::now();
  let sorg = &sdss.or(&glimpse);
  println!(
    "SDSS/GLIMPSE 'or' operation done in {} ms. In / Out sizes: {} or {} => {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    glimpse.entries.len(),
    sorg.entries.len()
  );
  s_or_g.assert_equals(&sorg);
  assert!(s_or_g.equals(&sorg));

  // Test iterator version
  let now = Instant::now();
  let mut builder = BMOCBuilderUnsafe::new(
    glimpse.get_depth_max().max(sdss.get_depth_max()),
    sdss.size(),
  );
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
  for HpxCell { depth, hash } in or_unchecked(
    LazyMOCIter::new(
      glimpse.get_depth_max(),
      glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
    ),
    LazyMOCIter::new(
      sdss.get_depth_max(),
      sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
    ),
  ) {
    builder.push(depth, hash, true);
  }
  let sorg = builder.to_bmoc();
  println!(
    "SDSS/GLIMPSE 'or' operation from iterators done in {} ms. In / Out sizes: {} or {} => {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    glimpse.entries.len(),
    sorg.entries.len()
  );
  assert!(s_or_g.equals(&sorg));
}

fn test_xor() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_xor_g = load_s_xor_g().unwrap().to_bmoc();

  let now = Instant::now();
  let sxorg = sdss.xor(&glimpse);
  println!(
    "SDSS/GLIMPSE 'xor' operation done in {} ms. In / Out sizes: {} xor {} => {}",
    now.elapsed().as_millis(),
    sdss.entries.len(),
    glimpse.entries.len(),
    sxorg.entries.len()
  );

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
  println!(
    "SDSS 'build' operation (from SDSS flat iterator, deep size: {}) done in {} ms",
    deep_size,
    now.elapsed().as_millis()
  );

  sdss.assert_equals(&sdss2);
  assert!(sdss.equals(&sdss2));
}

fn test_compress_glimpse() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let now = Instant::now();
  let compressed = glimpse.compress_lossy();
  println!(
    "GLIMPSE moc compression done in {} ms. Out size: {} bytes",
    now.elapsed().as_millis(),
    compressed.byte_size()
  );
  let now = Instant::now();
  let b64 = compressed.to_b64();
  println!(
    "GLIMPSE b64 done in {} ms. Out size: {} characters",
    now.elapsed().as_millis(),
    b64.len()
  );
  let now = Instant::now();
  let decompressed = compressed.self_decompress();
  println!(
    "GLIMPSE decompression done in {} ms. Out moc size: {} cells",
    now.elapsed().as_millis(),
    decompressed.size()
  );
  assert_eq!(decompressed.compress_lossy().to_b64(), b64);
  // println!("b64 MOC: {}", b64);

  // Iteartor version
  let now = Instant::now();
  let compressed_2: Vec<u8> = compress_unchecked(LazyMOCIter::new(
    glimpse.get_depth_max(),
    glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  ))
  .collect();
  println!(
    "GLIMPSE moc compression from iterator done in {} ms. Out size: {} bytes",
    now.elapsed().as_millis(),
    compressed_2.len()
  );
  let now = Instant::now();
  let b64_2 = encode(&compressed_2[..]);
  println!(
    "GLIMPSE b64 done in {} ms. Out size: {} characters",
    now.elapsed().as_millis(),
    b64_2.len()
  );
  assert_eq!(b64, b64_2);

  let now = Instant::now();
  let decompressed_2: Vec<HpxCell<u64>> = uncompress(&compressed_2[..]).collect();
  println!(
    "GLIMPSE decompression for iterator done in {} ms. Out moc size: {} cells",
    now.elapsed().as_millis(),
    decompressed_2.len()
  );
  assert_eq!(glimpse.iter().count(), decompressed_2.len());
  for (
    HpxCell { depth, hash },
    HpxCell {
      depth: depth2,
      hash: hash2,
    },
  ) in glimpse
    .into_iter()
    .map(|c| HpxCell::new(c.depth, c.hash))
    .zip(uncompress::<u64>(&compressed_2[..]))
  {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_compress_sdss() {
  let sdss = load_sdss().unwrap().to_bmoc();
  let now = Instant::now();
  let compressed = sdss.compress_lossy();
  println!(
    "SDSS moc compression done in {} ms. Out size: {} bytes",
    now.elapsed().as_millis(),
    compressed.byte_size()
  );
  let now = Instant::now();
  let b64 = compressed.to_b64();
  println!(
    "SDSS b64 done in {} ms. Out size: {} characters",
    now.elapsed().as_millis(),
    b64.len()
  );
  let now = Instant::now();
  let decompressed = compressed.self_decompress();
  println!(
    "SDSS decompression done in {} ms. Out moc size: {} cells",
    now.elapsed().as_millis(),
    decompressed.size()
  );

  assert_eq!(decompressed.compress_lossy().to_b64(), b64);
  // println!("b64 MOC: {}", b64);

  // Iteartor version
  let now = Instant::now();
  let compressed_2: Vec<u8> = compress_unchecked(LazyMOCIter::new(
    sdss.get_depth_max(),
    sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  ))
  .collect();
  println!(
    "SDSS moc compression from iterator done in {} ms. Out size: {} bytes",
    now.elapsed().as_millis(),
    compressed_2.len()
  );
  let now = Instant::now();
  let b64_2 = encode(&compressed_2[..]);
  println!(
    "SDSS b64 done in {} ms. Out size: {} characters",
    now.elapsed().as_millis(),
    b64_2.len()
  );
  assert_eq!(b64, b64_2);

  let now = Instant::now();
  let decompressed_2: Vec<HpxCell<u64>> = uncompress(&compressed_2[..]).collect();
  println!(
    "SDSS decompression for iterator done in {} ms. Out moc size: {} cells",
    now.elapsed().as_millis(),
    decompressed_2.len()
  );
  assert_eq!(sdss.iter().count(), decompressed_2.len());
  for (
    HpxCell { depth, hash },
    HpxCell {
      depth: depth2,
      hash: hash2,
    },
  ) in sdss
    .into_iter()
    .map(|c| HpxCell::new(c.depth, c.hash))
    .zip(uncompress::<u64>(&compressed_2[..]))
  {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_and_compressed_iterators() {
  let glimpse = load_glimpse().unwrap().to_bmoc();
  let sdss = load_sdss().unwrap().to_bmoc();
  let s_and_g = load_s_and_g().unwrap().to_bmoc();

  let compressed_glimpse: Vec<u8> = compress_unchecked(LazyMOCIter::new(
    glimpse.get_depth_max(),
    glimpse.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  ))
  .collect();

  let compressed_sdss: Vec<u8> = compress_unchecked(LazyMOCIter::new(
    sdss.get_depth_max(),
    sdss.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  ))
  .collect();

  // Decompress both MOCs while performing the AND operation and compressing on the fly the result
  let now = Instant::now();
  let compressed_and: Vec<u8> = compress_unchecked(and_unchecked(
    LazyMOCIter::new(
      glimpse.get_depth_max(),
      uncompress::<u64>(&compressed_glimpse[..]),
    ),
    LazyMOCIter::new(
      sdss.get_depth_max(),
      uncompress::<u64>(&compressed_sdss[..]),
    ),
  ))
  .collect();
  println!(
    "SDSS/GLIMPSE 'and' operation from compressed MOCs done in {} ms. \
            Compressed In / Compressed Out sizes: {} and {} => {}",
    now.elapsed().as_millis(),
    compressed_glimpse.len(),
    compressed_sdss.len(),
    compressed_and.len()
  );
}

fn test_ranges() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1: Vec<HpxCell<u64>> = moc
    .into_iter()
    .map(|c| HpxCell::new(c.depth, c.hash))
    .collect();
  let moc2 = LazyMOCIter::new(
    moc.get_depth_max(),
    moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  );
  // println!("N ranges: {} {}", moc1.len(), &moc2.to_range_iter().count());
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = moc2.to_range_iter().to_moc_iter().collect();
  println!(
    "To range/to cells/collect of size: {} in {} ms",
    moc1.len(),
    now.elapsed().as_millis()
  );
  assert_eq!(moc1.len(), moc3.len());
  for (
    HpxCell { depth, hash },
    HpxCell {
      depth: depth2,
      hash: hash2,
    },
  ) in moc1.into_iter().zip(&mut moc3.into_iter())
  {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_ranges_2() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1: Vec<HpxCell<u64>> = moc
    .into_iter()
    .map(|c| HpxCell::new(c.depth, c.hash))
    .collect();
  let moc2 = OwnedRangeMOC::from_it(
    moc.get_depth_max(),
    LazyMOCIter::new(
      moc.get_depth_max(),
      moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
    )
    .to_range_iter(),
  );
  // println!("N ranges: {} {}", moc1.len(), &moc2.to_range_iter().count());
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = moc2.into_range_moc_iter().to_moc_iter().collect();
  println!(
    "To range + collect of size: {} done in {} ms",
    moc1.len(),
    now.elapsed().as_millis()
  );
  assert_eq!(moc1.len(), moc3.len());
  for (
    HpxCell { depth, hash },
    HpxCell {
      depth: depth2,
      hash: hash2,
    },
  ) in moc1.into_iter().zip(&mut moc3.into_iter())
  {
    assert_eq!(depth, depth2);
    assert_eq!(hash, hash2);
  }
}

fn test_expand() {
  let moc = load_sdss().unwrap().to_bmoc();
  // let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc1 = OwnedMOC::from_it_unchecked(
    moc.get_depth_max(),
    moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  );
  let len = moc1.len();
  let now = Instant::now();
  let moc2: Vec<HpxCell<u64>> = moc1.expand().collect();
  println!(
    "Expand + collect of sizes: {} => {} at depth {} done in {} ms",
    len,
    moc2.len(),
    moc.get_depth_max(),
    now.elapsed().as_millis()
  );
}

fn test_expand_v2() {
  let moc = load_sdss().unwrap().to_bmoc();
  // let moc1: Vec<HpxCell<u64>> = moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)).collect();
  let moc1 = OwnedMOC::from_it_unchecked(
    moc.get_depth_max(),
    moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  );
  let moc2 = OwnedMOC::from_it_unchecked(
    moc.get_depth_max(),
    moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  );
  let len = moc1.len();
  let now = Instant::now();
  let moc3: Vec<HpxCell<u64>> = or(moc1.into_moc_iter(), moc2.external_border()).collect();
  println!(
    "Expand + collect of sizes: {} => {} at depth {} done in {} ms",
    len,
    moc3.len(),
    moc.get_depth_max(),
    now.elapsed().as_millis()
  );
}

fn test_expand_ranges() {
  let moc = load_sdss().unwrap().to_bmoc();
  let moc1 = OwnedMOC::<u64>::from_it_unchecked(
    moc.get_depth_max(),
    moc.into_iter().map(|c| HpxCell::new(c.depth, c.hash)),
  );
  let len = moc1.len();
  let ranges = OwnedRangeMOC::from_it(moc.get_depth_max(), moc1.into_moc_iter().to_range_iter());
  let now = Instant::now();
  let moc2: Vec<HpxCell<u64>> = ranges.expand().collect();
  println!(
    "Expand from ranges + collect of sizes: {} => {} at depth {} done in {} ms",
    len,
    moc2.len(),
    moc.get_depth_max(),
    now.elapsed().as_millis()
  );
}

/*
fn compute_dmin() {
  use cdshealpix::{
    compass_point::{Cardinal, Cardinal::E, CardinalSet},
    haversine_dist,
    nested::{get, internal_edge_northeast, path_along_cell_side, Layer},
    nside,
  };

  let mut east = CardinalSet::new();
  east.set(E, true);

  let n_segments = 1000;

  for depth in 0..29 {
    let mut dmin = 1.0;
    let mut imin = 0;
    let mut jmin = 0;
    let layer = get(depth);
    // let nside = nside(depth);
    for hash in internal_edge_northeast(0, depth).iter() {
      // let vertex_map = layer.vertices_map(*hash, east);
      // let vertex_east = vertex_map.get(E).unwrap();
      let edge_eastsouth =
        path_along_cell_side(depth, *hash, &Cardinal::E, &Cardinal::S, true, n_segments);
      let edge_northwest =
        path_along_cell_side(depth, *hash, &Cardinal::N, &Cardinal::W, true, n_segments);
      for (i, posi) in edge_eastsouth.iter().enumerate() {
        for (j, posj) in edge_northwest.iter().enumerate() {
          // let d = haversine_dist(vertex_east.0, vertex_east.1, pos.0, pos.1);
          let d = haversine_dist(posi.0, posi.1, posj.0, posj.1);
          if d < dmin {
            dmin = d;
            imin = i;
            jmin = j;
          }
        }
      }
    }
    println!(
      "depth: {}; dist: {}; imin: {}; ifrac: {}; jmin: {}; jfrac: {}",
      depth,
      dmin,
      imin,
      imin as f64 / n_segments as f64,
      jmin,
      jmin as f64 / n_segments as f64
    );
  }
}
depth: 0; dist: 0.8410686705679302; imin: 0; ifrac: 0; jmin: 0; jfrac: 0
depth: 1; dist: 0.3772364221970879; imin: 0; ifrac: 0; jmin: 265; jfrac: 0.265
depth: 2; dist: 0.1820336496604032; imin: 0; ifrac: 0; jmin: 275; jfrac: 0.275
depth: 3; dist: 0.08911455861353558; imin: 0; ifrac: 0; jmin: 287; jfrac: 0.287
depth: 4; dist: 0.04398973977023659; imin: 0; ifrac: 0; jmin: 289; jfrac: 0.289
depth: 5; dist: 0.021817370982437967; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288
depth: 6; dist: 0.010854009738967782; imin: 0; ifrac: 0; jmin: 287; jfrac: 0.287
depth: 7; dist: 0.005409888771503098; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288
depth: 8; dist: 0.002699583333229803; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288
depth: 9; dist: 0.0013481075907579773; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288
depth: 10; dist: 0.0006735241875226537; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288
depth: 11; dist: 0.00033659544768602315; imin: 0; ifrac: 0; jmin: 288; jfrac: 0.288


 */

fn compute_dmin() {
  use cdshealpix::{
    compass_point::{Cardinal, Cardinal::E, CardinalSet},
    haversine_dist, n_hash,
    nested::{get, internal_edge_northeast, path_along_cell_side, Layer},
    nside,
  };

  let mut east = CardinalSet::new();
  east.set(E, true);

  let mut percentage = 16703484 as f64 / nside(24) as f64; // 4164926 as f64 / nside(22) as f64;

  let n_segments = 10000; // DO NOT MODIFY!! Linked to iseq here after!!
  for depth in 0..30 {
    let mut dmin = 1.0;
    let mut imin = 0;
    let mut ihashmin = 0;
    let layer = get(depth);
    let nside = nside(depth);
    let internal_edge_northeast_cells = internal_edge_northeast(0, depth);
    if depth <= 24 {
      let (icell, iseg) = match depth {
        0 => (0, 0),
        1 => (0, 2646),
        2 => (1, 2750),
        3 => (3, 2873),
        4 => (8, 2892),
        5 => (20, 2876),
        6 => (46, 2870),
        7 => (99, 2878),
        8 => (210, 2880),
        9 => (439, 2882),
        10 => (909, 2882),
        11 => (1865, 2883),
        12 => (3806, 2883),
        13 => (7731, 2884),
        14 => (15653, 2884),
        15 => (31607, 2884),
        16 => (63694, 2884),
        17 => (128149, 2884),
        18 => (257501, 2884),
        19 => (516929, 2884),
        20 => (1036879, 2884),
        21 => (2078602, 2884),
        22 => (4164926, 2884),
        23 => (8342210, 2884),
        24 => (16703484, 2884),
        _ => unreachable!(),
      };
      let hash = internal_edge_northeast_cells[icell];
      let vertex_map = layer.vertices_map(hash, east);
      let vertex_east = vertex_map.get(E).unwrap();
      let edge_northwest =
        path_along_cell_side(depth, hash, &Cardinal::N, &Cardinal::W, true, n_segments);
      let pos = edge_northwest[iseg];
      dmin = haversine_dist(vertex_east.0, vertex_east.1, pos.0, pos.1);
      imin = iseg;
      ihashmin = icell;
    } else {
      for (ihash, hash) in internal_edge_northeast_cells
        .iter()
        .enumerate()
        .skip((nside as f64 * percentage) as usize)
      {
        let vertex_map = layer.vertices_map(*hash, east);
        let vertex_east = vertex_map.get(E).unwrap();
        let edge_northwest =
          path_along_cell_side(depth, *hash, &Cardinal::N, &Cardinal::W, true, n_segments);
        for (i, pos) in edge_northwest.iter().enumerate() {
          // let pos = edge_northwest[2884];
          let d = haversine_dist(vertex_east.0, vertex_east.1, pos.0, pos.1);
          if d < dmin {
            dmin = d;
            imin = i;
            ihashmin = ihash;
            percentage = ihashmin as f64 / nside as f64;
          }
        }
      }
    }
    println!(
      "depth: {}; dist: {:e}; imin: {}; ifrac: {}; ihash: {}; ihashfrac: {}",
      depth,
      dmin,
      imin,
      imin as f64 / n_segments as f64,
      ihashmin,
      ihashmin as f64 / nside as f64
    );
  }
  /*
  depth: 0; dist: 8.410686705679302e-1; imin: 0; ifrac: 0; ihash: 0; ihashfrac: 0
  depth: 1; dist: 3.772363172221682e-1; imin: 2646; ifrac: 0.2646; ihash: 0; ihashfrac: 0
  depth: 2; dist: 1.820336496604032e-1; imin: 2750; ifrac: 0.275; ihash: 1; ihashfrac: 0.25
  depth: 3; dist: 8.911454164164473e-2; imin: 2873; ifrac: 0.2873; ihash: 3; ihashfrac: 0.375
  depth: 4; dist: 4.3989734587906136e-2; imin: 2892; ifrac: 0.2892; ihash: 8; ihashfrac: 0.5
  depth: 5; dist: 2.1817362566999756e-2; imin: 2876; ifrac: 0.2876; ihash: 20; ihashfrac: 0.625
  depth: 6; dist: 1.0854009738967782e-2; imin: 2870; ifrac: 0.287; ihash: 46; ihashfrac: 0.71875
  depth: 7; dist: 5.409888146009395e-3; imin: 2878; ifrac: 0.2878; ihash: 99; ihashfrac: 0.7734375
  depth: 8; dist: 2.699583333229803e-3; imin: 2880; ifrac: 0.288; ihash: 210; ihashfrac: 0.8203125
  depth: 9; dist: 1.3481074890890777e-3; imin: 2882; ifrac: 0.2882; ihash: 439; ihashfrac: 0.857421875
  depth: 10; dist: 6.735240936621056e-4; imin: 2882; ifrac: 0.2882; ihash: 909; ihashfrac: 0.8876953125
  depth: 11; dist: 3.365953703468212e-4; imin: 2883; ifrac: 0.2883; ihash: 1865; ihashfrac: 0.91064453125
  depth: 12; dist: 1.6824522036126843e-4; imin: 2883; ifrac: 0.2883; ihash: 3806; ihashfrac: 0.92919921875
  depth: 13; dist: 8.410609068655501e-5; imin: 2884; ifrac: 0.2884; ihash: 7731; ihashfrac: 0.9437255859375
  depth: 14; dist: 4.204784323221455e-5; imin: 2884; ifrac: 0.2884; ihash: 15653; ihashfrac: 0.95538330078125
  depth: 15; dist: 2.1022283308073835e-5; imin: 2884; ifrac: 0.2884; ihash: 31607; ihashfrac: 0.964569091796875
  depth: 16; dist: 1.0510625670060442e-5; imin: 2884; ifrac: 0.2884; ihash: 63694; ihashfrac: 0.971893310546875
  depth: 17; dist: 5.255150320257332e-6; imin: 2884; ifrac: 0.2884; ihash: 128149; ihashfrac: 0.9776992797851563
  depth: 18; dist: 2.6275239729465538e-6; imin: 2884; ifrac: 0.2884; ihash: 257501; ihashfrac: 0.9822883605957031
  depth: 19; dist: 1.3137458638808036e-6; imin: 2884; ifrac: 0.2884; ihash: 516929; ihashfrac: 0.9859638214111328
  depth: 20; dist: 6.568678535571394e-7; imin: 2884; ifrac: 0.2884; ihash: 1036879; ihashfrac: 0.9888448715209961
  depth: 21; dist: 3.284323270983175e-7; imin: 2884; ifrac: 0.2884; ihash: 2078602; ihashfrac: 0.991154670715332
  depth: 22; dist: 1.642156595517884e-7; imin: 2884; ifrac: 0.2884; ihash: 4164926; ihashfrac: 0.9929957389831543
  depth: 23; dist: 8.21076709163242e-8; imin: 2884; ifrac: 0.2884; ihash: 8342210; ihashfrac: 0.9944689273834229
  depth: 24; dist: 4.105378528139296e-8; imin: 2884; ifrac: 0.2884; ihash: 16703484; ihashfrac: 0.9956052303314209
  depth: 25; dist: 2.0526876713226626e-8; imin: 2884; ifrac: 0.2884; ihash: 33436349; ihashfrac: 0.9964808523654938
  depth: 26; dist: 1.0263433216329513e-8; imin: 2884; ifrac: 0.2884; ihash: 66923451; ihashfrac: 0.9972371309995651
  depth: 27; dist: 5.131714858175969e-9; imin: 2884; ifrac: 0.2884; ihash: 133912273; ihashfrac: 0.9977241829037666
  depth: 28; dist: 2.5658567623093986e-9; imin: 2884; ifrac: 0.2884; ihash: 267902108; ihashfrac: 0.9980131238698959
  depth: 29; dist: 1.2829280665188905e-9; imin: 2883; ifrac: 0.2883; ihash: 536040904; ihashfrac: 0.9984539896249771
     */
}

pub fn main() {
  // compute_dmin();

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
