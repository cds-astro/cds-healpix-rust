//! Contains the code for MOC compression/decompression one can made on streamed MOCs
//! (i.e. based on uniq cells iterators ordered following the natural Z-order curve order).

// use std::cmp::{min, max, Ordering};
use super::{CheckedIterator, HasMaxDepth, HpxCell, HpxHash, MOCIterator};

///
/// # Panic
/// * this method panics if the input iterator is not ordered following the natural Z-order
/// curve ordering.
pub fn compress<H, T>(it: T) -> CompressMocIter<H, CheckedIterator<H, T>>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  let it = CheckedIterator::new(it);
  compress_unchecked(it)
}

pub fn compress_unchecked<H, T>(it: T) -> CompressMocIter<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  CompressMocIter::new_unchecked(it)
}

/// Uncompress on-the-fly, while iterating, the MOC compressed in the given slice.
/// For more information, see the `decompress` method in BMOCs.
pub fn uncompress<H: HpxHash>(compressed: &[u8]) -> UncompressMocIter<H> {
  UncompressMocIter::new(compressed)
}

/// Performs an `Compression` operation on-the-fly, on the input MOC, while iterating.
pub struct CompressMocIter<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  it: T,
  curr: Option<HpxCell<H>>,
  d: u8,
  h: H,
  buff: u32, // Go down worst case: from depth 0 to 29 => 3*30*2bits + 29*1 bit => 27 Bytes
  ibit: u8,
  is_going_down: bool,
}

impl<H, T> CompressMocIter<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  ///
  /// # WARNING
  /// * The input iterator **MUST** be ordered following the natural Z-order curve order.
  /// If it is not the case, the result will be wrong, without errors!
  /// * If you are not sure, use the `new` version.
  fn new_unchecked(mut it: T) -> CompressMocIter<H, T> {
    let depth_max = it.depth_max();
    let curr = it.next();
    CompressMocIter {
      it,
      curr,
      d: 0,
      h: H::zero(),
      buff: depth_max as u32,
      ibit: 8,
      is_going_down: true,
    }
  }

  fn push_0(&mut self) {
    self.ibit += 1;
  }
  /// n: number of 0s to be pushed
  fn push_0_x(&mut self, n: u8) {
    self.ibit += n;
  }
  fn push_1(&mut self) {
    self.buff |= 1_u32 << self.ibit;
    self.push_0();
  }
  /// Returns true if the first byte of the buffer is full
  fn is_byte_full(&self) -> bool {
    self.ibit > 7
  }
  /// Take the first byte if the buffer an replace its content with the second byte
  fn take_byte(&mut self) -> u8 {
    let r = self.buff as u8;
    self.buff >>= 8;
    self.ibit -= 8;
    r
  }

  /// Tells is there is a last byte to be taken
  fn has_last_byte(&self) -> bool {
    debug_assert!(self.ibit < 8);
    self.ibit > 0
  }

  /// Take the last byte (WARNING: to be called at the end!!)
  fn take_last_byte(&mut self) -> u8 {
    debug_assert!(self.ibit < 8);
    let r = self.buff as u8;
    self.buff >>= 8;
    self.ibit = 0;
    r
  }

  /// Returns true if the first byte of the buffer is full
  /*fn push_node_empty(&mut self) -> bool {
    self.push_node_empty_notest();
    self.is_byte_full()
  }
  fn push_node_empty_notest(&mut self) {
    self.push_1();
    self.push_0();
  }*/
  fn push_n_empty_nodes_notest(&mut self, n: H) {
    for _ in 0..n.to_u8().unwrap() {
      self.push_1();
      self.push_0();
    }
  }
  /// Returns true if the first byte of the buffer is full
  /*fn push_node_full(&mut self) -> bool {
    self.push_node_full_notest();
    self.is_byte_full()
  }*/
  fn push_node_full_notest(&mut self) {
    self.push_1();
    self.push_1();
  }

  /// Returns true if the first byte of the buffer is full
  fn push_node_partial(&mut self) -> bool {
    self.push_node_partial_notest();
    self.is_byte_full()
  }
  fn push_node_partial_notest(&mut self) {
    self.push_0();
  }
  /// Returns true if the first byte of the buffer is full
  /*fn push_leaf_empty(&mut self) -> bool {
    self.push_leaf_empty_notest();
    self.is_byte_full()
  }
  fn push_leaf_empty_notest(&mut self) {
    self.push_0();
  }*/
  fn push_n_empty_leaves_notest(&mut self, n: H) {
    self.push_0_x(n.to_u8().unwrap());
  }
  /// Returns true if the first byte of the buffer is full
  /*fn push_leaf_full(&mut self) -> bool {
    self.push_leaf_full_notest();
    self.is_byte_full()
  }*/
  fn push_leaf_full_notest(&mut self) {
    self.push_1();
  }
}

impl<H, T> HasMaxDepth for CompressMocIter<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl<H, T> Iterator for CompressMocIter<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  type Item = u8;

  // WE CAN PROBABLY SIMPLIFY THIS CODE!!
  fn next(&mut self) -> Option<Self::Item> {
    let three = H::one() | (H::one() << 1);
    if self.is_byte_full() {
      return Some(self.take_byte());
    }
    match &self.curr {
      Some(hpx) => {
        let (curr_d, curr_h) = (hpx.depth, hpx.hash);
        if !self.is_going_down {
          // go up (if needed)!
          let target_h = if self.d > curr_d {
            // case previous hash deeper that current hash
            let dd = self.d - curr_d;
            curr_h << ((dd << 1) as usize)
          } else {
            // case current hash deeper (or same depth) that previous hash, need to go up?
            let dd = curr_d - self.d;
            curr_h >> ((dd << 1) as usize)
          };
          let dd = ((63 - (self.h ^ target_h).leading_zeros()) >> 1) as u8;
          let target_d = if dd > self.d { 0 } else { self.d - dd };
          // - go up to common depth (if needed)
          while self.d > target_d {
            let n = three - (self.h & three);
            debug_assert!(n <= three);
            if self.d == self.depth_max() {
              self.push_n_empty_leaves_notest(n); // Writes max 3 0s
            } else {
              self.push_n_empty_nodes_notest(n); // Writes max 6 0s
            }
            self.h >>= 2;
            self.d -= 1;
            if self.is_byte_full() {
              return Some(self.take_byte());
            }
          }
          self.h += H::one();
        }
        // - go down
        self.is_going_down = true;
        let target_d = curr_d;
        while self.d < target_d {
          let dd = target_d - self.d;
          let target_h = curr_h >> ((dd << 1) as usize);
          let n = target_h - self.h;
          self.h = target_h << 2;
          self.d += 1;
          debug_assert!(n <= three);
          self.push_n_empty_nodes_notest(n); // Writes max 6 0s
          if self.push_node_partial() {
            return Some(self.take_byte());
          }
        }
        // - same level
        debug_assert_eq!(self.d, target_d);
        let n = curr_h - self.h;
        debug_assert!(
          n <= three,
          "n: {}, curr_h: {}, self.h: {}",
          n,
          curr_h,
          self.h
        );
        if self.d == self.depth_max() {
          self.push_n_empty_leaves_notest(n);
          self.push_leaf_full_notest()
        } else {
          self.push_n_empty_nodes_notest(n);
          self.push_node_full_notest();
        }
        self.h = curr_h;
        self.curr = self.it.next();
        self.is_going_down = false;
        self.next()
      }
      None => {
        // Just need to complete
        if self.d == 0 {
          // complete depth 0 till hash 12
          let twelve = three << 2_usize;
          let eleven = twelve - H::one();
          if self.h < eleven {
            let n = eleven - self.h;
            self.h = twelve;
            if self.d == self.depth_max() {
              self.push_n_empty_leaves_notest(n); // writes max 11 0s
            } else {
              self.push_n_empty_nodes_notest(n); // writes max 22 0s
            }
            self.next()
          } else if self.has_last_byte() {
            debug_assert!(
              self.d == 0 && self.h == eleven,
              "d: {}, h: {}",
              self.d,
              self.h
            );
            Some(self.take_last_byte())
          } else {
            debug_assert!(
              self.d == 0 && self.h == eleven,
              "d: {}, h: {}",
              self.d,
              self.h
            );
            None
          }
        } else {
          // go up to depth 0
          while {
            let n = three - (self.h & three);
            if self.d == self.depth_max() {
              self.push_n_empty_leaves_notest(n); // writes max 3 0s
            } else {
              self.push_n_empty_nodes_notest(n); // writes max 6 0s
            }
            self.h >>= 2;
            self.d -= 1;
            if self.is_byte_full() {
              return Some(self.take_byte());
            }
            self.d > 0
          } {} // do-while loop
          self.next()
        }
      }
    }
  }

  /*fn size_hint(&self) -> (usize, Option<usize>) {
    let size_hint = self.it.size_hint();
    if let Some(n) = size_hint.1 {
      (0, 4 + 3 * n)
    } else {
      (0, None)
    }
  }*/
}

/// Performs an `Uncompression` operation on-the-fly, while iterating.
pub struct UncompressMocIter<'a, H: HpxHash> {
  /// Slice containing the compressed MOC
  cmoc: &'a [u8],
  depth_max: u8,
  ibyte: usize,
  ibit: u8,
  depth: u8,
  hash: H,
}

impl<'a, H: HpxHash> UncompressMocIter<'a, H> {
  ///
  /// # WARNING
  /// * The input iterator **MUST** be ordered following the natural Z-order curve order.
  /// If it is not the case, the result will be wrong, without errors!
  /// * If you are not sure, use the `new` version.
  fn new(cmoc: &'a [u8]) -> UncompressMocIter<'a, H> {
    let depth_max = cmoc[0];
    UncompressMocIter {
      cmoc,
      depth_max,
      ibyte: 1,
      ibit: 0,
      depth: 0,
      hash: H::zero(),
    }
  }

  fn get(&mut self) -> bool {
    let r = (self.cmoc[self.ibyte] & (1_u8 << self.ibit)) != 0;
    self.ibyte += (self.ibit == 7) as usize;
    self.ibit += 1;
    self.ibit &= 7;
    r
  }

  fn get_current_cell(&self) -> HpxCell<H> {
    HpxCell {
      depth: self.depth,
      hash: self.hash,
    }
  }

  fn increment(&mut self) {
    self.go_up_if_needed();
    self.hash += H::one();
  }

  fn go_up_if_needed(&mut self) {
    let three = H::one() | (H::one() << 1);
    while self.hash & three == three && self.depth > 0 {
      self.hash >>= 2;
      self.depth -= 1;
    }
  }

  /* Both the above and this solution are equivalent.
     This one is used in NotMocIter (see op module).
     Test performances on a very large MOC?
  fn increment(&mut self) {
    self.hash += H::one();
    self.go_up_if_needed();
  }

  fn go_up_if_needed(&mut self, ) {
    let dd = self.depth.min(self.hash.trailing_zeros() as u8 / 2).min(self.depth);
    self.depth -= dd;
    self.hash >>= dd << 1;
  }*/

  fn go_down_of_1_depth(&mut self) {
    self.hash <<= 2;
    self.depth += 1;
  }
}

impl<'a, H> HasMaxDepth for UncompressMocIter<'a, H>
where
  H: HpxHash,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<'a, H: HpxHash> Iterator for UncompressMocIter<'a, H> {
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    let three = H::one() | (H::one() << 1);
    let twelve = three << 2_usize;
    while self.depth != 0 || self.hash != twelve {
      if self.get() {
        // bit = 1
        if self.depth == self.depth_max || self.get() {
          let res = Some(self.get_current_cell());
          self.increment();
          return res;
        } else {
          self.increment();
        }
      } else if self.depth == self.depth_max {
        // bit == 0
        self.increment();
      } else {
        // bit == 0
        debug_assert!(self.depth < self.depth_max);
        self.go_down_of_1_depth();
      }
    }
    None
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let n_remains = 8 * (self.cmoc.len() - self.ibyte) - (self.ibit as usize);
    (0, Some(n_remains))
  }
}
