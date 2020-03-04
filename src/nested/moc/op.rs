//! Contains the code of logical operations one can made on streamed MOCs
//! (i.e. based on uniq cells iterators ordered following the natural Z-order curve order).

use std::cmp::{min, max, Ordering};
use std::mem::replace;

use super::{HpxHash, HpxCell, HasMaxDepth, MOCIterator, CheckedIterator};


///
/// # Panic
/// * this method panics if at least one of the two input iterator is not ordered following
/// the natural Z-order curve ordering.
pub fn not<H, T>(it: T) -> NotMocIter<H, CheckedIterator<H, T>>
  where H: HpxHash,
        T: MOCIterator<H> {
  let it = CheckedIterator::new(it);
  not_unchecked(it)
}

pub fn not_unchecked<H, T>(it: T) -> NotMocIter<H, T>
  where H: HpxHash,
        T: MOCIterator<H> {
  NotMocIter::new_unchecked(it)
}


///
/// # Panic
/// * this method panics if at least one of the two input iterator is not ordered following
/// the natural Z-order curve ordering.
pub fn and<H, T1, T2>(left_it: T1, right_it: T2) 
  -> AndMocIter<H, CheckedIterator<H, T1>, CheckedIterator<H, T2>>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  let left_it = CheckedIterator::new(left_it);
  let right_it = CheckedIterator::new(right_it);
  and_unchecked(left_it, right_it)
}

pub fn and_unchecked<H, T1, T2>(left_it: T1, right_it: T2) -> AndMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  AndMocIter::new_unchecked(left_it, right_it)
}

///
/// # Panic
/// * this method panics if at least one of the two input iterator is not ordered following
/// the natural Z-order curve ordering.
pub fn or<H, T1, T2>(left_it: T1, right_it: T2)
  -> OrMocIter<H, CheckedIterator<H, T1>, CheckedIterator<H, T2>>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  let left_it = CheckedIterator::new(left_it);
  let right_it = CheckedIterator::new(right_it);
  or_unchecked(left_it, right_it)
}

pub fn or_unchecked<H, T1, T2>(left_it: T1, right_it: T2) -> OrMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  OrMocIter::new_unchecked(left_it, right_it)
}


/// Performs an `NOT` operation between two MOCs on-the-fly, while iterating.
/// # Remark:
/// * If you have an array in memory, the *NOT* operation is very fast and contains the number of 
/// ranges in the MOC plus 1.
/// * Operations on ranges is simpler and faster than operations on unique cells. But when the MOC
/// is stored in (compressed) uniq cells instead of ranges, we may account for the time needed to
/// convert from uniq cells to ranges.
pub struct NotMocIter<H, T>
  where H: HpxHash,
        T: MOCIterator<H> {
  it: T,
  curr: Option<HpxCell<H>>,
  curr_d: u8,
  curr_h: H,
}

impl <H, T> NotMocIter<H, T>
  where H: HpxHash,
        T: MOCIterator<H> {
  ///
  /// # WARNING
  /// * Each input iterator **MUST** be ordered following the natural Z-order curve order.
  /// If it is not the case, the result will be wrong, without errors!
  /// * If you are not sure, use the `new` version. 
  fn new_unchecked(mut it: T) -> NotMocIter<H, T> {
    let mut curr = it.next();
    NotMocIter {
      it,
      curr,
      curr_d: 0,
      curr_h: H::zero(),
    }
  }
  
  fn get_current_cell(&self) -> HpxCell<H> {
    HpxCell {
      depth: self.curr_d,
      hash: self.curr_h,
    }
  }
  
  fn next_cell(&mut self) {
    self.curr_h += H::one();
    self.go_up_if_needed();
  }
  
  fn go_up_if_needed(&mut self, ) {
    let dd = self.curr_d.min(self.curr_h.trailing_zeros() as u8 / 2).min(self.curr_d);
    self.curr_d -= dd;
    self.curr_h >>= dd << 1;
  }
  
}

impl<H, T> HasMaxDepth for NotMocIter<H, T>
  where H: HpxHash ,
        T: MOCIterator<H> {
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl <H, T> Iterator for NotMocIter<H, T>
  where H: HpxHash,
        T: MOCIterator<H> {
  type Item = HpxCell<H>;
  
  fn next(&mut self) -> Option<Self::Item> {
    match &self.curr {
      Some(c) => {
        // while the current cell contains the target cell, go down
        while self.curr_d < c.depth && self.curr_h == (c.hash >> (((c.depth - self.curr_d) << 1)) as usize) {
          self.curr_d += 1;
          self.curr_h <<= 2;
        }
        // if current cell and target cell are the same, continue
        if self.curr_d == c.depth && self.curr_h == c.hash {
          self.next_cell();
          self.curr = self.it.next();
          self.next()
        } else { // else return the current cell
          let c = self.get_current_cell();
          self.next_cell();
          Some(c)
        }
      },
      None => { 
        let twelve = (H::one() + (H::one() << 1)) << 2_usize; // Transform this in const when available in Rust
        if self.curr_h < twelve << ((self.curr_d << 1) as usize) {
          let c = self.get_current_cell();
          self.next_cell();
          Some(c)
        } else {
          None
        }
      },
    }
  }

  /*fn size_hint(&self) -> (usize, Option<usize>) {
    let size_hin = self.it.size_hint();
    if let Some(n) = size_hint.1 {
      (0, 3 * n + 12) // more than that: account for depth!!
    } else {
      (0, None)
    }
  }*/
}



/// Performs an `AND` operation between two MOCs on-the-fly, while iterating.
/// # Remark:
/// * If you have an array in memoryn the fastest *AND* operation consists in starting with 2 binary
/// searches to isolate the smallest common segment. This is not possible when working with iterators.
/// * Operations on ranges is simpler and faster than operations on unique cells. But when the MOC
/// is stored in (compressed) uniq cells instead of ranges, we may account for the time needed to
/// convert from uniq cells to ranges.  
pub struct AndMocIter<H, T1, T2> 
  where H: HpxHash,
       T1: MOCIterator<H>,
       T2: MOCIterator<H> {
  left_it: T1,
  right_it: T2,
  depth_max: u8,
  left: Option<HpxCell<H>>,
  right: Option<HpxCell<H>>,
}

impl <H, T1, T2> AndMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  ///
  /// # WARNING
  /// * Each input iterator **MUST** be ordered following the natural Z-order curve order.
  /// If it is not the case, the result will be wrong, without errors!
  /// * If you are not sure, use the `new` version. 
  fn new_unchecked(mut left_it: T1, mut right_it: T2) -> AndMocIter<H, T1, T2> {
    let depth_max = max(left_it.depth_max(), right_it.depth_max());
    let left = left_it.next();
    let right = right_it.next();
    AndMocIter {
      left_it,
      right_it,
      depth_max,
      left,
      right,
    }
  }
  
}

impl<H, T1, T2> HasMaxDepth for AndMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl <H, T1, T2> Iterator for AndMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {

  type Item = HpxCell<H>;
  
  fn next(&mut self) -> Option<Self::Item> {
    while let (Some(l), Some(r)) = (&self.left, &self.right) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          let hr_at_dl = r.hash.unsigned_shr(((r.depth - l.depth) << 1) as u32);
          match l.hash.cmp(&hr_at_dl) {
            Ordering::Less => self.left = self.left_it.next(),
            Ordering::Greater => self.right = self.right_it.next(),
            Ordering::Equal => { return replace(&mut self.right, self.right_it.next()); },
          }
        },
        Ordering::Greater => {
          let hl_at_dr = l.hash.unsigned_shr(((l.depth - r.depth) << 1) as u32);
          match hl_at_dr.cmp(&r.hash) {
            Ordering::Less => self.left = self.left_it.next(),
            Ordering::Greater => self.right = self.right_it.next(),
            Ordering::Equal => { return replace(&mut self.left, self.left_it.next()); },
          }
        },
        Ordering::Equal => {
          match l.hash.cmp(&r.hash) {
            Ordering::Less => self.left = self.left_it.next(),
            Ordering::Greater => self.right = self.right_it.next(),
            Ordering::Equal => {
              self.left = self.left_it.next();
              return replace(&mut self.right, self.right_it.next());
            },
          }
        },
      }
      /* This (equivalent) version seems to be slighlty faster, to be further tested!!
      if l.depth < r.depth {
        let hr_at_dl = r.hash.unsigned_shr(((r.depth - l.depth) << 1) as u32);
        if l.hash < hr_at_dl {
          self.left = self.left_it.next();
        } else if l.hash > hr_at_dl {
          self.right = self.right_it.next();
        } else {
          return replace(&mut self.right, self.right_it.next());
        }
      } else if l.depth > r.depth {
        let hl_at_dr = l.hash.unsigned_shr(((l.depth - r.depth) << 1) as u32);
        if hl_at_dr < r.hash {
          self.left = self.left_it.next();
        } else if hl_at_dr > r.hash {
          self.right = self.right_it.next();
        } else {
          return replace(&mut self.left, self.left_it.next());
        }
      } else {
        debug_assert_eq!(l.depth, r.depth);
        if l.hash < r.hash {
          self.left = self.left_it.next();
        } else if l.hash > r.hash  {
          self.right = self.right_it.next();
        } else {
          self.left = self.left_it.next();
          return replace(&mut self.right, self.right_it.next());
        }
      }*/
    }
    None
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let size_hint_l = self.left_it.size_hint();
    let size_hint_r = self.right_it.size_hint();
    if let (Some(n1), Some(n2)) = (size_hint_l.1, size_hint_r.1) {
      // each MOC may contains low res cells containing a lot of the other MOC cells
      (0, Some(1 + n1 + n2))
    } else {
      (0, None)
    }
  }

}

/// Performs an `OR` operation between two MOCs on-the-fly, while iterating.
pub struct OrMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  left_it: T1,
  right_it: T2,
  depth_max: u8,
  left: Option<HpxCell<H>>,
  right: Option<HpxCell<H>>,
  to_be_possibly_packed: Vec<HpxCell<H>>,
  flush_stack: u8,
}

impl <H, T1, T2> OrMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  ///
  /// # WARNING
  /// * Each input iterator **MUST** be ordered following the natural Z-order curve order.
  /// If it is not the case, the result will be wrong, without errors!
  /// * If you are not sure, use the `new` version. 
  fn new_unchecked(mut left_it: T1, mut right_it: T2) -> OrMocIter<H, T1, T2> {
    let depth_max = max(left_it.depth_max(), right_it.depth_max());
    let left = left_it.next();
    let right = right_it.next();
    OrMocIter {
      left_it,
      right_it,
      depth_max,
      left,
      right,
      to_be_possibly_packed: Default::default(),
      flush_stack: 0, // number of elements to be flushed in the to_be_possibly_packed stack
    }
  }

  pub fn depth_max(&self) -> u8 {
    self.depth_max
  }
  
  fn pack_stack(&mut self) {
    let h0 = H::zero();
    let h1 = H::one();
    let h2 = h1 << 1;
    let h3 = h2 | h1;
    while self.to_be_possibly_packed.len() >= 4 {
      let len = &self.to_be_possibly_packed.len();
      let one = &self.to_be_possibly_packed[len - 4];
      let two = &self.to_be_possibly_packed[len - 3];
      let thr = &self.to_be_possibly_packed[len - 2];
      let fou = &self.to_be_possibly_packed[len - 1];
      if one.hash & h3 == h0
        && one.depth == two.depth 
        && one.depth == thr.depth 
        && one.depth == fou.depth
        && one.hash + h1 == two.hash
        && one.hash + h2 == thr.hash
        && one.hash + h3 == fou.hash {
        let depth = one.depth - 1;
        let hash = one.hash >> 2;
        self.to_be_possibly_packed.truncate(len - 4);
        self.to_be_possibly_packed.push(HpxCell { depth, hash });
        if depth == 0 {
          debug_assert_eq!(self.to_be_possibly_packed.len(), 1);
          self.flush_stack = 1;
        }
      } else {
        break;
      }
    }
  }
  
  fn unpacked_next(&mut self) -> Option<HpxCell<H>> {
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    match (&self.left, &self.right) {
      (Some(l), Some(r)) =>
        match l.depth.cmp(&r.depth) {
          Ordering::Less => {
            let hr_at_dl = r.hash.unsigned_shr(((r.depth - l.depth) << 1) as u32);
            match l.hash.cmp(&hr_at_dl) {
              Ordering::Less => replace(&mut self.left, self.left_it.next()),
              Ordering::Greater => replace(&mut self.right, self.right_it.next()),
              Ordering::Equal => {
                // Consume right elements overlapped by l
                self.right = consume_overlapped_cells(&mut self.right_it, l);
                replace(&mut self.left, self.left_it.next())
              },
            }
          },
          Ordering::Greater => {
            let hl_at_dr = l.hash.unsigned_shr(((l.depth - r.depth) << 1) as u32);
            match hl_at_dr.cmp(&r.hash) {
              Ordering::Less => replace(&mut self.left, self.left_it.next()),
              Ordering::Greater => replace(&mut self.right, self.right_it.next()),
              Ordering::Equal => {
                // Consume left elements overlapped by r
                self.left = consume_overlapped_cells(&mut self.left_it, r);
                replace(&mut self.right, self.right_it.next())
              },
            }
          },
          Ordering::Equal => {
            match l.hash.cmp(&r.hash) {
              Ordering::Less => replace(&mut self.left, self.left_it.next()),
              Ordering::Greater => replace(&mut self.right, self.right_it.next()),
              Ordering::Equal => {
                self.left = self.left_it.next();
                replace(&mut self.right, self.right_it.next())
              },
            }
          },
        },
      (Some(_l), None) => replace(&mut self.left, self.left_it.next()),
      (None, Some(_r)) => replace(&mut self.right, self.right_it.next()),
      (None, None) => None,
    }
  }
}

impl<H, T1, T2> HasMaxDepth for OrMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl <H, T1, T2> Iterator for OrMocIter<H, T1, T2>
  where H: HpxHash,
        T1: MOCIterator<H>,
        T2: MOCIterator<H> {
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    // All this mess to handle cases in which the union of cells has to be merged 
    // in a cell of lower resolution
    let h0 = H::zero();
    let h1 = H::one();
    let h2 = h1 << 1;
    let h3 = h2 & h1;
    if self.flush_stack > 0 {
      self.flush_stack -= 1;
      Some(self.to_be_possibly_packed.remove(0))
    } else {
      let new = self.unpacked_next();
      match new {
        Some(curr) => 
          if self.to_be_possibly_packed.len() > 0 {
            // Unwrapping is safe since we first test that len() > 0
            let prev = self.to_be_possibly_packed.last().unwrap();
            if curr.depth == prev.depth && curr.hash == prev.hash + H::one() {
              let pack = curr.hash & h3 == h3;
              self.to_be_possibly_packed.push(curr);
              if pack {
                self.pack_stack();
              }
            } else if curr.depth > prev.depth && curr.hash == (prev.hash + h1) << (((curr.depth - prev.depth) << 1) as usize) {
              debug_assert_eq!(curr.hash & h3, h0);
              self.to_be_possibly_packed.push(curr);
            } else {
              let curr_d = curr.depth;
              let curr_h = curr.hash;
              self.to_be_possibly_packed.push(curr);
              self.flush_stack = if curr_d > 0 && curr_h & h3 == h0 {
                self.to_be_possibly_packed.len() - 1
              } else {
                self.to_be_possibly_packed.len()
              } as u8;
            }
            self.next()
          } else if curr.depth > 0 && curr.hash & h3 == h0 {
            self.to_be_possibly_packed.push(curr);
            self.next()
          } else {
            Some(curr)
          },
        None => 
          if self.to_be_possibly_packed.len() > 0 {
            self.flush_stack = self.to_be_possibly_packed.len() as u8;
            self.next()
          } else {
            None
          }, 
      }
    }
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let size_hint_l = self.left_it.size_hint();
    let size_hint_r = self.right_it.size_hint();
    match (size_hint_l.1, size_hint_r.1) {
      (Some(n1), Some(n2)) => (1 + min(n1, n2), Some(1 + n1 + n2)),
      (Some(n1), None) => (1 + n1, Some(1 + n1)),
      (None, Some(n2)) => (1 + n2, Some(1 + n2)),
      (None, None) => if self.left.is_some() || self.right.is_some() { (1, Some(1)) } else { (0, Some(0)) },
    }
  }
}

fn consume_overlapped_cells<H, T2>(it: &mut T2, low_res_cell: &HpxCell<H>) -> Option<HpxCell<H>>
  where H: HpxHash,
        T2: Iterator<Item=HpxCell<H>> {
  let mut curr = it.next();
  while let Some(c) = &curr {
    if c.depth <= low_res_cell.depth || low_res_cell.hash != c.hash.unsigned_shr(((c.depth - low_res_cell.depth) << 1) as u32) {
      break;
    } else {
      curr = it.next();
    }
  }
  curr
}

// XOR

// MINUS ??



#[cfg(test)]
mod tests {
  use super::*;
  use super::super::BasicMOCIter;
  
  #[test]
  pub fn test_not_empty() {
    for depth in 0..30 {
      let empty_it = std::iter::empty::<HpxCell<u64>>();
      let mut full_it = not_unchecked(BasicMOCIter::new(depth, empty_it));
      for hash in 0..12_u64 {
        assert_eq!(full_it.next(), Some(HpxCell { depth: 0, hash }));
      }
      assert_eq!(full_it.next(), None);
    }
  }
  
  #[test]
  pub fn test_not_full() {
    for depth in 0..30 {
      let mut empty_it = 
        not_unchecked(BasicMOCIter::new(depth, 
                   (0..12_u64).into_iter().map(|hash| HpxCell{depth: 0, hash}))
      );
      assert_eq!(empty_it.next(), None);
    }
  }

  #[test]
  pub fn test_and_full_full() {
    for depth in 0..30 {
      let mut full_1: Vec<HpxCell<u64>> = (0..12_u64).into_iter()
        .map(|hash| HpxCell{depth: 0, hash}).collect();
      let mut full_2: Vec<HpxCell<u64>> = (0..12_u64).into_iter()
        .map(|hash| HpxCell{depth: 0, hash}).collect();
      let mut full_it = and_unchecked(
        BasicMOCIter::new(depth, full_1.into_iter()),
        BasicMOCIter::new(depth, full_2.into_iter())
      );
      for hash in 0..12_u64 {
        assert_eq!(full_it.next(), Some(HpxCell { depth: 0, hash }));
      }
      assert_eq!(full_it.next(), None);
    }
  }

  #[test]
  pub fn test_and_empty_full() {
    for depth in 0..30 {
      let empty_it = std::iter::empty::<HpxCell<u64>>();
      let mut full: Vec<HpxCell<u64>> = (0..12_u64).into_iter()
        .map(|hash| HpxCell{depth: 0, hash}).collect();
      let mut emtpy_res = and_unchecked(
        BasicMOCIter::new(depth, empty_it),
        BasicMOCIter::new(depth, full.into_iter()),
      );
      assert_eq!(emtpy_res.next(), None);
    }
  }

  #[test]
  pub fn test_or_full_empty() {
    for depth in 0..30 {
      let empty_it = std::iter::empty::<HpxCell<u64>>();
      let mut full: Vec<HpxCell<u64>> = (0..12_u64).into_iter()
        .map(|hash| HpxCell{depth: 0, hash}).collect();
      let mut full_it = or_unchecked(
        BasicMOCIter::new(depth, empty_it),
        BasicMOCIter::new(depth, full.into_iter()),
      );
      for hash in 0..12_u64 {
        assert_eq!(full_it.next(), Some(HpxCell { depth: 0, hash }));
      }
      assert_eq!(full_it.next(), None);
    }
  }

  #[test]
  pub fn test_or_empty_empty() {
    for depth in 0..30 {
      let empty_it_1 = std::iter::empty::<HpxCell<u64>>();
      let empty_it_2 = std::iter::empty::<HpxCell<u64>>();
      let mut emtpy_res = and_unchecked(
        BasicMOCIter::new(depth, empty_it_1),
        BasicMOCIter::new(depth, empty_it_2),
      );
      assert_eq!(emtpy_res.next(), None);
    }
  }
  
}
