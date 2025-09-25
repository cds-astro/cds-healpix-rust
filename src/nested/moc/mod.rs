//! For efficients MOC operations, see the Rust code in [mocpy](https://github.com/cds-astro/mocpy).
//!
//! The purpose of this module is to test (and test the performances) of compression/decompression
//! and logical operations on a streaming mode (i.e. based on uniq cells iterators ordered following
//! the natural Z-order curve order).
//! The ideas are:
//! - those operations are slower than operations on ranges...
//! - ... and can't benefit from optimizations (like binary search of the common range for the
//!   AND operation), **BUT**...
//! - ... they may be faster than: uniqs to ranges + logical op on ranges + ranges to uniqs operations.
//! - do not need all the MOC to be loaded in memory (logical operations could also be implemented
//!   on streams of ordered ranges)
//! - works with compressed MOCs decompressed on-the-fly, and compressing the result on-the-fly
//!   too (with a very few amount of memory needed).

use num::{PrimInt, ToPrimitive, Unsigned};
use num_traits::NumAssign;

use std::cmp::{min, Ordering};
use std::fmt::{Debug, Display};
use std::mem::{self, replace};
use std::ops::{ShlAssign, ShrAssign};
use std::slice::Iter;
use std::vec::{IntoIter, Vec};

use super::{append_external_edge, internal_edge_sorted};

pub mod compressed;
pub mod op;

use compressed::*;
use op::*;

/// Trait defining an HEALPix Hash number.
/// In practice u16, u32 and u64 implement this trait.
pub trait HpxHash:
  NumAssign
  + PrimInt
  + ToPrimitive
  + Unsigned
  + ShlAssign<u8>
  + ShrAssign<u8>
  + Copy
  + Display
  + Debug
{
  const N_BYTES: u8 = mem::size_of::<Self>() as u8;
  const N_BITS: u8 = Self::N_BYTES << 3;
  // 4 bits reserved for base cells
  // + 1 bit for the sentinel bit
  // + 1 bit for the flag in BMOCs
  const MAX_DEPTH: u8 = (((mem::size_of::<Self>() as u8) << 3) - 6) >> 1;
}
impl<H> HpxHash for H where
  H: NumAssign
    + PrimInt
    + ToPrimitive
    + Unsigned
    + ShlAssign<u8>
    + ShrAssign<u8>
    + Copy
    + Display
    + Debug
{
}

/// Represents an HEALPix cell with its depth and hash value.
#[derive(Debug, PartialEq, Eq)]
pub struct HpxCell<H: HpxHash> {
  pub depth: u8,
  pub hash: H,
}

impl<H: HpxHash> HpxCell<H> {
  pub fn new(depth: u8, hash: H) -> HpxCell<H> {
    HpxCell { depth, hash }
  }

  /// Transform this cell in a Range
  pub fn range(&self) -> HpxRange<H> {
    let twice_dd = ((H::MAX_DEPTH - self.depth) << 1) as usize;
    HpxRange::new(self.hash << twice_dd, (self.hash + H::one()) << twice_dd)
  }
}

/// Type defining a HEALPix Range, the value of from and to are the value at `H::MAX_DEPTH`, i.e.
/// depth 29 for `u64`.
/// # Remark:
/// I first wanted to declare:
/// > pub type HpxRange<H: HpxHash> = Range<H>;
///
/// but it was not possible to define a new `Iterator` on it.
pub struct HpxRange<H: HpxHash> {
  /// The lower bound of the range (inclusive).
  pub from: H,
  /// The upper bound of the range (exclusive).
  pub to: H,
}

impl<H: HpxHash> HpxRange<H> {
  pub fn new(from: H, to: H) -> HpxRange<H> {
    HpxRange { from, to }
  }

  /// # Args:
  /// * `depth_max`: the maximum depth of cell in the range (the RangeMOC depth max).
  /// * `twice_dd`: twice the depth difference between DEPTH_MAX and RangeMOC `depth_max`
  /// * `range_len_min`: the length of a range of one cell (depends on the RangeMOC depth max).
  /// * `mask`: a mask allowing to know if the first element of the range if a cell of `depth_max`
  /// # Remark:
  /// We can deduce `twice_dd` and `range_len_min` from `H` and `depth_max`, but we pass them to avoid
  /// recomputing them (even if it is fast)
  fn next_with_knowledge(
    &mut self,
    depth_max: u8,
    twice_dd: usize,
    range_len_min: H,
    mask: H,
  ) -> Option<HpxCell<H>> {
    let len = self.to - self.from;
    if len < H::one() {
      None
    } else if len == range_len_min || self.from & mask != H::zero() {
      // A range of 1 cell at depth_max
      let c = HpxCell::new(depth_max, self.from >> twice_dd);
      self.from += range_len_min;
      Some(c)
    } else {
      let dd_max_from_len = (H::N_BITS - 1 - len.leading_zeros() as u8) >> 1;
      let dd_max_from_low = (self.from.trailing_zeros() as u8) >> 1;
      let dd = min(H::MAX_DEPTH, min(dd_max_from_len, dd_max_from_low));
      let twice_dd = (dd << 1) as usize;
      let c = HpxCell::new(H::MAX_DEPTH - dd, self.from >> twice_dd);
      self.from += H::one() << twice_dd;
      Some(c)
    }
  }
}

impl<H: HpxHash> Iterator for HpxRange<H> {
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    let len = self.to - self.from;
    if len < H::one() {
      None
    } else {
      let dd_max_from_len = (H::N_BITS - 1 - len.leading_zeros() as u8) >> 1;
      let dd_max_from_low = (self.from.trailing_zeros() as u8) >> 1;
      let dd = min(H::MAX_DEPTH, min(dd_max_from_len, dd_max_from_low));
      let twice_dd = (dd << 1) as usize;
      let c = HpxCell::new(H::MAX_DEPTH - dd, self.from >> twice_dd);
      self.from += H::one() << twice_dd;
      Some(c)
    }
  }

  /* fn size_hint(&self) -> (usize, Option<usize>) { => NEED TO KNOW THE MOC DEPTH!!
    // Worst case: 3 + 3 (at d - 1) + 3 (at d - 2) + ... + 3 (at d - 2) + 3 (at d - 1) + 3
    //             6 * (1 + 4 + 16 + ...)
    // 2^0 + 2^1 + 2^2 + ... + 2^k = 2^(k+1) - 1
    // Solve 6 * n * (2^(k+1) - 1) = range_len with n: range size at the MOC depth
    // =>
    let n = ?
    let dd_max_from_len = H::N_BITS - 1 - (self.end - self.start).leading_zeros();
    let dd_max = max(MAX_DEPTH, dd_max_from_len);
    (1, Some())
  }*/
}

/// Returns the maximum depth a struct implementing this trait contains.
pub trait HasMaxDepth {
  fn depth_max(&self) -> u8;
}

/// Here we define a MOC iterator as a simple sequence of HEALPix cells having the following properties:
/// * the list contains mutually exclusive cells (i.e. no cell is a sub-cell or a super-cell of another cell)
/// * the list is sorted following the natural z-order curve order!
pub trait MOCIterator<H: HpxHash>: Sized + HasMaxDepth + Iterator<Item = HpxCell<H>> {
  fn to_checked(self) -> CheckedIterator<H, Self> {
    CheckedIterator::new(self)
  }

  fn to_range_iter(self) -> RangesFromMOCIterator<H, Self> {
    RangesFromMOCIterator::new(self)
  }

  fn not(self) -> NotMocIter<H, Self> {
    not_unchecked(self)
  }

  fn and<T>(self, rhs: T) -> AndMocIter<H, Self, T>
  where
    T: MOCIterator<H>,
  {
    and_unchecked(self, rhs)
  }

  fn or<T>(self, rhs: T) -> OrMocIter<H, Self, T>
  where
    T: MOCIterator<H>,
  {
    or_unchecked(self, rhs)
  }

  // xor
  // minus

  fn compress(self) -> CompressMocIter<H, Self> {
    compress_unchecked(self)
  }

  // expand(radius) // WILL TAKE MEMORY, CAN'T BE STREAMED!!
}
impl<H: HpxHash, T> MOCIterator<H> for T where T: Sized + HasMaxDepth + Iterator<Item = HpxCell<H>> {}

/// A very basic, unchecked, MOC iterator.
/// We call it lazy because if the stored iterator is the result of a operation,
/// the operation will be performed only as we iterate over the LazyIterator: this
/// is an important point to be aware of when performing benches.
/// For an structure actually performing the operation, see `OwnedMOC`.
pub struct LazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxCell<H>>,
{
  depth_max: u8,
  it: T,
}

impl<H, T> LazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxCell<H>>,
{
  pub fn new(depth_max: u8, it: T) -> LazyMOCIter<H, T> {
    LazyMOCIter { depth_max, it }
  }
}

impl<H, T> HasMaxDepth for LazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxCell<H>>,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<H, T> Iterator for LazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxCell<H>>,
{
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next()
  }
}

pub struct FlatLazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = H>,
{
  depth_max: u8,
  it: T,
}

impl<H, T> FlatLazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = H>,
{
  pub fn new(depth_max: u8, it: T) -> FlatLazyMOCIter<H, T> {
    FlatLazyMOCIter { depth_max, it }
  }
}

impl<H, T> HasMaxDepth for FlatLazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = H>,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<H, T> Iterator for FlatLazyMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = H>,
{
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next().map(|hash| HpxCell {
      depth: self.depth_max,
      hash,
    })
  }
}

/// A very basic, unchecked, MOC.
pub struct OwnedMOC<H>
where
  H: HpxHash,
{
  depth_max: u8,
  moc: Vec<HpxCell<H>>,
}

impl<H> OwnedMOC<H>
where
  H: HpxHash,
{
  pub fn new_unchecked(depth_max: u8, moc: Vec<HpxCell<H>>) -> OwnedMOC<H> {
    OwnedMOC { depth_max, moc }
  }

  /// Internally, a `collect()` is performed.
  pub fn from_it_unchecked<I>(depth_max: u8, it: I) -> OwnedMOC<H>
  where
    I: Iterator<Item = HpxCell<H>>,
  {
    OwnedMOC {
      depth_max,
      moc: it.collect(),
    }
  }

  pub fn into_moc_iter(self) -> LazyMOCIter<H, IntoIter<HpxCell<H>>> {
    LazyMOCIter::new(self.depth_max, self.moc.into_iter())
  }

  pub fn len(&self) -> usize {
    self.moc.len()
  }

  pub fn is_empty(&self) -> bool {
    self.moc.is_empty()
  }
}

type LazyMoc = LazyMOCIter<u64, IntoIter<HpxCell<u64>>>;
type MergeMoc = MergeIter<u64, FlatLazyMOCIter<u64, IntoIter<u64>>>;

impl OwnedMOC<u64> {
  pub fn expand(self) -> OrMocIter<u64, LazyMoc, MergeMoc> {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.moc.len()); // constant to be adjusted
    for HpxCell { depth, hash } in &self.moc {
      append_external_edge(*depth, *hash, self.depth_max - *depth, &mut ext);
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    ext.dedup();
    let depth = self.depth_max;
    self
      .into_moc_iter()
      .or(MergeIter::new(FlatLazyMOCIter::new(depth, ext.into_iter())))
  }

  pub fn external_border(self) -> MinusOnSingleDepthCellsIter<u64, IntoIter<u64>, IntoIter<u64>> {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.moc.len()); // constant to be adjusted
    let mut int: Vec<u64> = Vec::with_capacity(4 * self.moc.len()); // constant to be adjusted
    for HpxCell { depth, hash } in &self.moc {
      append_external_edge(*depth, *hash, self.depth_max - *depth, &mut ext);
      int.append(&mut internal_edge_sorted(
        *depth,
        *hash,
        self.depth_max - *depth,
      ));
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    ext.dedup();
    MinusOnSingleDepthCellsIter::new(self.depth_max, ext.into_iter(), int.into_iter())
  }
}

/// Here we define a MOC as a simple sequence of HEALPix ranges having the following properties:
/// * teh range bounds are expressed at the MAX_DEPTH
/// * the list contains mutually exclusive ranges
/// * the list is sorted from the smaller to the largest bounds  
pub trait RangeMOCIterator<H: HpxHash>: Sized + HasMaxDepth + Iterator<Item = HpxRange<H>> {
  /// Transforms this iterator over ranges into a struct implementing the trait `MOCIterator`.
  fn to_moc_iter(self) -> MOCIteratorFromRanges<H, Self> {
    MOCIteratorFromRanges::new(self)
  }

  // not
  // and
  // or
  // xor
  // minus
  // ...
}
impl<H: HpxHash, T> RangeMOCIterator<H> for T where
  T: Sized + HasMaxDepth + Iterator<Item = HpxRange<H>>
{
}

/// A very basic, unchecked, Range MOC iterator.
/// We call it lazy because if the stored iterator is the result of a operation,
/// the operation will be performed only as we iterate over the LazyIterator: this
/// is an important point to be aware of when performing benches.
/// For an structure actually performing the operation, see `OwnedRangeMOC`.
pub struct LazyRangeMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxRange<H>>,
{
  depth_max: u8,
  it: T,
}

impl<H, T> LazyRangeMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxRange<H>>,
{
  pub fn new(depth_max: u8, it: T) -> LazyRangeMOCIter<H, T> {
    LazyRangeMOCIter { depth_max, it }
  }
}

impl<H, T> HasMaxDepth for LazyRangeMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxRange<H>>,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<H, T> Iterator for LazyRangeMOCIter<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxRange<H>>,
{
  type Item = HpxRange<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.it.next()
  }
}

pub struct LazyRangeMOCVecIter<'a, H>
where
  H: HpxHash,
{
  depth_max: u8,
  iter: Iter<'a, HpxRange<H>>,
}

impl<'a, H> LazyRangeMOCVecIter<'a, H>
where
  H: HpxHash,
{
  pub fn new(depth_max: u8, iter: Iter<'a, HpxRange<H>>) -> LazyRangeMOCVecIter<'a, H> {
    LazyRangeMOCVecIter { depth_max, iter }
  }
}

impl<H> HasMaxDepth for LazyRangeMOCVecIter<'_, H>
where
  H: HpxHash,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

impl<H> Iterator for LazyRangeMOCVecIter<'_, H>
where
  H: HpxHash,
{
  type Item = HpxRange<H>;

  fn next(&mut self) -> Option<Self::Item> {
    self.iter.next().map(|HpxRange { from, to }| HpxRange {
      from: *from,
      to: *to,
    })
  }
}

/// A very basic, unchecked, MOC.
pub struct OwnedRangeMOC<H>
where
  H: HpxHash,
{
  depth_max: u8,
  ranges: Vec<HpxRange<H>>,
}

impl<H> OwnedRangeMOC<H>
where
  H: HpxHash,
{
  pub fn new_unchecked(depth_max: u8, ranges: Vec<HpxRange<H>>) -> OwnedRangeMOC<H> {
    OwnedRangeMOC { depth_max, ranges }
  }

  /// Internally, a `collect()` is performed.
  pub fn from_it<I>(depth_max: u8, it: I) -> OwnedRangeMOC<H>
  where
    I: Iterator<Item = HpxRange<H>>,
  {
    OwnedRangeMOC {
      depth_max,
      ranges: it.collect(),
    }
  }

  pub fn into_range_moc_iter(self) -> LazyRangeMOCIter<H, IntoIter<HpxRange<H>>> {
    LazyRangeMOCIter::new(self.depth_max, self.ranges.into_iter())
  }

  pub fn range_moc_iter(&self) -> LazyRangeMOCVecIter<'_, H> {
    LazyRangeMOCVecIter::new(self.depth_max, self.ranges.iter())
  }
}

type LazyRangeMoc = LazyRangeMOCIter<u64, IntoIter<HpxRange<u64>>>;

impl OwnedRangeMOC<u64> {
  pub fn expand(self) -> OrMocIter<u64, MOCIteratorFromRanges<u64, LazyRangeMoc>, MergeMoc> {
    let mut ext: Vec<u64> = Vec::with_capacity(10 * self.ranges.len()); // constant to be adjusted
    for HpxCell { depth, hash } in MOCIteratorFromRanges::new(self.range_moc_iter()) {
      append_external_edge(depth, hash, self.depth_max - depth, &mut ext);
    }
    ext.sort_unstable(); // parallelize with rayon? It is the slowest part!!
    ext.dedup();
    let depth = self.depth_max;
    MOCIteratorFromRanges::new(self.into_range_moc_iter())
      .or(MergeIter::new(FlatLazyMOCIter::new(depth, ext.into_iter())))
  }
}

impl<H> HasMaxDepth for OwnedRangeMOC<H>
where
  H: HpxHash,
{
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}

/// Transforms a `MOCIterator` into a `RangeMOCIterator`
pub struct RangesFromMOCIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  it: T,
  curr: Option<HpxRange<H>>,
}

impl<H, T> RangesFromMOCIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  fn new(mut it: T) -> RangesFromMOCIterator<H, T> {
    let curr = it.next().map(|c| c.range());
    RangesFromMOCIterator { it, curr }
  }
}

impl<H, T> HasMaxDepth for RangesFromMOCIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl<H, T> Iterator for RangesFromMOCIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  type Item = HpxRange<H>;

  fn next(&mut self) -> Option<Self::Item> {
    if let (Some(r), Some(c)) = (&mut self.curr, &self.it.next()) {
      let nr = c.range();
      if r.to == nr.from {
        r.to = nr.to;
        self.next()
      } else {
        replace(&mut self.curr, Some(nr))
      }
    } else {
      self.curr.take()
    }
  }
}

/// Transforms a `RangeMOCIterator` into a `MOCIterator`.
pub struct MOCIteratorFromRanges<H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  it: R,
  curr: Option<HpxRange<H>>,
  depth_max: u8,
  twice_dd: usize,
  range_len_min: H,
  mask: H,
}

impl<H, R> MOCIteratorFromRanges<H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  fn new(mut it: R) -> MOCIteratorFromRanges<H, R> {
    let curr = it.next();
    let depth_max = it.depth_max();
    let twice_dd = ((H::MAX_DEPTH - depth_max) << 1) as usize;
    let range_len_min = H::one() << twice_dd;
    let mask = (H::one() << 1 | H::one()) << twice_dd;
    MOCIteratorFromRanges {
      it,
      curr,
      depth_max,
      twice_dd,
      range_len_min,
      mask,
    }
  }
}

impl<H, R> HasMaxDepth for MOCIteratorFromRanges<H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl<H, R> Iterator for MOCIteratorFromRanges<H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(c) = &mut self.curr {
      // let res = c.next();
      let res = c.next_with_knowledge(self.depth_max, self.twice_dd, self.range_len_min, self.mask);
      if res.is_none() {
        self.curr = self.it.next();
        self.next()
      } else {
        res
      }
    } else {
      None
    }
  }
}

/// Transforms a ref on `RangeMOCIterator` into a `MOCIterator`.
pub struct MOCIteratorFromRefRanges<'a, H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  it: &'a mut R,
  curr: Option<HpxRange<H>>,
  depth_max: u8,
  twice_dd: usize,
  range_len_min: H,
  mask: H,
}

impl<H, R> MOCIteratorFromRefRanges<'_, H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  pub fn new(it: &mut R) -> MOCIteratorFromRefRanges<'_, H, R> {
    let curr = it.next();
    let depth_max = it.depth_max();
    let twice_dd = ((H::MAX_DEPTH - depth_max) << 1) as usize;
    let range_len_min = H::one() << twice_dd;
    let mask = (H::one() << 1 | H::one()) << twice_dd;
    MOCIteratorFromRefRanges {
      it,
      curr,
      depth_max,
      twice_dd,
      range_len_min,
      mask,
    }
  }
}

impl<H, R> HasMaxDepth for MOCIteratorFromRefRanges<'_, H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl<H, R> Iterator for MOCIteratorFromRefRanges<'_, H, R>
where
  H: HpxHash,
  R: RangeMOCIterator<H>,
{
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    if let Some(c) = &mut self.curr {
      // let res = c.next();
      let res = c.next_with_knowledge(self.depth_max, self.twice_dd, self.range_len_min, self.mask);
      if res.is_none() {
        self.curr = self.it.next();
        self.next()
      } else {
        res
      }
    } else {
      None
    }
  }
}

/// Iterator decorator made to ensure that the decorated iterator returns NESTED HEALPix cells
/// following the natural Z-order curve order.
/// If it is not the case, a call to the `next` method will `panic`!
pub struct CheckedIterator<H, T>
where
  H: HpxHash,
  T: Iterator<Item = HpxCell<H>>,
{
  it: T,
  curr: Option<HpxCell<H>>,
}

impl<H, T> CheckedIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  fn new(mut it: T) -> CheckedIterator<H, T> {
    let curr = it.next();
    CheckedIterator { it, curr }
  }
}

impl<H, T> HasMaxDepth for CheckedIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  fn depth_max(&self) -> u8 {
    self.it.depth_max()
  }
}

impl<H, T> Iterator for CheckedIterator<H, T>
where
  H: HpxHash,
  T: MOCIterator<H>,
{
  type Item = HpxCell<H>;

  fn next(&mut self) -> Option<Self::Item> {
    let prev = replace(&mut self.curr, self.it.next());
    if let (Some(l), Some(r)) = (&prev, &self.curr) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          let hr_at_dl = r.hash.unsigned_shr(((r.depth - l.depth) << 1) as u32);
          assert!(l.hash < hr_at_dl)
        }
        Ordering::Greater => {
          let hl_at_dr = l.hash.unsigned_shr(((l.depth - r.depth) << 1) as u32);
          assert!(hl_at_dr < r.hash)
        }
        Ordering::Equal => {
          debug_assert_eq!(l.depth, r.depth);
          assert!(l.hash < r.hash)
        }
      }
    }
    prev
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let mut sh = self.it.size_hint();
    sh.0 += 1;
    sh.1 = sh.1.map(|v| v + 1);
    sh
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_basic() {
    assert_eq!(u16::N_BYTES, 2);
    assert_eq!(u16::N_BITS, 16);
    assert_eq!(u16::MAX_DEPTH, 5);

    assert_eq!(u32::N_BYTES, 4);
    assert_eq!(u32::N_BITS, 32);
    assert_eq!(u32::MAX_DEPTH, 13);

    assert_eq!(u64::N_BYTES, 8);
    assert_eq!(u64::N_BITS, 64);
    assert_eq!(u64::MAX_DEPTH, 29);

    assert_eq!(u128::N_BYTES, 16);
    assert_eq!(u128::N_BITS, 128);
    assert_eq!(u128::MAX_DEPTH, 61);
  }

  #[test]
  fn test_range() {
    let c = HpxCell::new(1_u8, 10_u64);
    let mut r = c.range();
    assert_eq!(r.next(), Some(c));
    assert_eq!(r.next(), None);
  }
}
