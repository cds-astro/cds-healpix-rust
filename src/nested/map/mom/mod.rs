//! Multi-order Map, i.e. list of non-overlapping `(UNIQ, VALUE)`.
// IDEE:
// * use zuniq
// * use array of tuple, sorted according to zuniq
// * for each img pixel, get sky position, get zuniq, binay_search in map

use std::{
  cmp::Ordering,
  fmt::{Debug, Display},
  mem,
};
use std::ops::AddAssign;

use num::PrimInt;

use super::{
  HHash,
  skymap::{SkyMap, SkyMapValue}
};

pub mod impls;

/// `ZUniqHHashT` stands for `Z-curve ordered Uniq Healpix Hash Type`.
pub trait ZUniqHashT:
// 'static mean that Idx does not contains any reference
'static + HHash + PrimInt + Send + Sync + Debug + AddAssign + Display + Clone
{
  /// Number of unused bits (without including the sentinel bit).
  const N_RESERVED_BITS: u8 = 2;
  /// Number of bits needed to code a sub-cell relative index.
  const DIM: u8 = 2;
  /// Number of base cells, i.e. number of cell at depth 0.
  const N_D0_CELLS: u8 = 12;
  /// Number of bits needed to code the base cell index.
  const N_D0_BITS: u8 = 4;
  /// Number of bytes available in the 'ZUniqHHashT'.
  const N_BYTES: u8 = mem::size_of::<Self>() as u8;
  /// Number of bits available in the 'zuniq'.
  const N_BITS: u8 = Self::N_BYTES << 3;
  /// Maximum depth that can be coded on this 'ZUniqHHashT'.
  const MAX_DEPTH: u8 = (Self::N_BITS - (Self::N_RESERVED_BITS + Self::N_D0_BITS)) / Self::DIM;
  /// MASK to retrieve the bits corresponding to the deepest layer.
  /// Must equals 3.
  const LAST_LAYER_MASK: Self; // = Self:: Self::one().unsigned_shl(Self::DIM as u32) - Self::one();

  /// Must return 2.
  fn two() -> Self;
  /// Must return 3.
  fn three() -> Self;
  fn mult_by_dim<T: PrimInt>(v: T) -> T {
    v.unsigned_shl(1)
  }
  fn div_by_dim<T: PrimInt>(v: T) -> T {
    v.unsigned_shr(1)
  }
  fn shift(delta_depth: u8) -> u8 {
    Self::mult_by_dim(delta_depth)
  }
  fn delta_with_depth_max(depth: u8) -> u8 {
    Self::MAX_DEPTH - depth
  }
  fn shift_from_depth_max(depth: u8) -> u8 {
    Self::shift(Self::delta_with_depth_max(depth))
  }
  
  /// Returns the `depth` and `hash value` this `zuniq` contains.
  fn from_zuniq(zuniq: Self) -> (u8, Self) {
    let n_trailing_zero = zuniq.trailing_zeros() as u8;
    let delta_depth = Self::div_by_dim(n_trailing_zero);
    let depth = Self::MAX_DEPTH - delta_depth;
    let hash = zuniq >> (n_trailing_zero + 1) as usize;
    (depth, hash)
  }

  fn depth_from_zuniq(zuniq: Self) -> u8 {
    let n_trailing_zero = zuniq.trailing_zeros() as u8;
    let delta_depth = Self::div_by_dim(n_trailing_zero);
    Self::MAX_DEPTH - delta_depth
  }

  /// Transforms a `depth` and `hash value` tuple into a `zuniq`.
  fn to_zuniq(depth: u8, hash: Self) -> Self {
    let zuniq = (hash << 1) | Self::one();
    zuniq.unsigned_shl(Self::shift_from_depth_max(depth) as u32)
  }

  fn are_overlapping(zuniq_l: Self, zuniq_r: Self) -> bool {
    let (depth_l, hash_l) = Self::from_zuniq(zuniq_l);
    let (depth_r, hash_r) = Self::from_zuniq(zuniq_r);
    Self::are_overlapping_cells(depth_l, hash_l, depth_r, hash_r)
  }
  
  fn are_overlapping_cells(depth_l: u8, hash_l: Self, depth_r: u8, hash_r: Self) -> bool {
    match depth_l.cmp(&depth_r) {
      Ordering::Equal   => hash_l == hash_r,
      Ordering::Less    => hash_l == hash_r.unsigned_shr(((depth_r - depth_l) << 1) as u32),
      Ordering::Greater => hash_r == hash_l.unsigned_shr(((depth_l - depth_r) << 1) as u32),
    }
  }

}

impl ZUniqHashT for u32 {
  const LAST_LAYER_MASK: Self = 3;
  fn two() -> Self {
    2
  }
  fn three() -> Self {
    Self::LAST_LAYER_MASK
  }
}
impl ZUniqHashT for u64 {
  const LAST_LAYER_MASK: Self = 3;
  fn two() -> Self {
    2
  }
  fn three() -> Self {
    Self::LAST_LAYER_MASK
  }
}

/// `MOM` stands for **M**ulti **O**rder healpix **M**aps.
/// Here, it consists in a list of HEALPix cells (hash) at various depth (HEALPixordes) 
/// with a value attached to each cell.
/// All cells in a given MOM:
/// * must be non-overlapping, 
/// * have the same type of value attached to them
/// * are sorted following the z-order curve order.
/// This last property allow for streaming processing and operations.
/// In practice, we use the `zuniq` index: it encodes both the depth and the cell hash value
/// and is built in such a way that the natural order follow the z-order curve order.  
pub trait Mom<'a>: Sized {
  /// Type of the HEALPix zuniq hash value (mainly `u32` or `u64`).
  type ZUniqHType: ZUniqHashT;
  /// Type of the iterator iterating on the skymap values.
  type ValueType: SkyMapValue + 'a;
  /// Type of iterator iterating on all (sorted!) zuniq values.
  type ZuniqIt: Iterator<Item = Self::ZUniqHType>;
  /// Type of the iterator iterating on the MOM values.
  type ValuesIt: Iterator<Item = &'a Self::ValueType>;
  /// Type of iterator iterating on all (sorted!) entries.
  /// # Remark
  /// We could have defined `Iterator<Item = &'a (Self::ZUniqHType, Self::ValueType)>;`
  /// but it would have limited the possible implementations (e.g. using 2 vectors and the `zip` operation.
  type EntriesIt: Iterator<Item = (Self::ZUniqHType, &'a Self::ValueType)>;
  /// Type of iterator iterating on all (sorted!) owned entries.
  type OwnedEntriesIt: Iterator<Item = (Self::ZUniqHType, Self::ValueType)>;

  /// Largest depth the MOM may contain.
  fn depth_max(&self) -> u8;

  /// Returns the number of elements in the `mom`.
  fn len(&self) -> usize;
  
  /// Returns the entry, if any, containing the given HEALPix cell hash computed at the `Mom`
  /// maximum depth.
  fn get_cell_containing(&'a self, zuniq_at_depth_max: Self::ZUniqHType) -> Result<Option<(Self::ZUniqHType, &'a Self::ValueType)>, String> {
    self.check_zuniq_depth_is_depth_max(zuniq_at_depth_max)
      .map(|_| self.get_cell_containing_unsafe(zuniq_at_depth_max))
  }

  /// Same as `get_cell_containing` without checking that `zuniq_at_depth_max` depth is the MOM
  /// maximum depth.
  fn get_cell_containing_unsafe(&'a self, hash_at_depth_max: Self::ZUniqHType) -> Option<(Self::ZUniqHType, &'a Self::ValueType)>;

  /// Returns all entries overlapped by the HEALPix cell of given `zuniq` hash value.
  fn get_overlapped_cells(&'a self, zuniq: Self::ZUniqHType) -> Vec<(Self::ZUniqHType, &'a Self::ValueType)>;

  /// Returns all HEALPix zuniq hash, ordered following the z-order curve.
  fn zuniqs(&'a self) -> Self::ZuniqIt;

  /// Returns all values associated with HEALPix cells, ordered by increasing cell hash number.
  fn values(&'a self) -> Self::ValuesIt;
  
  /// Returns all entries, i.e. HEALPix zuniq hash / value tuples, ordered following the z-order curve.
  fn entries(&'a self) -> Self::EntriesIt;
  
  fn owned_entries(self) -> Self::OwnedEntriesIt;

  /// Check if the gieven `zuniq` depth is the MOM maximum depth.
  fn check_zuniq_depth_is_depth_max(&self, zuniq_at_depth_max: Self::ZUniqHType) -> Result<(), String> {
    let depth = Self::ZUniqHType::depth_from_zuniq(zuniq_at_depth_max);
    if depth == self.depth_max() {
      Ok(())
    } else {
      Err(format!("Wrong depth for zuniq {}. Expected: {}. Actual: {}", zuniq_at_depth_max, self.depth_max(), depth))
    }
  }

  /// # Warning
  /// * assumes that the order of the elements returned by `zuniqs()` is the same as the one returned
  ///   by `entries()`.
  fn check_is_mom(&'a self) -> Result<(), String> {
    let mut it = self.zuniqs();
    if let Some(mut l) = it.next() {
      let (mut depth_l, mut hash_l) = Self::ZUniqHType::from_zuniq(l);
      for r in it {
        if depth_l < self.depth_max() {
          return Err(format!("Element has a larger depth than MOM maximum depth. Elem: {}; Depth: {}; Mom max depth: {}", l, depth_l, self.depth_max()));
        }
        let (depth_r, hash_r) =  Self::ZUniqHType::from_zuniq(r);
        if l >= r {
          return Err(format!("The MOM is not ordered: {} >= {}", l, r));
        } else if Self::ZUniqHType::are_overlapping_cells(depth_l, hash_l, depth_r, hash_r) {
          return Err(format!("Overlapping elements in the MOM: {} and {}. I.e. depth: {}; hash: {} and depth: {}, hash: {}.", l, r, depth_l, hash_l, depth_r, hash_r));
        }
        l = r;
        depth_l = depth_r;
        hash_l = hash_r;
      }
    }
    Ok(())
  }

  /// # Params
  /// * `M`: merger function, i.e. function applied on the 4 values of 4 sibling cells
  /// (i.e. the 4 cells belonging to a same direct parent cell).
  /// The function decide whether value are merge (and how they are merged) or not returning
  /// either `Some` or `None`.
  fn from_skymap_ref<'s, S, M>(skymap: &'s S, merger: M) -> Self
    where
      S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
      M: Fn(&Self::ValueType, &Self::ValueType, &Self::ValueType, &Self::ValueType) -> Option<Self::ValueType>,
      Self::ValueType: 's;

  /// # Params
  /// * `M`: returns `Ok` if merge succeed, else return `Err` with original elements (in the same order).
  fn from_skymap<'s, S, M>(skymap: S, merger: M) -> Self
    where
      S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
      M: Fn(Self::ValueType, Self::ValueType, Self::ValueType, Self::ValueType) -> Result<Self::ValueType, (Self::ValueType, Self::ValueType, Self::ValueType, Self::ValueType)>,
      Self::ValueType: 's;
  
  fn merge<S,O,M>(self, rhs: Self, split: S, op: O, merge: M) -> Self
    where
      // Split a parent cell into for siblings
      S: Fn(Self::ValueType) -> (Self::ValueType, Self::ValueType, Self::ValueType, Self::ValueType),
      // Merge the values of a same cell from both MOMs
      O: Fn(Self::ValueType, Self::ValueType) -> Self::ValuesIt,
      // Performs a post-merge operation?
      M: Fn(Self::ValueType, Self::ValueType, Self::ValueType, Self::ValueType) -> Self::ValueType,
  {
    // get both owned iterators
    todo!()
  }

  
}



/*

pub struct Mom {
  pub depth_max: u8,
  pub elems: MomElems,
}

pub enum FitsMomElems {
  U64U8(Vec<(u64, u8)>),
  U64I16(Vec<(u64, i16)>),
  U64I32(Vec<(u64, i32)>),
  U64I64(Vec<(u64, i64)>),
  U64F32(Vec<(u64, f32)>),
  U64F64(Vec<(u64, f64)>),
}

*/
