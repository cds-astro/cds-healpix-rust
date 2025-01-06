//! Multi-order Map, i.e. list of non-overlapping `(UNIQ, VALUE)`.
// IDEE:
// * use zuniq
// * use array of tuple, sorted according to zuniq
// * for each img pixel, get sky position, get zuniq, binay_search in map

use std::{
  cmp::Ordering,
  fmt::{Debug, Display},
  mem,
  ops::AddAssign,
};

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

/// Enum defining the 3 possibilities whe merging two nullable inputs:
/// * Left: non-null LHS, null RHS 
/// * Right: null LHS, non-null RHS
/// * Both: non-null LHS, non-null RHS
pub enum LhsRhsBoth<V> {
  Left(V),
  Right(V),
  Both(V, V),
}

/// `MOM` stands for **M**ulti **O**rder healpix **M**aps.
/// Here, it consists in a list of HEALPix cells (hash) at various depth (HEALPixordes) 
/// with a value attached to each cell.
/// All cells in a given MOM:
/// * must be non-overlapping, 
/// * have the same type of value attached to them
/// * are sorted following the z-order curve order.
/// 
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
  ///   (i.e. the 4 cells belonging to a same direct parent cell).
  ///   The function decide whether value are merge (and how they are merged) or not returning
  ///   either `Some` or `None`.
  fn from_skymap_ref<'s, S, M>(skymap: &'s S, merger: M) -> Self
    where
      S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
      M: Fn(u8, Self::ZUniqHType, [&Self::ValueType; 4]) -> Option<Self::ValueType>,
      Self::ValueType: 's;

  /// # Params
  /// * `M`: returns `Ok` if merge succeed, else return `Err` with original elements (in the same order).
  fn from_skymap<'s, S, M>(skymap: S, merger: M) -> Self
    where
      S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
      M: Fn(u8, Self::ZUniqHType, [Self::ValueType; 4]) -> Result<Self::ValueType, [Self::ValueType; 4]>,
      Self::ValueType: 's;
  
  /// Merge two input MOMs to produce an output MOM.
  /// # Params
  /// * `lhs`: the Left Hand Side MOM
  /// * `rhs`: the right Hand side MOM
  /// * `split`: the operator splitting a lhs or rhs cell into its four direct siblings
  ///     + the depth and hash value of the parent cell are provided in input.
  ///     + together with the value of the cell to be split
  /// * `op`: merge/mixe/join the values of a same cell from both input MOMs (or decide what to do
  ///         when one of the MOM cell has a value and the other does not).
  ///     + it allows for various types of JOIN (inner, left, right, full).
  /// * `merge`: possibly performs a post-merge operation, considering the (four) siblings of the join operation.
  ///     + the depth and hash value of the parent cell (if cells are merged) are provided in input,
  ///     + together with the 4 values to be possibly merged in the parent cell.
  ///     + result if either `Ok` with the merged value if a merge occurs, or `Err` with the input values if 
  ///       the merge does not occurs.
  /// * `L`: type of the left MOM
  /// * `R`: type of the right MOM
  fn merge<'s, L, R, S, O, M>(lhs: L, rhs: R, split: S, op: O, merge: M) -> Self
    where
      L: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
      R: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
      S: Fn(u8, Self::ZUniqHType, Self::ValueType) -> [Self::ValueType; 4],
      O: Fn(LhsRhsBoth<Self::ValueType>) -> Option<Self::ValueType>,
      M: Fn(u8, Self::ZUniqHType, [Self::ValueType; 4]) -> Result<Self::ValueType, [Self::ValueType; 4]>;

  
}


#[cfg(test)]
mod tests {
  use mapproj::pseudocyl::mol::Mol;
  use crate::nested::map::img::to_mom_png_file;
  use crate::nested::map::mom::LhsRhsBoth;

  use super::{
    Mom,
    impls::zvec::MomVecImpl,
    super::{
      img::{ColorMapFunctionType, PosConversion},
      skymap::SkyMapEnum,
    }
  };

  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_mom_diff_spec() {
    let path = "test/resources/skymap/skymap.2mass.depth6.fits";
    let skymap_2mass = SkyMapEnum::from_fits_file(path).unwrap();
    let path = "test/resources/skymap/skymap.gaiadr3.depth6.fits";
    let skymap_gaia = SkyMapEnum::from_fits_file(path).unwrap();

    let chi2_merger = |_depth: u8, _hash: u64, [n0, n1, n2, n3]: [&i32; 4]| -> Option<i32> {
      let mu0 = *n0 as f64;
      // let sig0 = mu0.sqrt();
      let mu1 = *n1 as f64;
      // let sig1 = mu1.sqrt();
      let mu2 = *n2 as f64;
      // let sig2 = mu2.sqrt();
      let mu3 = *n3 as f64;
      // let sig3 = mu3.sqrt();

      let sum = mu0 + mu1 + mu2 + mu3;
      let weighted_var_inv = 1.0 / mu0 + 1.0 / mu1 + 1.0 / mu2 + 1.0 / mu3;
      let weighted_mean = 4.0 / weighted_var_inv;
      let chi2_of_3dof = sum - 4.0 * weighted_mean;
      // chi2 3 dof:
      // 90.0% =>  6.251
      // 95.0% =>  7.815
      // 97.5% =>  9.348
      // 99.0% => 11.345
      // 99.9% => 16.266
      if chi2_of_3dof < 16.266 {
        Some(*n0 + *n1 + *n2 + *n3)
      } else {
        None
      }
    };

    match (skymap_2mass, skymap_gaia) {
      (SkyMapEnum::ImplicitU64I32(iskymap_2mass), SkyMapEnum::ImplicitU64I32(iskymap_gaia_dr3)) => {
        let mom_2mass = MomVecImpl::from_skymap_ref(&iskymap_2mass, chi2_merger);
        let mom_gaia = MomVecImpl::from_skymap_ref(&iskymap_gaia_dr3, chi2_merger);
        // Transform count MOMs into PDF MOMs.
        let n_tot_2mass = mom_2mass.values().sum::<i32>() as f64;
        let n_tot_gaia = mom_gaia.values().sum::<i32>() as f64;
        let mom_2mass = MomVecImpl::from(mom_2mass, |_, n| n as f64 / n_tot_2mass);
        let mom_gaia = MomVecImpl::from(mom_gaia, |_, n| n as f64 / n_tot_gaia);

        to_mom_png_file::<'_, _, Mol, _>(
          &mom_2mass,
          (1600, 800),
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/mom_2mass.png",
          false,
        )
          .unwrap();
        to_mom_png_file::<'_, _, Mol, _>(
          &mom_gaia,
          (1600, 800),
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::LinearLog), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/mom_gaia.png",
          false,
        )
          .unwrap();
        
        
        /*let split= |n: i32| {
          let n_div_4 = n >> 2;
          // Preserve the total value by possibly adding 1 to the 3 last cells
          let mut a = (n_div_4 & 1 == 1);
          let b = (n_div_4 & 2 == 2);
          let c = (n_div_4 & 3 == 3);
          a |= b;
          let a = a as i32;
          let b = b as i32;
          let c = c as i32;
          assert!(a == 0 || a == 1);
          assert!(b == 0 || b == 1);
          assert!(c == 0 || c == 1);
          (n_div_4, n_div_4 + a, n_div_4 + b, n_div_4 + c)
        };
        */
        let split = |_depth: u8, _hash: u64, pdf: f64| {
          let pdf_div_4 = pdf / 4.0;
          [pdf_div_4, pdf_div_4, pdf_div_4, pdf_div_4]
        };
        let op = |lrb: LhsRhsBoth<f64>| {
          match lrb {
            LhsRhsBoth::Left(l) => Some(l),
            LhsRhsBoth::Right(r) => Some(r),
            LhsRhsBoth::Both(l, r) => Some(l - r),
          }
        };
        let merge = |_depth: u8, _hash: u64, [v1, v2, v3, v4]: [f64; 4]| {
          let sum = v1 + v2 + v3 + v4;
          if sum.abs() < 1e-4 {
            Ok(sum)
          } else {
            Err([v1, v2, v3, v4])
          }
        };
        let mom = MomVecImpl::merge(mom_2mass, mom_gaia, split, op, merge);
        to_mom_png_file::<'_, _, Mol, _>(
          &mom,
          (1600, 800),
          None,
          None,
          None,
          Some(PosConversion::EqMap2GalImg),
          None,
          Some(ColorMapFunctionType::Linear), //Some(ColorMapFunctionType::LinearSqrt)
          "test/resources/skymap/mreged_mom.png",
          false,
        )
          .unwrap();
      }
      _ => assert!(false)
    }
  }
}