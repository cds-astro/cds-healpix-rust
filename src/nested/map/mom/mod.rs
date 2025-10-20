//! Multi-Order healpix Map, i.e. list of non-overlapping `(UNIQ, VALUE)`.
//! We could also call this a `Valued-MOC` (`VMOC`), i.e. a list of non-overlapping HEALPix cells, at
//! various oder, having a value (or a quantity) associated to each cell.

use std::{
  cmp::Ordering,
  error::Error,
  fmt::{Debug, Display},
  fs::File,
  io::{BufWriter, Write},
  ops::{AddAssign, RangeInclusive},
  path::Path,
};

use chrono::{SecondsFormat, Utc};
use colorous::Gradient;
use mapproj::CanonicalProjection;
use num::{Integer, PrimInt};
use num_traits::{AsPrimitive, Float, FloatConst, FromPrimitive, ToBytes};

#[cfg(not(target_arch = "wasm32"))]
use super::img::to_mom_png_file;
use super::{
  skymap::{SkyMap, SkyMapValue},
  HHash,
};
use crate::nested::{
  map::{
    fits::{
      error::FitsError,
      write::{
        write_final_padding, write_keyword_record, write_primary_hdu, write_str_keyword_record,
        write_str_mandatory_keyword_record, write_uint_mandatory_keyword_record,
      },
    },
    img::{to_mom_png, ColorMapFunctionType, PosConversion, Val},
  },
  n_hash,
};

pub mod impls;

/// `ZUniqHHashT` stands for `Z-curve ordered Uniq Healpix Hash Type`.
pub trait ZUniqHashT:
// 'static mean that Idx does not contains any reference
'static + HHash + PrimInt + Send + Sync + Debug + AddAssign + Display + Clone + ToBytes
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
  const N_BYTES: u8 = size_of::<Self>() as u8;
  /// Number of bits available in the 'zuniq'.
  const N_BITS: u8 = Self::N_BYTES << 3;
  /// Maximum depth that can be coded on this 'ZUniqHHashT'.
  const MAX_DEPTH: u8 = (Self::N_BITS - (Self::N_RESERVED_BITS + Self::N_D0_BITS)) / Self::DIM;
  /// MASK to retrieve the bits corresponding to the deepest layer.
  /// Must equals 3.
  const LAST_LAYER_MASK: Self; // = Self:: Self::one().unsigned_shl(Self::DIM as u32) - Self::one();
  /// The datatype writen in FITS
  const FITS_DATATYPE: &'static str;
  // The TFORM compliant FITS keyword value.
  // const FITS_TFORM: &'static str;

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
      Ordering::Equal => hash_l == hash_r,
      Ordering::Less => hash_l == hash_r.unsigned_shr(((depth_r - depth_l) << 1) as u32),
      Ordering::Greater => hash_r == hash_l.unsigned_shr(((depth_l - depth_r) << 1) as u32),
    }
  }
}

impl ZUniqHashT for u32 {
  const LAST_LAYER_MASK: Self = 3;
  const FITS_DATATYPE: &'static str = "u32";
  // const FITS_TFORM: &'static str = "J";
  fn two() -> Self {
    2
  }
  fn three() -> Self {
    Self::LAST_LAYER_MASK
  }
}
impl ZUniqHashT for u64 {
  const LAST_LAYER_MASK: Self = 3;
  const FITS_DATATYPE: &'static str = "u64";
  // const FITS_TFORM: &'static str = "K";
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
#[allow(clippy::len_without_is_empty, clippy::type_complexity)]
pub trait Mom<'a>: Sized {
  /// Byte size, in memory, of the zuniq type.
  const Z_SIZE: usize = size_of::<Self::ZUniqHType>();
  /// Byte size, in memory, of the value type.
  const V_SIZE: usize = size_of::<Self::ValueType>();
  /// Byte size, in memory, of the zuniq plus the value type.
  const ZV_SIZE: usize = Self::Z_SIZE + Self::V_SIZE;
  /// Type of the HEALPix zuniq hash value (mainly `u32` or `u64`).
  type ZUniqHType: ZUniqHashT;
  /// Type of the iterator iterating on the skymap values.
  type ValueType: SkyMapValue + 'a;
  /// Type of the iterator iterating over a part of the entries.
  type OverlappedEntries: Iterator<Item = (Self::ZUniqHType, &'a Self::ValueType)>;
  /// Type of the iterator iterating over a part of the entries in which the values are copied.
  type OverlappedEntriesCopy: Iterator<Item = (Self::ZUniqHType, Self::ValueType)>;
  /// Type of iterator iterating on all (sorted!) zuniq values.
  type ZuniqIt: Iterator<Item = Self::ZUniqHType>;
  /// Type of the iterator iterating on the MOM values.
  type ValuesIt: Iterator<Item = &'a Self::ValueType>;
  /// Owned values, probably performs a clone so avoid using it on non-copy values.
  type ValuesCopyIt: Iterator<Item = Self::ValueType>;
  /// Type of iterator iterating on all (sorted!) entries.
  /// # Remark
  /// We could have defined `Iterator<Item = &'a (Self::ZUniqHType, Self::ValueType)>;`
  /// but it would have limited the possible implementations (e.g. using 2 vectors and the `zip` operation.
  type EntriesIt: Iterator<Item = (Self::ZUniqHType, &'a Self::ValueType)>;
  /// Type of iterator iterating on all (sorted!) copied entries.
  type EntriesCopyIt: Iterator<Item = (Self::ZUniqHType, Self::ValueType)>;
  /// Type of iterator iterating on all (sorted!) owned entries.
  type OwnedEntriesIt: Iterator<Item = (Self::ZUniqHType, Self::ValueType)>;

  /// Largest depth the MOM may contain.
  fn depth_max(&self) -> u8;

  /// Returns the number of elements in the `mom`.
  fn len(&self) -> usize;

  /// Returns the entry, if any, containing the given HEALPix cell hash computed at the `Mom`
  /// maximum depth.
  fn get_cell_containing(
    &'a self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Result<Option<(Self::ZUniqHType, &'a Self::ValueType)>, String> {
    self
      .check_zuniq_depth_is_depth_max(zuniq_at_depth_max)
      .map(|_| self.get_cell_containing_unsafe(zuniq_at_depth_max))
  }

  /// Returns the entry, if any, containing the given HEALPix cell hash computed at the `Mom`
  /// maximum depth, and making a copy of the value.
  fn get_copy_of_cell_containing(
    &'a self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Result<Option<(Self::ZUniqHType, Self::ValueType)>, String> {
    self
      .check_zuniq_depth_is_depth_max(zuniq_at_depth_max)
      .map(|_| self.get_copy_of_cell_containing_unsafe(zuniq_at_depth_max))
  }

  /// Same as `get_cell_containing` without checking that `zuniq_at_depth_max` depth is the MOM
  /// maximum depth.
  fn get_cell_containing_unsafe(
    &'a self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Option<(Self::ZUniqHType, &'a Self::ValueType)>;

  /// Same as `get_copy_of_cell_containing` without checking that `zuniq_at_depth_max` depth is the MOM
  /// maximum depth.
  fn get_copy_of_cell_containing_unsafe(
    &'a self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Option<(Self::ZUniqHType, Self::ValueType)>;

  /// Returns all entries, in the z-order curve order, overlapped by the HEALPix cell of given `zuniq` hash value.
  fn get_overlapped_cells(&'a self, zuniq: Self::ZUniqHType) -> Self::OverlappedEntries;

  /// Returns all entries, in the z-order curve order, overlapped by the HEALPix cell of given `zuniq` hash value,
  /// making a copy of the value.
  fn get_copy_of_overlapped_cells(&'a self, zuniq: Self::ZUniqHType)
    -> Self::OverlappedEntriesCopy;

  /// Returns all HEALPix zuniq hash, ordered following the z-order curve.
  fn zuniqs(&'a self) -> Self::ZuniqIt;

  /// Returns all values associated with HEALPix cells, ordered by increasing cell hash number.
  fn values(&'a self) -> Self::ValuesIt;

  /// Returns a copy of the values associated with HEALPix cells, ordered by increasing cell hash number.
  /// Avoid using this method on non-copy values.
  fn values_copy(&'a self) -> Self::ValuesCopyIt;

  /// Returns all entries, i.e. HEALPix zuniq hash / value tuples, ordered following the z-order curve.
  fn entries(&'a self) -> Self::EntriesIt;

  /// Returns a copy of all entries, i.e. HEALPix zuniq hash / value tuples, ordered following the z-order curve.
  fn entries_copy(&'a self) -> Self::EntriesCopyIt;

  /// Returns owned entries, i.e. HEALPix zuniq hash / value tuples, ordered following the z-order curve.
  fn owned_entries(self) -> Self::OwnedEntriesIt;

  /// Check if the given `zuniq` depth is the MOM maximum depth.
  fn check_zuniq_depth_is_depth_max(
    &self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Result<(), String> {
    let depth = Self::ZUniqHType::depth_from_zuniq(zuniq_at_depth_max);
    if depth == self.depth_max() {
      Ok(())
    } else {
      Err(format!(
        "Wrong depth for zuniq {}. Expected: {}. Actual: {}",
        zuniq_at_depth_max,
        self.depth_max(),
        depth
      ))
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
        let (depth_r, hash_r) = Self::ZUniqHType::from_zuniq(r);
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
  /// * `depth`; the HEALPix depth of the hash values provide as first parameter of the given iterator tuples.
  /// * `sorted_entries_it`: an iterator over `(hash, values)` tuples. The iterator **must** be sorted
  /// * `merger`: returns `Ok` if merge succeed, else return `Err` with original elements (in the same order).
  /// # Warning
  /// * the caller has to ensure that the iterator is sorted. In case of doubt, create a decorator
  /// similar to the ones in `nested::iter::checksorted` and call `from_hpx_sorted_entries_fallible` instead.
  fn from_hpx_sorted_entries<I, M>(depth: u8, sorted_entries_it: I, merger: M) -> Self
  where
    I: Iterator<Item = (Self::ZUniqHType, Self::ValueType)>,
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>;

  /// # Params
  /// * `depth`; the HEALPix depth of the hash values provide as first parameter of the given iterator tuples.
  /// * `sorted_entries_it`: a fallible iterator over `(hash, values)` tuples. The iterator **must** be sorted
  /// * `merger`: returns `Ok` if merge succeed, else return `Err` with original elements (in the same order).
  /// # Warning
  /// * the caller has to ensure that the iterator is sorted. In case of doubt, create a decorator
  /// similar to the ones in `nested::iter::checksorted`.
  fn from_hpx_sorted_entries_fallible<E, I, M>(
    depth: u8,
    sorted_entries_it: I,
    merger: M,
  ) -> Result<Self, E>
  where
    E: Error,
    I: Iterator<Item = Result<(Self::ZUniqHType, Self::ValueType), E>>,
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>;

  /// # Params
  /// * `M`: merger function, i.e. function applied on the 4 values of 4 sibling cells
  ///   (i.e. the 4 cells belonging to a same direct parent cell).
  ///   The function decide whether values are merge (and how they are merged) or not returning
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
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>,
    Self::ValueType: 's;

  /// Merge two input MOMs to produce an output MOM.
  /// # Params
  /// * `lhs`: the Left Hand Side MOM
  /// * `rhs`: the right Hand side MOM
  /// * `split`: the operator splitting a lhs or rhs cell into its four direct siblings
  ///   + the depth and hash value of the parent cell are provided in input.
  ///   + together with the value of the cell to be split
  /// * `op`: merge/mixe/join the values of a same cell from both input MOMs (or decide what to do
  ///   when one of the MOM cell has a value and the other does not).
  ///   + it allows for various types of JOIN (inner, left, right, full).
  /// * `merge`: possibly performs a post-merge operation, considering the (four) siblings of the join operation.
  ///   + the depth and hash value of the parent cell (if cells are merged) are provided in input,
  ///   + together with the 4 values to be possibly merged in the parent cell.
  ///   + result if either `Ok` with the merged value if a merge occurs, or `Err` with the input values if
  ///     the merge does not occurs.
  /// * `L`: type of the left MOM
  /// * `R`: type of the right MOM
  fn merge<'s, L, R, S, O, M>(lhs: L, rhs: R, split: S, op: O, merge: M) -> Self
  where
    L: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
    R: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
    S: Fn(u8, Self::ZUniqHType, Self::ValueType) -> [Self::ValueType; 4],
    O: Fn(LhsRhsBoth<Self::ValueType>) -> Option<Self::ValueType>,
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>;
}

pub trait ViewableMom<'a>: Mom<'a>
where
  <Self as Mom<'a>>::ValueType: SkyMapValue + Val + 'a,
{
  /// # Params
  /// * `size`: the `(X, Y)` number of pixels in the image;
  /// * `proj`: a projection, if different from Mollweide;
  /// * `proj_center`: the `(lon, lat)` coordinates of the center of the projection, in radians,
  ///                      if different from `(0, 0)`;
  /// * `proj_bounds`: the `(X, Y)` bounds of the projection, if different from the default values
  ///                  which depends on the projection. For unbounded projections, de default value
  ///                  is `(-PI..PI, -PI..PI)`.
  /// * `path`: the path of th PNG file to be written.
  /// * `view`: set to true to visualize the saved image.
  #[cfg(not(target_arch = "wasm32"))]
  fn to_mom_png_file<P: CanonicalProjection, W: AsRef<Path>>(
    &'a self,
    img_size: (u16, u16),
    proj: Option<P>,
    proj_center: Option<(f64, f64)>,
    proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
    pos_convert: Option<PosConversion>,
    color_map: Option<Gradient>,
    color_map_func_type: Option<ColorMapFunctionType>,
    path: W,
    view: bool,
  ) -> Result<(), Box<dyn Error>> {
    to_mom_png_file(
      self,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
      path,
      view,
    )
  }

  fn to_mom_png<P: CanonicalProjection, W: Write>(
    &'a self,
    img_size: (u16, u16),
    proj: Option<P>,
    proj_center: Option<(f64, f64)>,
    proj_bounds: Option<(RangeInclusive<f64>, RangeInclusive<f64>)>,
    pos_convert: Option<PosConversion>,
    color_map: Option<Gradient>,
    color_map_func_type: Option<ColorMapFunctionType>,
    writer: W,
  ) -> Result<(), Box<dyn Error>> {
    to_mom_png(
      self,
      img_size,
      proj,
      proj_center,
      proj_bounds,
      pos_convert,
      color_map,
      color_map_func_type,
      writer,
    )
  }
}
/// Implement `ViewableMom` for all `Mom` having a `ValueType` implementing `Val`.
impl<'a, V: SkyMapValue + Val + 'a, M: Mom<'a, ValueType = V>> ViewableMom<'a> for M {}

pub trait WritableMom<'a>: Mom<'a>
where
  <Self as Mom<'a>>::ValueType: SkyMapValue + ToBytes + 'a,
{
  fn write_all_entries<W: Write>(&'a self, mut writer: W) -> Result<usize, FitsError> {
    // TODO: make a method trying first 'entries', adn 'entries_copy' only if 'entries' no implemented
    for (z, v) in self.entries_copy() {
      writer
        .write_all(z.to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(v.to_le_bytes().as_ref()))
        .map_err(FitsError::Io)?;
    }
    Ok(self.len() * Self::ZV_SIZE)
  }

  fn write_all_bintable_entries<W: Write>(&'a self, mut writer: W) -> Result<usize, FitsError> {
    // TODO: make a method trying first 'entries', adn 'entries_copy' only if 'entries' no implemented
    for (z, v) in self.entries_copy() {
      // let (d, h) = Self::ZUniqHType::from_zuniq(z);
      writer
        // .write_u8(d)
        .write_all(z.to_be_bytes().as_ref())
        .and_then(|()| writer.write_all(v.to_be_bytes().as_ref()))
        .map_err(FitsError::Io)?;
    }
    Ok(self.len() * Self::ZV_SIZE)
  }

  #[cfg(not(target_arch = "wasm32"))]
  fn to_fits_file<P: AsRef<Path>, T: AsRef<str>>(
    &'a self,
    path: P,
    value_name: T,
  ) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits(BufWriter::new(file), value_name))
  }

  /// # Params
  /// * `value_name`: the name of the value so we know the king of quantity this map stores.
  fn to_fits<T: AsRef<str>, W: Write>(
    &'a self,
    mut writer: W,
    value_name: T,
  ) -> Result<(), FitsError> {
    // Prepare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS1  ", Self::ZV_SIZE as u64); // Len of data axis 1 = number of bytes in zuniq + in values
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.len() as u64); // Len of data axis 2 = n_hash(depth)+1
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..35].copy_from_slice(b"PRODTYPE= 'HEALPIX MULTI ORDER MAP'"); // Product type
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth_max() as u64);
    write_str_keyword_record(
      it.next().unwrap(),
      b"ZUNIQ_DT",
      Self::ZUniqHType::FITS_DATATYPE,
    );
    write_str_keyword_record(
      it.next().unwrap(),
      b"VALUE_DT",
      Self::ValueType::FITS_DATATYPE,
    );
    write_str_keyword_record(it.next().unwrap(), b"VAL_NAME", value_name.as_ref());
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    write_keyword_record(
      it.next().unwrap(),
      b"DATE    ",
      format!(
        "'{}'",
        Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true)
      )
      .as_str(),
    );
    write_keyword_record(
      it.next().unwrap(),
      b"CREATOR ",
      format!(
        "'Rust crate {} {}'",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
      )
      .as_str(),
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..]).map_err(FitsError::Io)?;
    // Write the data part
    self
      .write_all_entries(&mut writer)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }

  #[cfg(not(target_arch = "wasm32"))]
  fn to_fits_bintable_file<P: AsRef<Path>, T: AsRef<str>>(
    &'a self,
    path: P,
    value_name: T,
  ) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits_bintable(BufWriter::new(file), value_name))
  }

  /// Write the MOM in a bintable with 3 columns:
  /// * healpix order (depth)
  /// * healpix cell number (hash value)
  /// * the associated value
  ///
  /// This is supposedly less efficient to work with (to be tested), and required extra space
  /// (`depth` coded on a `u8` instead of a sentinel bit in the cell number, useless ITS PRimary header
  /// of 2880 bytes), but has a broder compatibility since it can be read from a all tools supporting
  /// FITS BINTABLE reading.
  // Allow for other, user custom, FIST cards in input parameters?
  fn to_fits_bintable<T: AsRef<str>, W: Write>(
    &'a self,
    mut writer: W,
    value_colname: T,
  ) -> Result<(), FitsError> {
    self
      .write_bintable_fits_header(&mut writer, value_colname)
      .and_then(|()| self.write_all_bintable_entries(&mut writer))
      .and_then(|n_data_bytes| write_final_padding(&mut writer, n_data_bytes))
  }

  fn write_bintable_fits_header<T: AsRef<str>, W: Write>(
    &'a self,
    mut writer: W,
    value_colname: T,
  ) -> Result<(), FitsError> {
    write_primary_hdu(&mut writer)?;
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    // Write BINTABLE specific keywords in the buffer
    it.next().unwrap()[0..20].copy_from_slice(b"XTENSION= 'BINTABLE'");
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8");
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS1  ", Self::ZV_SIZE as u64);
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.len() as u64);
    it.next().unwrap()[0..30].copy_from_slice(b"PCOUNT  =                    0");
    it.next().unwrap()[0..30].copy_from_slice(b"GCOUNT  =                    1");
    it.next().unwrap()[0..30].copy_from_slice(b"TFIELDS =                    2");
    it.next().unwrap()[0..20].copy_from_slice(b"TTYPE1  = 'ZUNIQ   '");
    write_str_mandatory_keyword_record(
      it.next().unwrap(),
      b"TFORM1  ",
      Self::ZUniqHType::FITS_TFORM,
    );
    it.next().unwrap()[0..20].copy_from_slice(b"TTYPE2  = 'VALUE   '");
    write_str_mandatory_keyword_record(
      it.next().unwrap(),
      b"TFORM2  ",
      Self::ValueType::FITS_TFORM,
    );
    it.next().unwrap()[0..23].copy_from_slice(b"PRODTYPE= 'HEALPIX MOM'");
    it.next().unwrap()[0..20].copy_from_slice(b"COORDSYS= 'CEL     '");
    it.next().unwrap()[0..20].copy_from_slice(b"EXTNAME = 'xtension'");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"MAXORDER", self.depth_max() as u64);
    write_str_mandatory_keyword_record(it.next().unwrap(), b"VALNAME ", value_colname.as_ref());
    write_keyword_record(
      it.next().unwrap(),
      b"CREATOR ",
      format!(
        "'Rust crate {} {}'",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
      )
      .as_str(),
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..]).map_err(FitsError::Io)
  }
}
/// Implement `WritableMom` for all `Mom` having a `ValueType` implementing `ToBytes`.
impl<'a, V: SkyMapValue + ToBytes + 'a, M: Mom<'a, ValueType = V>> WritableMom<'a> for M {}

/// Provide a merger merging density values based on a chi-square criterion.
/// # Params
/// * `chi2_of_3dof_threshold`: threshold on the value of the chi square distribution with 3
/// degrees of freedom below which we consider the 4 values of 4 sibling cells as coming
/// from the same normal distribution which mean and variance comes from a poisson distribution.
/// Here a few typical values corresponding the the given completeness:
/// * Completeness = 90.0% =>  6.251
/// * Completeness = 95.0% =>  7.815
/// * Completeness = 97.5% =>  9.348
/// * Completeness = 99.0% => 11.345
/// * Completeness = 99.9% => 16.266
pub fn new_chi2_density_merger<
  Z: ZUniqHashT,
  V: SkyMapValue + Float + FloatConst + FromPrimitive,
>(
  chi2_of_3dof_threshold: f64,
) -> impl Fn(u8, Z, [V; 4]) -> Result<V, [V; 4]> {
  let merger_ref = new_chi2_density_ref_merger(chi2_of_3dof_threshold);
  move |depth: u8, hash: Z, values: [V; 4]| -> Result<V, [V; 4]> {
    merger_ref(
      depth,
      hash,
      [&values[0], &values[1], &values[2], &values[3]],
    )
    .ok_or(values)
  }
}

pub fn new_chi2_density_ref_merger<
  Z: ZUniqHashT,
  V: SkyMapValue + Float + FloatConst + FromPrimitive,
>(
  chi2_of_3dof_threshold: f64,
) -> impl Fn(u8, Z, [&V; 4]) -> Option<V> {
  let threshold = V::from_f64(chi2_of_3dof_threshold).unwrap();
  let four = V::from_f64(4.0).unwrap();
  move |depth: u8, _hash: Z, [n0, n1, n2, n3]: [&V; 4]| -> Option<V> {
    // With Poisson distribution:
    // * s_i = Surface of cell i = s (all cell have the same surface at a given depth)
    // * mu_i = Number of source in cell i / Surface of cell i = Density in cell i
    // * sigma_i = sqrt(Number of source in cell i) / Surface cell i = sqrt(mu_i / s)
    // weight_i = 1 / sigma_i^2 = s / mu_i
    // mu_e = weighted_mean = ( sum_{1=1}^4 weight_i * mu_i ) / ( sum_{1=1}^4 weight_i )
    //                      = 4 / ( sum_{1=1}^4 1/mu_i )
    // V_e^{-1} = s * sum_{1=1}^4 1/mu_i
    // Applying Pineau et al. 2017:
    // => sum_{1=1}^4 (mu_i - mu_e)^2 / sigma_i^2 = ... = s * [(sum_{1=1}^4 mu_i) - 4 * mu_e]
    // let four = V::from_f64(4.0).unwrap();
    let one_over_s = V::from_u64(n_hash(depth + 1).unsigned_shr(2)).unwrap() / V::PI();

    let mu0 = *n0;
    let mu1 = *n1;
    let mu2 = *n2;
    let mu3 = *n3;

    let sum = mu0 + mu1 + mu2 + mu3;
    let weighted_var_inv = V::one() / mu0.max(one_over_s)
      + V::one() / mu1.max(one_over_s)
      + V::one() / mu2.max(one_over_s)
      + V::one() / mu3.max(one_over_s);
    // let weighted_var_inv = 1.0 / mu0 + 1.0 / mu1 + 1.0 / mu2 + 1.0 / mu3;
    let weighted_mean = four / weighted_var_inv;
    let chi2_of_3dof = (sum - four * weighted_mean) / one_over_s;
    if chi2_of_3dof < threshold {
      Some(sum / four)
    } else {
      None
    }
  }
}

/// Provide a merger merging counts based on a chi-square criterion.
/// # Params
/// * `chi2_of_3dof_threshold`: threshold on the value of the chi square distribution with 3
/// degrees of freedom below which we consider the 4 values of 4 sibling cells as coming
/// from the same normal distribution which mean and variance comes from a poisson distribution.
/// Here a few typical values corresponding the the given completeness:
/// * Completeness = 90.0% =>  6.251
/// * Completeness = 95.0% =>  7.815
/// * Completeness = 97.5% =>  9.348
/// * Completeness = 99.0% => 11.345
/// * Completeness = 99.9% => 16.266
pub fn new_chi2_count_merger<Z: ZUniqHashT, V: SkyMapValue + Integer + AsPrimitive<f64>>(
  chi2_of_3dof_threshold: f64,
) -> impl Fn(u8, Z, [V; 4]) -> Result<V, [V; 4]> {
  let merger_ref = new_chi2_count_ref_merger(chi2_of_3dof_threshold);
  move |depth: u8, hash: Z, values: [V; 4]| -> Result<V, [V; 4]> {
    merger_ref(
      depth,
      hash,
      [&values[0], &values[1], &values[2], &values[3]],
    )
    .ok_or(values)
  }
}

pub fn new_chi2_count_ref_merger<Z: ZUniqHashT, V: SkyMapValue + Integer + AsPrimitive<f64>>(
  chi2_of_3dof_threshold: f64,
) -> impl Fn(u8, Z, [&V; 4]) -> Option<V> {
  move |_depth: u8, _hash: Z, [n0, n1, n2, n3]: [&V; 4]| -> Option<V> {
    let mu0 = (*n0).as_();
    let mu1 = (*n1).as_();
    let mu2 = (*n2).as_();
    let mu3 = (*n3).as_();

    let sum = mu0 + mu1 + mu2 + mu3;
    let weighted_var_inv = 1.0 / mu0 + 1.0 / mu1 + 1.0 / mu2 + 1.0 / mu3;
    let weighted_mean = 4.0 / weighted_var_inv;
    let chi2_of_3dof = sum - 4.0 * weighted_mean;

    if chi2_of_3dof < chi2_of_3dof_threshold {
      Some(*n0 + *n1 + *n2 + *n3)
    } else {
      None
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::nested::map::img::to_mom_png_file;
  use crate::nested::map::mom::{new_chi2_count_ref_merger, LhsRhsBoth};
  use mapproj::pseudocyl::mol::Mol;

  use super::{
    super::{
      img::{ColorMapFunctionType, PosConversion},
      skymap::SkyMapEnum,
    },
    impls::zvec::MomVecImpl,
    Mom,
  };

  #[test]
  #[cfg(not(target_arch = "wasm32"))]
  fn test_mom_diff_spec() {
    let path = "test/resources/skymap/skymap.2mass.depth6.fits";
    let skymap_2mass = SkyMapEnum::from_fits_file(path).unwrap();
    let path = "test/resources/skymap/skymap.gaiadr3.depth6.fits";
    let skymap_gaia = SkyMapEnum::from_fits_file(path).unwrap();
    match (skymap_2mass, skymap_gaia) {
      (SkyMapEnum::ImplicitU64I32(iskymap_2mass), SkyMapEnum::ImplicitU64I32(iskymap_gaia_dr3)) => {
        let chi2_merger = new_chi2_count_ref_merger(16.266);
        let mom_2mass = MomVecImpl::from_skymap_ref(&iskymap_2mass, chi2_merger);
        let chi2_merger = new_chi2_count_ref_merger(16.266);
        let mom_gaia = MomVecImpl::from_skymap_ref(&iskymap_gaia_dr3, chi2_merger);
        // Transform count MOMs into PDF MOMs.
        let n_tot_2mass = mom_2mass.values().sum::<i32>() as f64;
        let n_tot_gaia = mom_gaia.values().sum::<i32>() as f64;
        let mom_2mass = MomVecImpl::from_map(mom_2mass, |_, n| n as f64 / n_tot_2mass);
        let mom_gaia = MomVecImpl::from_map(mom_gaia, |_, n| n as f64 / n_tot_gaia);

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
        let op = |lrb: LhsRhsBoth<f64>| match lrb {
          LhsRhsBoth::Left(l) => Some(l),
          LhsRhsBoth::Right(r) => Some(r),
          LhsRhsBoth::Both(l, r) => Some(l - r),
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
      _ => panic!(),
    }
  }
}
