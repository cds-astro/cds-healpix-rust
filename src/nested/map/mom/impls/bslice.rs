use super::super::{
  super::{
    fits::{
      error::FitsError,
      read::{
        check_keyword_and_parse_uint_val, check_keyword_and_str_val, check_keyword_and_val,
        get_str_val_no_quote, next_36_chunks_of_80_bytes, parse_uint_val,
      },
    },
    skymap::{SkyMap, SkyMapValue},
  },
  LhsRhsBoth, Mom, ZUniqHashT,
};
use crate::nested::map::mom::WritableMom;
use chrono::{DateTime, Utc};
use log::debug;
use memmap2::{Mmap, MmapOptions};
use num_traits::{FromBytes, ToBytes};
use std::io::{BufWriter, Write};
use std::{
  array::TryFromSliceError,
  cmp::Ordering,
  convert::{TryFrom, TryInto},
  fs::File,
  iter::{Empty, Map},
  marker::PhantomData,
  ops::Range,
  path::Path,
  slice::{from_raw_parts, ChunksExact},
  str,
  time::SystemTime,
};

/// Defines the type of ZUniq Hash values that can be read/write from/to FITS files.
pub trait Z4FITS:
  ZUniqHashT + ToBytes + FromBytes<Bytes: for<'a> TryFrom<&'a [u8], Error = TryFromSliceError>>
{
}
/// Auto implement it for all types implementing the required traits.
impl<T> Z4FITS for T
where
  T: ZUniqHashT + ToBytes + FromBytes,
  <T as FromBytes>::Bytes: for<'a> TryFrom<&'a [u8], Error = TryFromSliceError>,
{
}

/// Defines the type of values that can be read/write from/to FITS files.
/// Only static lifetime supported here since we build/create the values from the FITS file bytes.
pub trait V4FITS:
  SkyMapValue
  + ToBytes
  + FromBytes<Bytes: for<'a> TryFrom<&'a [u8], Error = TryFromSliceError>>
  + 'static
{
}
/// Auto implement it for all types implementing the required traits.
impl<T> V4FITS for T
where
  T: SkyMapValue + ToBytes + FromBytes + 'static,
  <T as FromBytes>::Bytes: for<'a> TryFrom<&'a [u8], Error = TryFromSliceError>,
{
}

pub enum FITSMom {
  U32U32(FitsMMappedCIndex<u32, u32>),
  U32F32(FitsMMappedCIndex<u32, f32>),
  U32F64(FitsMMappedCIndex<u32, f64>),
  U64U32(FitsMMappedCIndex<u64, u32>),
  U64F32(FitsMMappedCIndex<u64, f32>),
  U64F64(FitsMMappedCIndex<u64, f64>),
}
impl FITSMom {
  // TODO: make a method loading everything from a reader a aking a zvec object!
  #[cfg(not(target_arch = "wasm32"))]
  pub fn from_fits_file<P: AsRef<Path>>(path: P) -> Result<Self, FitsError> {
    let mut file = File::open(path)?;
    let mut raw_header = [b' '; 2880];
    // Parse header
    let mut raw_cards_it = next_36_chunks_of_80_bytes(&mut file, &mut raw_header)?;
    // Parse mandatory, well ordered keywords
    let (n_bytes_per_val, n_rows) =
      check_keyword_and_val(raw_cards_it.next().unwrap(), b"SIMPLE ", b"T")
        .and_then(|()| check_keyword_and_val(raw_cards_it.next().unwrap(), b"BITPIX ", b"8"))
        .and_then(|()| check_keyword_and_val(raw_cards_it.next().unwrap(), b"NAXIS  ", b"2"))
        .and_then(|()| {
          check_keyword_and_parse_uint_val::<u64>(raw_cards_it.next().unwrap(), b"NAXIS1  ")
        })
        .and_then(|n_bytes| {
          check_keyword_and_parse_uint_val::<u64>(raw_cards_it.next().unwrap(), b"NAXIS2  ")
            .map(move |n_rows| (n_bytes, n_rows))
        })?;
    // Parse other keywords
    let mut end_found = false;
    let mut prodtype_found = false;
    let mut dtendian_found = false;
    let mut hpxoder: Option<u8> = None;
    let mut datatype_z: Option<String> = None;
    let mut datatype_v: Option<String> = None;
    let mut value_name: Option<String> = None;
    let mut date: Option<SystemTime> = None;
    for kw_record in &mut raw_cards_it {
      match &kw_record[0..8] {
        b"EXTEND  " => check_keyword_and_val(kw_record, b"EXTEND  ", b"F"),
        b"PRODTYPE" => {
          check_keyword_and_str_val(kw_record, b"PRODTYPE", b"HEALPIX MULTI ORDER MAP")
            .map(|()| prodtype_found = true)
        }
        b"HPXORDER" => parse_uint_val::<u8>(kw_record).map(|v| hpxoder = Some(v)),
        b"ZUNIQ_DT" => get_str_val_no_quote(kw_record)
          .map(|v| datatype_z = Some(String::from_utf8_lossy(v).to_string())),
        b"VALUE_DT" => get_str_val_no_quote(kw_record)
          .map(|v| datatype_v = Some(String::from_utf8_lossy(v).to_string())),
        b"VAL_NAME" => get_str_val_no_quote(kw_record)
          .map(|v| value_name = Some(String::from_utf8_lossy(v).to_string())),
        b"DTENDIAN" => check_keyword_and_str_val(kw_record, b"DTENDIAN", b"LITTLE")
          .map(|()| dtendian_found = true),
        b"DATE    " => get_str_val_no_quote(kw_record).map(|v| {
          date = unsafe { str::from_utf8_unchecked(v) }
            .parse::<DateTime<Utc>>()
            .ok()
            .map(|dt| dt.into())
        }),
        b"CREATOR " => continue,
        b"END     " => {
          end_found = true;
          break;
        }
        _ => {
          debug!("Ignored FITS card: {}", unsafe {
            str::from_utf8_unchecked(kw_record)
          });
          continue;
        }
      }?;
    }
    // Check keywords
    if !end_found {
      return Err(FitsError::new_custom(String::from(
        "'END' keyword not found in the first 36 primary header cards.",
      )));
    }
    if !(prodtype_found & dtendian_found) {
      return Err(FitsError::new_custom(String::from(
        "One of the HEALPIX MULTI ORDER MAP mandatory cards is missing in the FITS header!",
      )));
    }
    let value_name = value_name.ok_or_else(|| {
      FitsError::new_custom(String::from(
        "'VAL_NAME' keyword not found in the primary header cards.",
      ))
    })?;
    let depth = hpxoder.ok_or_else(|| {
      FitsError::new_custom(String::from(
        "'HPXORDER' keyword not found in the primary header cards.",
      ))
    })?;
    let n_bytes_data = n_bytes_per_val * n_rows;
    let mmap = unsafe {
      MmapOptions::new()
        .offset(2880)
        .len(n_bytes_data as usize)
        .map(&file)
        .map_err(FitsError::Io)?
    };
    match (
      datatype_z.as_ref().map(|v| v.as_str()),
      datatype_v.as_ref().map(|v| v.as_str()),
    ) {
      (Some(<u32 as ZUniqHashT>::FITS_DATATYPE), Some(<u32 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U32U32(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(<u32 as ZUniqHashT>::FITS_DATATYPE), Some(<f32 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U32F32(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(<u32 as ZUniqHashT>::FITS_DATATYPE), Some(<f64 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U32F64(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(<u64 as ZUniqHashT>::FITS_DATATYPE), Some(<u32 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U64U32(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(<u64 as ZUniqHashT>::FITS_DATATYPE), Some(<f32 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U64F32(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(<u64 as ZUniqHashT>::FITS_DATATYPE), Some(<f64 as SkyMapValue>::FITS_DATATYPE)) => {
        Ok(FITSMom::U64F64(FitsMMappedCIndex::new(
          date, value_name, depth, n_rows, mmap,
        )))
      }
      (Some(zs), Some(vs)) => Err(FitsError::UnexpectedValue {
        keyword: String::from("ZUNIQ_DT or VALUE_DT"),
        expected: String::from("One of: u32, u64, f32, f64"),
        actual: format!("{} and {}", zs, vs),
      }),
      _ => Err(FitsError::new_custom(String::from(
        "FITS card ZUNIQ_DT or VALUE_DT is missing",
      ))),
    }
  }

  pub fn to_fits_bintable<W: Write>(&self, writer: W) -> Result<(), FitsError> {
    match &self {
      FITSMom::U32U32(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
      FITSMom::U32F32(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
      FITSMom::U32F64(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
      FITSMom::U64U32(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
      FITSMom::U64F32(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
      FITSMom::U64F64(e) => e.get_mom().to_fits_bintable(writer, &e.value_name),
    }
  }

  #[cfg(not(target_arch = "wasm32"))]
  pub fn to_fits_bintable_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FitsError> {
    File::create(path)
      .map_err(FitsError::Io)
      .and_then(|file| self.to_fits_bintable(BufWriter::new(file)))
  }
}

#[derive(Debug)]
pub struct FitsMMappedCIndex<Z: Z4FITS, V: V4FITS> {
  fits_creation_date: Option<SystemTime>,
  value_name: String,
  depth: u8,
  n_rows: u64,
  mmap: Mmap,
  _phantom_z: PhantomData<Z>,
  _phantom_v: PhantomData<V>,
}

impl<Z: Z4FITS, V: V4FITS> FitsMMappedCIndex<Z, V> {
  fn new(
    fits_creation_date: Option<SystemTime>,
    value_name: String,
    depth: u8,
    n_rows: u64,
    mmap: Mmap,
  ) -> Self {
    assert_eq!(
      n_rows * (size_of::<Z>() + size_of::<V>()) as u64,
      mmap.len() as u64
    );
    Self {
      fits_creation_date,
      value_name,
      depth,
      n_rows,
      mmap,
      _phantom_z: PhantomData,
      _phantom_v: PhantomData,
    }
  }

  pub fn get_value_name(&self) -> &str {
    self.value_name.as_str()
  }
  pub fn get_fits_creation_date(&self) -> Option<&SystemTime> {
    self.fits_creation_date.as_ref()
  }

  pub fn get_mom(&self) -> MomSliceImpl<Z, V> {
    MomSliceImpl::new(
      self.depth,
      &self.mmap[0..self.n_rows as usize * (size_of::<Z>() + size_of::<V>())],
    )
  }
  /*
  pub fn into_in_mem_mom(self) -> MomVecImpl<Z, V> {
    self.into()
  }*/
}

/// Implementation of a MOM from a slice of bytes containing ordered `(zuniq, values)` tuples.
/// The purpose of this implementation is to be able to directly read entries from a slice of u8,
/// as stored in a file.   
#[derive(Debug)]
pub struct MomSliceImpl<'b, Z, V>
where
  Z: Z4FITS,
  V: V4FITS,
{
  depth: u8,
  bytes: &'b [u8],
  _z_type: PhantomData<Z>,
  _v_type: PhantomData<V>,
}

impl<'b, Z, V> MomSliceImpl<'b, Z, V>
where
  Z: Z4FITS,
  V: V4FITS,
{
  pub fn new(depth: u8, bytes: &'b [u8]) -> Self {
    Self {
      depth,
      bytes,
      _z_type: PhantomData,
      _v_type: PhantomData,
    }
  }

  pub(crate) fn get_z(&self, index: usize) -> Z {
    let ptr = self.bytes.as_ptr();
    Z::from_le_bytes(
      &unsafe {
        from_raw_parts(
          ptr.add(index * <Self as Mom>::ZV_SIZE),
          <Self as Mom>::Z_SIZE,
        )
      }
      .try_into()
      .unwrap(),
    )
  }

  pub(crate) fn get_entry(&self, index: usize) -> (Z, V) {
    let ptr = self.bytes.as_ptr();
    let sub_slice = unsafe {
      from_raw_parts(
        ptr.add(index * <Self as Mom>::ZV_SIZE),
        <Self as Mom>::ZV_SIZE,
      )
    };
    let (z, v) = sub_slice.split_at(<Self as Mom>::Z_SIZE);
    (
      Z::from_le_bytes(&z.try_into().unwrap()),
      V::from_le_bytes(&v.try_into().unwrap()),
    )
  }

  fn get_entries(&self, range: Range<usize>) -> Map<ChunksExact<u8>, fn(&[u8]) -> (Z, V)> {
    let range = range.start * <Self as Mom>::ZV_SIZE..range.end * <Self as Mom>::ZV_SIZE;
    self.bytes[range]
      .chunks_exact(<Self as Mom>::ZV_SIZE)
      .map(|bytes| {
        (
          Z::from_le_bytes(&bytes[0..<Self as Mom>::Z_SIZE].try_into().unwrap()),
          V::from_le_bytes(&bytes[<Self as Mom>::Z_SIZE..].try_into().unwrap()),
        )
      })
  }

  fn binary_search(&self, zuniq_at_depth_max: Z) -> Result<usize, usize> {
    // Adapted from the binary search algo in std::slice
    let mut size = self.len();
    if size == 0 {
      return Err(0);
    }
    let mut base = 0usize;
    while size > 1 {
      let half = size >> 1;
      let mid = base + half;
      let cmp = self.get_z(mid).cmp(&zuniq_at_depth_max);
      base = if cmp == Ordering::Greater { base } else { mid };
      size -= half;
    }
    let cmp = self.get_z(base).cmp(&zuniq_at_depth_max);
    if cmp == Ordering::Equal {
      Ok(base)
    } else {
      let result = base + (cmp == Ordering::Less) as usize;
      Err(result)
    }
  }
}

impl<'a, Z, V> Mom<'a> for MomSliceImpl<'a, Z, V>
where
  Z: Z4FITS,
  V: V4FITS,
{
  type ZUniqHType = Z;
  type ValueType = V;
  type OverlappedEntries = Empty<(Z, &'a V)>;
  type OverlappedEntriesCopy = Map<ChunksExact<'a, u8>, fn(&'a [u8]) -> (Z, V)>;
  type ZuniqIt = Map<ChunksExact<'a, u8>, fn(&'a [u8]) -> Z>;
  type ValuesIt = Empty<&'a V>;
  type ValuesCopyIt = Map<ChunksExact<'a, u8>, fn(&'a [u8]) -> V>;
  type EntriesIt = Empty<(Z, &'a V)>;
  type EntriesCopyIt = Map<ChunksExact<'a, u8>, fn(&'a [u8]) -> (Z, V)>;
  type OwnedEntriesIt = Self::EntriesCopyIt;

  fn depth_max(&self) -> u8 {
    self.depth
  }

  fn len(&self) -> usize {
    self.bytes.len() / <Self as Mom>::ZV_SIZE
  }

  /// WARNING: not implemented!
  fn get_cell_containing_unsafe(
    &'a self,
    _hash_at_depth_max: Self::ZUniqHType,
  ) -> Option<(Self::ZUniqHType, &'a Self::ValueType)> {
    unimplemented!("Unable to get references while directly reading from bytes. Look at the get_copy_of_cell_containing_unsafe method or the MomVecImpl implementation.")
  }

  fn get_copy_of_cell_containing_unsafe(
    &'a self,
    zuniq_at_depth_max: Self::ZUniqHType,
  ) -> Option<(Self::ZUniqHType, Self::ValueType)> {
    match self.binary_search(zuniq_at_depth_max) {
      Ok(i) => Some(self.get_entry(i)),
      Err(i) => {
        if i > 0 {
          // if array len is 0, i will be 0 so we do not enter here.
          let e = self.get_entry(i - 1);
          if Z::are_overlapping(zuniq_at_depth_max, e.0) {
            return Some(e);
          }
        }
        if i < self.len() {
          let e = self.get_entry(i);
          if Z::are_overlapping(zuniq_at_depth_max, e.0) {
            return Some(e);
          }
        }
        None
      }
    }
  }

  /// WARNING: not implemented!
  fn get_overlapped_cells(&'a self, _zuniq: Self::ZUniqHType) -> Self::OverlappedEntries {
    unimplemented!("Unable to get references while directly reading from bytes. Look at the get_copy_of_overlapped_cells method or the MomVecImpl implementation.")
  }

  fn get_copy_of_overlapped_cells(
    &'a self,
    zuniq: Self::ZUniqHType,
  ) -> Self::OverlappedEntriesCopy {
    let mut range = match self.binary_search(zuniq) {
      Ok(i) => i..i + 1,
      Err(i) => i..i,
    };
    while range.start - 1 > 0 && Z::are_overlapping(zuniq, self.get_z(range.start - 1)) {
      range.start -= 1;
    }
    while range.end < self.len() && Z::are_overlapping(zuniq, self.get_z(range.end)) {
      range.end += 1;
    }
    self.get_entries(range)
  }

  fn zuniqs(&'a self) -> Self::ZuniqIt {
    self
      .bytes
      .chunks_exact(<Self as Mom>::ZV_SIZE)
      .map(|bytes| Z::from_le_bytes(&bytes[0..<Self as Mom>::Z_SIZE].try_into().unwrap()))
  }

  /// WARNING: not implemented!
  fn values(&'a self) -> Self::ValuesIt {
    unimplemented!("Unable to get references while directly reading from bytes. Look at the owned_entries method or the MomVecImpl implementation.")
  }

  fn values_copy(&'a self) -> Self::ValuesCopyIt {
    self
      .bytes
      .chunks_exact(<Self as Mom>::ZV_SIZE)
      .map(|bytes| V::from_le_bytes(&bytes[<Self as Mom>::Z_SIZE..].try_into().unwrap()))
  }

  /// WARNING: not implemented!
  fn entries(&'a self) -> Self::EntriesIt {
    unimplemented!("Unable to get references while directly reading from bytes. Look at the owned_entries method or the MomVecImpl implementation.")
  }

  fn entries_copy(&'a self) -> Self::EntriesCopyIt {
    self
      .bytes
      .chunks_exact(<Self as Mom>::ZV_SIZE)
      .map(|bytes| {
        (
          Z::from_le_bytes(&bytes[0..<Self as Mom>::Z_SIZE].try_into().unwrap()),
          V::from_le_bytes(&bytes[<Self as Mom>::Z_SIZE..].try_into().unwrap()),
        )
      })
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    self
      .bytes
      .chunks_exact(<Self as Mom>::ZV_SIZE)
      .map(|bytes| {
        (
          Z::from_le_bytes(&bytes[0..<Self as Mom>::Z_SIZE].try_into().unwrap()),
          V::from_le_bytes(&bytes[<Self as Mom>::Z_SIZE..].try_into().unwrap()),
        )
      })
  }

  /// WARNING: not implemented, see `MomVecImpl`!
  fn from_skymap_ref<'s, S, M>(_skymap: &'s S, _merger: M) -> Self
  where
    S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
    M: Fn(u8, Self::ZUniqHType, [&Self::ValueType; 4]) -> Option<Self::ValueType>,
    Self::ValueType: 's,
  {
    unimplemented!("Unable to create this object from a skymap ref. Look at MomVecImpl instead.")
  }

  /// WARNING: not implemented, see `MomVecImpl`!
  fn from_skymap<'s, S, M>(_skymap: S, _merger: M) -> Self
  where
    S: SkyMap<'s, HashType = Self::ZUniqHType, ValueType = Self::ValueType>,
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>,
    Self::ValueType: 's,
  {
    unimplemented!("Unable to create this object from a skymap. Look at MomVecImpl instead.")
  }

  /// WARNING: not implemented, see `MomVecImpl`!
  fn merge<'s, L, R, S, O, M>(_lhs: L, _rhs: R, _split: S, _op: O, _merge: M) -> Self
  where
    L: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
    R: Mom<'s, ZUniqHType = Self::ZUniqHType, ValueType = Self::ValueType>,
    S: Fn(u8, Self::ZUniqHType, Self::ValueType) -> [Self::ValueType; 4],
    O: Fn(LhsRhsBoth<Self::ValueType>) -> Option<Self::ValueType>,
    M: Fn(
      u8,
      Self::ZUniqHType,
      [Self::ValueType; 4],
    ) -> Result<Self::ValueType, [Self::ValueType; 4]>,
  {
    unimplemented!(
      "Unable to create this object from the merge operation. Look at MomVecImpl instead."
    )
  }
}
