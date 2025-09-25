//! Defines a HEALPix Sampled Cumulative Index for a HEALPix sorted files.
//!
//! This index is well-balanced (contrary to the simple HEALPix Cumulative Index) but requires
//! binary search and has fuzzy borders (i.e. hash value are not on low level cell boundaries).
//!
//! Here, we sample the sorted file and simply stored, each `N` (the sample parameter) rows,
//! the value of the HEALPix index of the row at a given -- constant -- depth,
//! together with value of the starting row byte in the file.
//!
//! For formats of fixed row byte length, and knowing the value `N` used to build the index,
//! it is not necessary to store the starting byte.
//!
//! In practice, `N` in often chosen according to the desired size of the index in bytes
//! (or, equivalently, its number of elements).
//!
//! `N`, is the sampling size or parameter.
//! The index stores the HEALPix values of rows at row indices: 0, N, 2N, 3N, ...
//! The last value is the value of the last row in the file.
//!
//! # Remark:
//!A better approach could consist in to store such an index in
//! a [bstree-file-readonly](https://github.com/cds-astro/cds-bstree-file-readonly-rust) file since
//! the latter is organized to minimize query time (avoiding distant memory access in the first
//! part of the binary-search process).

// hpx type: u32 / u64
// value type: cumul bytes or implicit (index * sampling_step)

use std::{
  cmp::Ordering,
  fs::File,
  io::{Error as IoError, Write},
  marker::PhantomData,
  ops::Range,
  path::Path,
  ptr::slice_from_raw_parts,
  time::SystemTime,
};

use chrono::{DateTime, SecondsFormat, Utc};
use log::debug;
use memmap2::{Mmap, MmapOptions};

use crate::{
  depth_from_n_hash_unsafe, n_hash,
  nested::{
    children,
    map::fits::{
      error::FitsError,
      read::{
        check_keyword_and_parse_uint_val, check_keyword_and_str_val, check_keyword_and_val,
        get_str_val_no_quote, next_36_chunks_of_80_bytes, parse_uint_val,
      },
      write::{
        write_final_padding, write_keyword_record, write_str_keyword_record,
        write_uint_mandatory_keyword_record,
      },
    },
    parent,
  },
};

/// Possible datatype of the HEALPix index value used as Key in the index.
pub trait UInt: 'static + Sized + Ord + PartialOrd + PartialEq + Eq {
  const FITS_DATATYPE_KEY: &'static str;
  fn to_u64(self) -> u64;
  fn to_u64_range(range: Range<Self>) -> Range<u64>;
  fn from_u64(val: u64) -> Self;
  fn from_u64_range(range: Range<u64>) -> Range<Self>;
}

impl UInt for u32 {
  const FITS_DATATYPE_KEY: &'static str = "u32";
  fn to_u64(self) -> u64 {
    self as u64
  }
  fn to_u64_range(range: Range<Self>) -> Range<u64> {
    let Range { start, end } = range;
    start.to_u64()..end.to_u64()
  }
  fn from_u64(val: u64) -> Self {
    val as u32
  }
  fn from_u64_range(range: Range<u64>) -> Range<Self> {
    let Range { start, end } = range;
    UInt::from_u64(start)..UInt::from_u64(end)
  }
}

impl UInt for u64 {
  const FITS_DATATYPE_KEY: &'static str = "u64";
  fn to_u64(self) -> u64 {
    self
  }
  fn to_u64_range(range: Range<Self>) -> Range<u64> {
    range
  }
  fn from_u64(val: u64) -> Self {
    val
  }
  fn from_u64_range(range: Range<u64>) -> Range<Self> {
    range
  }
}

/// Defines common methods for all types of Sampled HEALPix Cumulative Index (SHCI).
pub trait SHCIndex<T: UInt> {
  /// HEALPix depth of the keys of the index.
  fn depth(&self) -> u8;
  /// Total number of indexed rows, i.e. number fo rows in the indexed file.
  fn n_rows(&self) -> u64;
  /// Number of rows between two index values.
  fn sampling_step(&self) -> u64;
  /// Number of elements in the index (basically `n_rows/sampling_size`, `+1` if `n_rows % sampling_size != 0`).
  fn len(&self) -> usize;

  /// Returns the index of the first row indexed by the SHCI element `i`.
  /// In practice, returns  `i * n` if `i` is lower than the index `len`, else returns `n_rows`.
  fn starting_row(&self, i: usize) -> u64 {
    if i < self.len() {
      (i as u64) * self.sampling_step()
    } else {
      self.n_rows()
    }
  }

  /// Returns the starting and ending index of rows indexed by the SHCI elements in the given range.
  fn row_range(&self, Range { start, end }: Range<usize>) -> Range<u64> {
    self.starting_row(start)..self.starting_row(end)
  }

  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth(&self, hash: u64) -> Range<usize>;

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<usize>;

  /// For the fits serialization. Must write the `self.len()` key, in ascending order,
  /// and in little-endian, in the given writer.
  /// Exactly `self.len() * size_of::<T>` bytes must be written.
  fn write_all_keys<W: Write>(&self, writer: W) -> Result<usize, IoError>;
}

/// In an `Implicit` SHCIndex, we return records numbers without having to store values.
/// For the last index, we have to store in the metadata the total number of rows.
pub trait ImplicitSHCIndex<T: UInt>: SHCIndex<T> {
  /// The input hash value in the range are at the index depth.
  fn get_recnos_containing_hash_at_index_depth(&self, hash: u64) -> Range<u64> {
    self.row_range(self.get_indices_of_elems_containing_hash_at_index_depth(hash))
  }

  /// The input hash value in the range are at the index depth.
  fn get_recnos_containing_hash_range_at_index_depth(&self, hash_range: Range<u64>) -> Range<u64> {
    self.row_range(self.get_indices_of_elems_containing_hash_range_at_index_depth(hash_range))
  }

  fn get_recnos_containing_hash(&self, depth: u8, hash: u64) -> Range<u64> {
    match depth.cmp(&self.depth()) {
      Ordering::Less => {
        self.get_recnos_containing_hash_range_at_index_depth(children(hash, self.depth() - depth))
      }
      Ordering::Equal => self.get_recnos_containing_hash_at_index_depth(hash),
      Ordering::Greater => {
        self.get_recnos_containing_hash_at_index_depth(parent(hash, depth - self.depth()))
      }
    }
  }

  fn get_recnos_containing_hash_range(&self, depth: u8, hash_range: Range<u64>) -> Range<u64> {
    match depth.cmp(&self.depth()) {
      Ordering::Less => {
        let twice_dd = (self.depth() - depth) << 1;
        let start = hash_range.start << twice_dd;
        let end = hash_range.start << twice_dd;
        self.get_recnos_containing_hash_range_at_index_depth(start..end)
      }
      Ordering::Equal => self.get_recnos_containing_hash_range_at_index_depth(hash_range),
      Ordering::Greater => {
        let twice_dd = (depth - self.depth()) << 1;
        let start = hash_range.start >> twice_dd;
        let end = hash_range.start >> twice_dd;
        if start == end {
          self.get_recnos_containing_hash_at_index_depth(start)
        } else {
          self.get_recnos_containing_hash_range_at_index_depth(start..end)
        }
      }
    }
  }

  fn to_fits<W: Write>(
    &self,
    mut writer: W,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    let indexed_file_last_modif_date = indexed_file_last_modif_date.map(DateTime::<Utc>::from);
    // Perpare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS1  ", size_of::<T>() as u64); // Len of data axis 1 = number of bytes in dataype
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.len() as u64);
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..31].copy_from_slice(b"PRODTYPE= 'HEALPIX SIMPLE CUMUL INDEX'"); // Product type
    it.next().unwrap()[0..20].copy_from_slice(b"VAL_SCHM= 'IMPLICIT'");
    it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth() as u64); // Must be NSIDE in case of RING
    write_str_keyword_record(it.next().unwrap(), b"DATATYPE", T::FITS_DATATYPE_KEY);
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NTOTROWS", self.n_rows());
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"SAMPSTEP", self.sampling_step());
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    if let Some(indexed_file_name) = indexed_file_name {
      write_str_keyword_record(it.next().unwrap(), b"IDXF_NAM", indexed_file_name);
    }
    if let Some(indexed_file_len) = indexed_file_len {
      write_keyword_record(
        it.next().unwrap(),
        b"IDXF_LEN",
        indexed_file_len.to_string().as_str(),
      );
    }
    if let Some(indexed_file_mdr5) = indexed_file_md5 {
      write_str_keyword_record(it.next().unwrap(), b"IDXF_MD5", unsafe {
        str::from_utf8_unchecked(&indexed_file_mdr5)
      });
    }
    if let Some(indexed_file_last_modif_date) = indexed_file_last_modif_date {
      write_str_keyword_record(
        it.next().unwrap(),
        b"IDXF_LMD",
        indexed_file_last_modif_date
          .to_rfc3339_opts(SecondsFormat::Secs, true)
          .as_str(),
      );
    }
    if let Some(indexed_colname_lon) = indexed_colname_lon {
      write_str_keyword_record(it.next().unwrap(), b"IDXC_LON", indexed_colname_lon);
    }
    if let Some(indexed_colname_lat) = indexed_colname_lat {
      write_str_keyword_record(it.next().unwrap(), b"IDXC_LAT", indexed_colname_lat);
    }
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
      .write_all_keys(&mut writer)
      .map_err(FitsError::Io)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }
}

/// In an `Explicit` SHCIndex, we return pointers -- or byte offsets -- stored as values
/// associated to each Key of the index.
pub trait ExplicitSHCIndex<T: UInt>: SHCIndex<T> {
  /// Returns the value associated to key of given index `i`.
  fn value(&self, i: usize) -> u64;

  fn get_pointers_containing_hash_at_index_depth(&self, hash: u64) -> Range<u64> {
    let Range { start, end } = self.get_indices_of_elems_containing_hash_at_index_depth(hash);
    self.value(start)..self.value(end)
  }

  /// The input hash value in the range are at the index depth.
  fn get_pointers_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<u64> {
    let Range { start, end } =
      self.get_indices_of_elems_containing_hash_range_at_index_depth(hash_range);
    self.value(start)..self.value(end)
  }

  /// For the fits serialization. Must write the values (byte pointers), in ascending order,
  /// and in little-endian, in the given writer.
  /// Exactly `self.len() * size_of::<u64>` bytes must be written.
  fn write_all_values<W: Write>(&self, writer: W) -> Result<usize, IoError>;

  fn to_fits<W: Write>(
    &self,
    mut writer: W,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    let indexed_file_last_modif_date = indexed_file_last_modif_date.map(DateTime::<Utc>::from);
    // Perpare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    write_uint_mandatory_keyword_record(
      it.next().unwrap(),
      b"NAXIS1  ",
      (size_of::<T>() + size_of::<u64>()) as u64,
    ); // Len of data axis 1 = number of bytes in dataype
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.len() as u64);
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..31].copy_from_slice(b"PRODTYPE= 'HEALPIX SIMPLE CUMUL INDEX'"); // Product type
    it.next().unwrap()[0..20].copy_from_slice(b"VAL_SCHM= 'EXPLICIT'");
    it.next().unwrap()[0..24].copy_from_slice(b"DATASCHM= 'COL_ORIENTED'"); // MANDATORY. Means first keys, then values
    it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth() as u64); // Must be NSIDE in case of RING
    write_str_keyword_record(it.next().unwrap(), b"DATATYPE", T::FITS_DATATYPE_KEY);
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NTOTROWS", self.n_rows());
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"SAMPSTEP", self.sampling_step());
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    if let Some(indexed_file_name) = indexed_file_name {
      write_str_keyword_record(it.next().unwrap(), b"IDXF_NAM", indexed_file_name);
    }
    if let Some(indexed_file_len) = indexed_file_len {
      write_keyword_record(
        it.next().unwrap(),
        b"IDXF_LEN",
        indexed_file_len.to_string().as_str(),
      );
    }
    if let Some(indexed_file_mdr5) = indexed_file_md5 {
      write_str_keyword_record(it.next().unwrap(), b"IDXF_MD5", unsafe {
        str::from_utf8_unchecked(&indexed_file_mdr5)
      });
    }
    if let Some(indexed_file_last_modif_date) = indexed_file_last_modif_date {
      write_str_keyword_record(
        it.next().unwrap(),
        b"IDXF_LMD",
        indexed_file_last_modif_date
          .to_rfc3339_opts(SecondsFormat::Secs, true)
          .as_str(),
      );
    }
    if let Some(indexed_colname_lon) = indexed_colname_lon {
      write_str_keyword_record(it.next().unwrap(), b"IDXC_LON", indexed_colname_lon);
    }
    if let Some(indexed_colname_lat) = indexed_colname_lat {
      write_str_keyword_record(it.next().unwrap(), b"IDXC_LAT", indexed_colname_lat);
    }
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
      .write_all_keys(&mut writer)
      .and_then(|n_bytes_written| {
        self
          .write_all_values(&mut writer)
          .map(|n_val_bytes_written| n_bytes_written + n_val_bytes_written)
      })
      .map_err(FitsError::Io)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }
}

/// In an `implicit` SHCIndex, the key is a row number which is computed from the index of
/// the value (knowing `n`, the sample size). For the last index, we have to store in the metadata
/// the total number of rows.
pub struct OwnedImplicitSHCIndex<T: UInt> {
  /// Sampling size, i.e. number of rows between two successive indexed values.
  /// Remark: could be a u32, but u64 chosen to avoid cast in various operations.
  n: u64,
  /// Total number of rows in the original files.
  n_rows: u64,
  /// Healpix depth of the index (in practice, 13 for u32 and 29 for u64).
  depth: u8,
  /// Sampled HEALPix hash values
  sampled_hash: Box<[T]>,
}
impl<T: UInt> OwnedImplicitSHCIndex<T> {
  // from_iter_of_positions ...
}
impl<T: UInt> SHCIndex<T> for OwnedImplicitSHCIndex<T> {
  fn depth(&self) -> u8 {
    self.depth
  }
  fn n_rows(&self) -> u64 {
    self.n_rows
  }
  fn sampling_step(&self) -> u64 {
    self.n
  }
  fn len(&self) -> usize {
    self.sampled_hash.len()
  }

  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth(&self, hash: u64) -> Range<usize> {
    BorrowedImplicitSHCIndex::get_indices_of_elems_containing_hash_at_index_depth_gen(
      &self.sampled_hash,
      hash,
    )
  }

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<usize> {
    BorrowedImplicitSHCIndex::get_indices_of_elems_containing_hash_range_at_index_depth_gen(
      &self.sampled_hash,
      hash_range,
    )
  }

  fn write_all_keys<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    /*for k in self.sampled_hash.iter() {
      writer.write_all(k.as_le_bytes())?
    }*/
    let pointer = self.sampled_hash.as_ref();
    let n_bytes = pointer.len() * size_of::<T>();
    let bytes = unsafe { std::slice::from_raw_parts(pointer.as_ptr() as *const u8, n_bytes) };
    writer.write_all(bytes).map(|()| n_bytes)
  }
}
impl<T: UInt> ImplicitSHCIndex<T> for OwnedImplicitSHCIndex<T> {}

/// The hash value at index `i` is the value of the row at index `i * n`.
/// If `n_rows % n = 0`, the index contains `n_rows / n` elements.
/// Else it contains `(n_rows / n) + 1` elements.
/// The last element of the index contains the hash value of the last rows.
/// So that:
/// * `hash[0]` is the smaller possible value in all rows.
/// * `hash[(n_rows / n) + 1]` or `hash[(n_rows / n)]` is the largest possible value in all rows.
pub struct BorrowedImplicitSHCIndex<'a, T: UInt> {
  /// Sampling size, i.e. number of rows between two successive indexed values.
  /// Remark: could be a u32, but u64 chosen to avoid cast in various operations.
  n: u64,
  /// Total number of rows in the original files.
  n_rows: u64,
  /// Healpix depth of the index (in practice, 13 for u32 and 29 for u64).
  depth: u8,
  /// Sampled HEALPix hash values
  sampled_hash: &'a [T],
}

impl<'a, T: UInt> BorrowedImplicitSHCIndex<'a, T> {
  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth_gen(
    sampled_hash: &[T],
    hash: u64,
  ) -> Range<usize> {
    let len = sampled_hash.len();
    let hash = UInt::from_u64(hash);
    match sampled_hash.binary_search(&hash) {
      Ok(mut i) => {
        let mut j = i;
        // Explore the left side
        if i > 0 {
          i -= 1;
          while i > 0 && sampled_hash[i] == hash {
            i -= 1;
          }
        }
        // Explore the right side
        if j < len {
          j += 1;
          while j < len && sampled_hash[j] == hash {
            j += 1;
          }
        }
        // Return result
        i..j
      }
      Err(i) => {
        let start = if i > 0 {
          i - 1
        } else {
          debug_assert_eq!(i, 0);
          i
        };
        start..i
      }
    }
  }

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth_gen(
    sampled_hash: &[T],
    hash_range: Range<u64>,
  ) -> Range<usize> {
    let len = sampled_hash.len();
    let rstart = UInt::from_u64(hash_range.start);
    let rend = UInt::from_u64(hash_range.end);
    let start = match sampled_hash.binary_search(&rstart) {
      Ok(mut i) => {
        if i > 0 {
          i -= 1;
          while i > 0 && sampled_hash[i] == rstart {
            i -= 1;
          }
          i
        } else {
          debug_assert_eq!(i, 0);
          i
        }
      }
      Err(i) => {
        if i > 0 {
          i - 1
        } else {
          i
        }
      }
    };
    let end = match sampled_hash.binary_search(&rend) {
      Ok(mut j) => {
        if j < len {
          j += 1;
          while j < len && sampled_hash[j] == rend {
            j += 1;
          }
        }
        j
      }
      Err(j) => j,
    };
    start..end
  }
}

impl<'a, T: UInt> SHCIndex<T> for BorrowedImplicitSHCIndex<'a, T> {
  fn depth(&self) -> u8 {
    self.depth
  }
  fn n_rows(&self) -> u64 {
    self.n_rows
  }
  fn sampling_step(&self) -> u64 {
    self.n
  }
  fn len(&self) -> usize {
    self.sampled_hash.len()
  }

  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth(&self, hash: u64) -> Range<usize> {
    Self::get_indices_of_elems_containing_hash_at_index_depth_gen(self.sampled_hash, hash)
  }

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<usize> {
    Self::get_indices_of_elems_containing_hash_range_at_index_depth_gen(
      self.sampled_hash,
      hash_range,
    )
  }

  fn write_all_keys<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    /*for k in self.sampled_hash.iter() {
      writer.write_all(k.as_le_bytes())?
    }*/
    let pointer = self.sampled_hash;
    let n_bytes = pointer.len() * size_of::<T>();
    let bytes = unsafe { std::slice::from_raw_parts(pointer.as_ptr() as *const u8, n_bytes) };
    writer.write_all(bytes).map(|()| n_bytes)
  }
}
impl<'a, T: UInt> ImplicitSHCIndex<T> for BorrowedImplicitSHCIndex<'a, T> {}

pub struct OwnedExplicitSHCIndex<T: UInt> {
  index: OwnedImplicitSHCIndex<T>,
  /// Byte number associated to each index value.
  byte_offsets: Box<[u64]>,
}
impl<T: UInt> OwnedExplicitSHCIndex<T> {
  // from_iter_of_positions ...
  /*pub fn new() -> Self {
    TOTO
  }*/
}
impl<T: UInt> SHCIndex<T> for OwnedExplicitSHCIndex<T> {
  fn depth(&self) -> u8 {
    self.index.depth()
  }
  fn n_rows(&self) -> u64 {
    self.index.n_rows()
  }
  fn sampling_step(&self) -> u64 {
    self.index.sampling_step()
  }
  fn len(&self) -> usize {
    self.index.len()
  }

  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth(&self, hash: u64) -> Range<usize> {
    self
      .index
      .get_indices_of_elems_containing_hash_at_index_depth(hash)
  }

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<usize> {
    self
      .index
      .get_indices_of_elems_containing_hash_range_at_index_depth(hash_range)
  }

  fn write_all_keys<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    self.index.write_all_keys(writer)
  }
}

impl<T: UInt> ExplicitSHCIndex<T> for OwnedExplicitSHCIndex<T> {
  fn value(&self, i: usize) -> u64 {
    self.byte_offsets[i]
  }

  fn write_all_values<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    let pointer = self.byte_offsets.as_ref();
    let n_bytes = pointer.len() * size_of::<u64>();
    let bytes = unsafe { std::slice::from_raw_parts(pointer.as_ptr() as *const u8, n_bytes) };
    writer.write_all(bytes).map(|()| n_bytes)
  }
}

pub struct BorrowedExplicitSHCIndex<'a, T: UInt> {
  index: BorrowedImplicitSHCIndex<'a, T>,
  /// Byte number associated to each index value.
  byte_offsets: &'a [u64],
}

impl<'a, T: UInt> SHCIndex<T> for BorrowedExplicitSHCIndex<'a, T> {
  fn depth(&self) -> u8 {
    self.index.depth()
  }
  fn n_rows(&self) -> u64 {
    self.index.n_rows()
  }
  fn sampling_step(&self) -> u64 {
    self.index.sampling_step()
  }
  fn len(&self) -> usize {
    self.index.len()
  }

  /// Returns the indices, in the Sampled Cumulative Index, of the element possibly
  /// containing the given hash, assuming its depth is the same as the index depth.
  fn get_indices_of_elems_containing_hash_at_index_depth(&self, hash: u64) -> Range<usize> {
    self
      .index
      .get_indices_of_elems_containing_hash_at_index_depth(hash)
  }

  /// The input hash value in the range are at the index depth.
  fn get_indices_of_elems_containing_hash_range_at_index_depth(
    &self,
    hash_range: Range<u64>,
  ) -> Range<usize> {
    self
      .index
      .get_indices_of_elems_containing_hash_range_at_index_depth(hash_range)
  }

  fn write_all_keys<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    self.index.write_all_keys(writer)
  }
}

impl<'a, T: UInt> ExplicitSHCIndex<T> for BorrowedExplicitSHCIndex<'a, T> {
  fn value(&self, i: usize) -> u64 {
    self.byte_offsets[i]
  }

  fn write_all_values<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    let pointer = self.byte_offsets;
    let n_bytes = pointer.len() * size_of::<u64>();
    let bytes = unsafe { std::slice::from_raw_parts(pointer.as_ptr() as *const u8, n_bytes) };
    writer.write_all(bytes).map(|()| n_bytes)
  }
}

/// Enum we get when reading a FITS file.
#[derive(Debug)]
pub enum FITSSCIndex {
  ImplicitU32(FitsMMappedSCIndex<u32, Implicit>),
  ImplicitU64(FitsMMappedSCIndex<u64, Implicit>),
  ExplicitU32(FitsMMappedSCIndex<u32, Explicit>),
  ExplicitU64(FitsMMappedSCIndex<u64, Explicit>),
}

impl FITSSCIndex {
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
    let depth = depth_from_n_hash_unsafe(n_rows - 1);
    // Parse other keywords
    // - prepare variables
    let mut end_found = false;
    let mut prodtype_found = false;
    let mut val_schm_found = false;
    let mut is_val_schm_explicit = false;
    let mut dataschm_found = false;
    let mut ordering_found = false;
    let mut dtendian_found = false;
    let mut hpxoder: Option<u8> = None;
    let mut datatype: Option<String> = None;
    let mut n_tot_rows: Option<u64> = None;
    let mut sampling_step: Option<u64> = None;
    let mut indexed_file_name: Option<String> = None;
    let mut indexed_file_len: Option<u64> = None;
    let mut indexed_file_md5: Option<String> = None;
    let mut indexed_file_last_modif_date: Option<SystemTime> = None;
    let mut indexed_colname_lon: Option<String> = None;
    let mut indexed_colname_lat: Option<String> = None;
    let mut date: Option<SystemTime> = None;
    // - do parse
    for kw_record in &mut raw_cards_it {
      match &kw_record[0..8] {
        b"EXTEND  " => check_keyword_and_val(kw_record, b"EXTEND  ", b"F"),
        b"PRODTYPE" => {
          check_keyword_and_str_val(kw_record, b"PRODTYPE", b"HEALPIX SIMPLE CUMUL INDEX")
            .map(|()| prodtype_found = true)
        }
        b"VAL_SCHM" => check_keyword_and_str_val(kw_record, b"VAL_SCHM", b"IMPLICIT")
          .map(|()| val_schm_found = true)
          .or_else(|_err| {
            check_keyword_and_str_val(kw_record, b"VAL_SCHM", b"EXPLICIT").map(|()| {
              val_schm_found = true;
              is_val_schm_explicit = true
            })
          }),
        b"DATASCHM" => check_keyword_and_str_val(kw_record, b"DATASCHM", b"COL_ORIENTED")
          .map(|()| dataschm_found = true),
        b"ORDERING" => check_keyword_and_str_val(kw_record, b"ORDERING", b"NESTED")
          .map(|()| ordering_found = true),
        b"HPXORDER" => parse_uint_val::<u8>(kw_record).map(|v| hpxoder = Some(v)),
        b"DATATYPE" => get_str_val_no_quote(kw_record)
          .map(|v| datatype = Some(String::from_utf8_lossy(v).to_string())),
        b"NTOTROWS" => parse_uint_val::<u64>(kw_record).map(|v| n_tot_rows = Some(v)),
        b"SAMPSTEP" => parse_uint_val::<u64>(kw_record).map(|v| sampling_step = Some(v)),
        b"DTENDIAN" => check_keyword_and_str_val(kw_record, b"DTENDIAN", b"LITTLE")
          .map(|()| dtendian_found = true),
        b"IDXF_NAM" => get_str_val_no_quote(kw_record)
          .map(|v| indexed_file_name = Some(String::from_utf8_lossy(v).to_string())),
        b"IDXF_LEN" => parse_uint_val::<u64>(kw_record).map(|v| indexed_file_len = Some(v)),
        b"IDXF_MD5" => get_str_val_no_quote(kw_record)
          .map(|v| indexed_file_md5 = Some(String::from_utf8_lossy(v).to_string())),
        b"IDXF_LMD" => get_str_val_no_quote(kw_record).map(|v| {
          indexed_file_last_modif_date = unsafe { str::from_utf8_unchecked(v) }
            .parse::<DateTime<Utc>>()
            .ok()
            .map(|dt| dt.into())
        }),
        b"IDXC_LON" => get_str_val_no_quote(kw_record)
          .map(|v| indexed_colname_lon = Some(String::from_utf8_lossy(v).to_string())),
        b"IDXC_LAT" => get_str_val_no_quote(kw_record)
          .map(|v| indexed_colname_lat = Some(String::from_utf8_lossy(v).to_string())),
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
    if !(prodtype_found & val_schm_found & ordering_found & dtendian_found) {
      return Err(FitsError::new_custom(String::from(
        "One of the HEALPIX SIMPLE CUMULATIVE INDEX mandatory cards is missing on the FITS header!",
      )));
    }
    match hpxoder {
      Some(order) => {
        if order == depth {
          Ok(())
        } else {
          Err(FitsError::new_custom(String::from(
            "Number of rows does not match the value of HPXORDER",
          )))
        }
      }
      None => Err(FitsError::new_custom(String::from(
        "The HPXORDER card is missing.",
      ))),
    }?;
    // Check info
    let n_tot_rows = n_tot_rows
      .ok_or_else(|| FitsError::new_custom(String::from("FITS card NTOTROWS is missing")))?;
    let sampling_step = sampling_step
      .ok_or_else(|| FitsError::new_custom(String::from("FITS card SAMPSTEP is missing")))?;
    // Map de data
    let n_bytes_data = n_bytes_per_val * n_rows;
    let mmap = unsafe {
      MmapOptions::new()
        .offset(2880)
        .len(n_bytes_data as usize)
        .map(&file)
        .map_err(FitsError::Io)?
    };
    // Build the right struct
    if is_val_schm_explicit {
      match datatype.as_deref() {
        Some(u32::FITS_DATATYPE_KEY) => Ok(FITSSCIndex::ExplicitU32(FitsMMappedSCIndex::new(
          n_tot_rows,
          sampling_step,
          date,
          indexed_file_name,
          indexed_file_len,
          indexed_file_md5,
          indexed_file_last_modif_date,
          indexed_colname_lon,
          indexed_colname_lat,
          depth,
          mmap,
        ))),
        Some(u64::FITS_DATATYPE_KEY) => Ok(FITSSCIndex::ExplicitU64(FitsMMappedSCIndex::new(
          n_tot_rows,
          sampling_step,
          date,
          indexed_file_name,
          indexed_file_len,
          indexed_file_md5,
          indexed_file_last_modif_date,
          indexed_colname_lon,
          indexed_colname_lat,
          depth,
          mmap,
        ))),
        Some(s) => Err(FitsError::UnexpectedValue {
          keyword: "DATATYPE".to_string(),
          expected: format!(
            "One of: {}, {}.",
            u32::FITS_DATATYPE_KEY,
            u64::FITS_DATATYPE_KEY
          ),
          actual: s.to_string(),
        }),
        None => Err(FitsError::new_custom(String::from(
          "FITS card DATATYPE is missing",
        ))),
      }
    } else {
      match datatype.as_deref() {
        Some(u32::FITS_DATATYPE_KEY) => Ok(FITSSCIndex::ImplicitU32(FitsMMappedSCIndex::new(
          n_tot_rows,
          sampling_step,
          date,
          indexed_file_name,
          indexed_file_len,
          indexed_file_md5,
          indexed_file_last_modif_date,
          indexed_colname_lon,
          indexed_colname_lat,
          depth,
          mmap,
        ))),
        Some(u64::FITS_DATATYPE_KEY) => Ok(FITSSCIndex::ImplicitU64(FitsMMappedSCIndex::new(
          n_tot_rows,
          sampling_step,
          date,
          indexed_file_name,
          indexed_file_len,
          indexed_file_md5,
          indexed_file_last_modif_date,
          indexed_colname_lon,
          indexed_colname_lat,
          depth,
          mmap,
        ))),
        Some(s) => Err(FitsError::UnexpectedValue {
          keyword: "DATATYPE".to_string(),
          expected: format!(
            "One of: {}, {}.",
            u32::FITS_DATATYPE_KEY,
            u64::FITS_DATATYPE_KEY
          ),
          actual: s.to_string(),
        }),
        None => Err(FitsError::new_custom(String::from(
          "FITS card DATATYPE is missing",
        ))),
      }
    }
  }
}

mod seal {
  pub trait Sealed {}
}

/// Used as a marker trait to avoid having to make two distinct `FitsMMappedSCIndex`
/// structures with the same code except when returning the borrowed index.
pub trait IndexType: seal::Sealed {}
#[derive(Debug)]
pub struct Implicit;
impl seal::Sealed for Implicit {}
impl IndexType for Implicit {}
#[derive(Debug)]
pub struct Explicit;
impl seal::Sealed for Explicit {}
impl IndexType for Explicit {}

pub trait BorrowedSHCIndexProvider<'a, T: UInt> {
  type BorrowedSHCIndexProviderType: 'a + SHCIndex<T>;

  fn get_hscindex(&'a self) -> Self::BorrowedSHCIndexProviderType;
}

#[derive(Debug)]
pub struct FitsMMappedSCIndex<T: UInt, I: IndexType> {
  n_tot_rows: u64,
  sampling_size: u64,
  fits_creation_date: Option<SystemTime>,
  indexed_file_name: Option<String>,
  indexed_file_len: Option<u64>,
  indexed_file_md5: Option<String>,
  indexed_file_last_modif_date: Option<SystemTime>,
  indexed_colname_lon: Option<String>,
  indexed_colname_lat: Option<String>,
  depth: u8,
  mmap: Mmap,
  _phantom_t: PhantomData<T>,
  _phantom_i: PhantomData<I>,
}
impl<T: UInt, I: IndexType> FitsMMappedSCIndex<T, I> {
  /// Private, only meant to be called from FITS reader.
  #[allow(clippy::too_many_arguments)]
  fn new(
    n_tot_rows: u64,
    sampling_size: u64,
    fits_creation_date: Option<SystemTime>,
    indexed_file_name: Option<String>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<String>,
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<String>,
    indexed_colname_lat: Option<String>,
    depth: u8,
    mmap: Mmap,
  ) -> Self {
    assert_eq!(
      (n_hash(depth) + 1) * size_of::<T>() as u64,
      mmap.len() as u64
    );
    Self {
      n_tot_rows,
      sampling_size,
      fits_creation_date,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
      depth,
      mmap,
      _phantom_t: PhantomData,
      _phantom_i: PhantomData,
    }
  }

  pub fn get_fits_creation_date(&self) -> Option<&SystemTime> {
    self.fits_creation_date.as_ref()
  }
  pub fn get_indexed_file_name(&self) -> Option<&String> {
    self.indexed_file_name.as_ref()
  }
  pub fn get_indexed_file_len(&self) -> Option<u64> {
    self.indexed_file_len
  }
  pub fn get_indexed_file_md5(&self) -> Option<&String> {
    self.indexed_file_md5.as_ref()
  }
  pub fn get_indexed_file_last_modif_date(&self) -> Option<&SystemTime> {
    self.indexed_file_last_modif_date.as_ref()
  }
  pub fn get_indexed_colname_lon(&self) -> Option<&String> {
    self.indexed_colname_lon.as_ref()
  }
  pub fn get_indexed_colname_lat(&self) -> Option<&String> {
    self.indexed_colname_lat.as_ref()
  }
}

impl<'a, T: UInt> BorrowedSHCIndexProvider<'a, T> for FitsMMappedSCIndex<T, Implicit> {
  type BorrowedSHCIndexProviderType = BorrowedImplicitSHCIndex<'a, T>;

  fn get_hscindex(&self) -> Self::BorrowedSHCIndexProviderType {
    let offset = self.mmap.as_ptr().align_offset(align_of::<T>());
    if offset != 0 {
      // I assume we never enter here, but the assumption had to be tested!
      panic!("Unable to work from MMap, it is not well aligned!!");
    }
    let len = self.mmap.len() / size_of::<T>();
    BorrowedImplicitSHCIndex {
      n: self.sampling_size,
      n_rows: self.n_tot_rows,
      depth: self.depth,
      sampled_hash: unsafe { &*slice_from_raw_parts(self.mmap.as_ptr() as *const T, len) },
    }
  }
}

impl<'a, T: UInt> BorrowedSHCIndexProvider<'a, T> for FitsMMappedSCIndex<T, Explicit> {
  type BorrowedSHCIndexProviderType = BorrowedExplicitSHCIndex<'a, T>;

  fn get_hscindex(&self) -> Self::BorrowedSHCIndexProviderType {
    let len = self.mmap.len() / (size_of::<T>() + size_of::<u64>());
    let (key, val) = self.mmap.split_at(len * size_of::<T>());
    let offset = key.as_ptr().align_offset(align_of::<T>());
    if offset != 0 {
      // I assume we never enter here, but the assumption had to be tested!
      panic!("Unable to work from MMap, it is not well aligned!!");
    }
    let index = BorrowedImplicitSHCIndex {
      n: self.sampling_size,
      n_rows: self.n_tot_rows,
      depth: self.depth,
      sampled_hash: unsafe { &*slice_from_raw_parts(key.as_ptr() as *const T, len) },
    };
    let byte_offsets = unsafe { &*slice_from_raw_parts(val.as_ptr() as *const u64, len) };
    BorrowedExplicitSHCIndex {
      index,
      byte_offsets,
    }
  }
}

/*
/// The result of reading a FITS file, containing a memory map on the data, and from which we can obtain
/// an actual Healpix Cumulative Index object.
#[derive(Debug)]
pub struct FitsMMappedSCIndexImplicit<T: UInt> {
  n_tot_rows: u64,
  sampling_size: u64,
  fits_creation_date: Option<SystemTime>,
  indexed_file_name: Option<String>,
  indexed_file_len: Option<u64>,
  indexed_file_md5: Option<String>,
  indexed_file_last_modif_date: Option<SystemTime>,
  indexed_colname_lon: Option<String>,
  indexed_colname_lat: Option<String>,
  depth: u8,
  mmap: Mmap,
  _phantom: PhantomData<T>,
}
impl<T: UInt> FitsMMappedSCIndexImplicit<T> {
  /// Private, only meant to be called from FITS reader.
  #[allow(clippy::too_many_arguments)]
  fn new(
    n_tot_rows: u64,
    sampling_size: u64,
    fits_creation_date: Option<SystemTime>,
    indexed_file_name: Option<String>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<String>,
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<String>,
    indexed_colname_lat: Option<String>,
    depth: u8,
    mmap: Mmap,
  ) -> Self {
    assert_eq!(
      (n_hash(depth) + 1) * size_of::<T>() as u64,
      mmap.len() as u64
    );
    Self {
      n_tot_rows,
      sampling_size,
      fits_creation_date,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
      depth,
      mmap,
      _phantom: PhantomData,
    }
  }

  pub fn get_fits_creation_date(&self) -> Option<&SystemTime> {
    self.fits_creation_date.as_ref()
  }
  pub fn get_indexed_file_name(&self) -> Option<&String> {
    self.indexed_file_name.as_ref()
  }
  pub fn get_indexed_file_len(&self) -> Option<u64> {
    self.indexed_file_len
  }
  pub fn get_indexed_file_md5(&self) -> Option<&String> {
    self.indexed_file_md5.as_ref()
  }
  pub fn get_indexed_file_last_modif_date(&self) -> Option<&SystemTime> {
    self.indexed_file_last_modif_date.as_ref()
  }
  pub fn get_indexed_colname_lon(&self) -> Option<&String> {
    self.indexed_colname_lon.as_ref()
  }
  pub fn get_indexed_colname_lat(&self) -> Option<&String> {
    self.indexed_colname_lat.as_ref()
  }

  /// Get the actual Healpix Cumulative Index on the MMapped data.
  pub fn get_hscindex(&self) -> BorrowedImplicitSHCIndex<'_, T> {
    let offset = self.mmap.as_ptr().align_offset(align_of::<T>());
    if offset != 0 {
      // I assume we never enter here, but the assumption had to be tested!
      panic!("Unable to work from MMap, it is not well aligned!!");
    }
    let len = self.mmap.len() / size_of::<T>();
    BorrowedImplicitSHCIndex {
      n: self.sampling_size,
      n_rows: self.n_tot_rows,
      depth: self.depth,
      sampled_hash: unsafe { &*slice_from_raw_parts(self.mmap.as_ptr() as *const T, len) },
    }
  }
}

#[derive(Debug)]
pub struct FitsMMappedSCIndexExplicit<T: UInt> {
  n_tot_rows: u64,
  sampling_size: u64,
  fits_creation_date: Option<SystemTime>,
  indexed_file_name: Option<String>,
  indexed_file_len: Option<u64>,
  indexed_file_md5: Option<String>,
  indexed_file_last_modif_date: Option<SystemTime>,
  indexed_colname_lon: Option<String>,
  indexed_colname_lat: Option<String>,
  depth: u8,
  mmap: Mmap,
  _phantom: PhantomData<T>,
}
impl<T: UInt> FitsMMappedSCIndexExplicit<T> {
  /// Private, only meant to be called from FITS reader.
  #[allow(clippy::too_many_arguments)]
  fn new(
    n_tot_rows: u64,
    sampling_size: u64,
    fits_creation_date: Option<SystemTime>,
    indexed_file_name: Option<String>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<String>,
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<String>,
    indexed_colname_lat: Option<String>,
    depth: u8,
    mmap: Mmap,
  ) -> Self {
    assert_eq!(
      (n_hash(depth) + 1) * size_of::<T>() as u64,
      mmap.len() as u64
    );
    Self {
      n_tot_rows,
      sampling_size,
      fits_creation_date,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
      depth,
      mmap,
      _phantom: PhantomData,
    }
  }

  pub fn get_fits_creation_date(&self) -> Option<&SystemTime> {
    self.fits_creation_date.as_ref()
  }
  pub fn get_indexed_file_name(&self) -> Option<&String> {
    self.indexed_file_name.as_ref()
  }
  pub fn get_indexed_file_len(&self) -> Option<u64> {
    self.indexed_file_len
  }
  pub fn get_indexed_file_md5(&self) -> Option<&String> {
    self.indexed_file_md5.as_ref()
  }
  pub fn get_indexed_file_last_modif_date(&self) -> Option<&SystemTime> {
    self.indexed_file_last_modif_date.as_ref()
  }
  pub fn get_indexed_colname_lon(&self) -> Option<&String> {
    self.indexed_colname_lon.as_ref()
  }
  pub fn get_indexed_colname_lat(&self) -> Option<&String> {
    self.indexed_colname_lat.as_ref()
  }

  /// Get the actual Healpix Cumulative Index on the MMapped data.
  pub fn get_hscindex(&self) -> BorrowedExplicitSHCIndex<'_, T> {
    let len = self.mmap.len() / (size_of::<T>() + size_of::<u64>());
    let (key, val) = self.mmap.split_at(len * size_of::<T>());
    let offset = key.as_ptr().align_offset(align_of::<T>());
    if offset != 0 {
      // I assume we never enter here, but the assumption had to be tested!
      panic!("Unable to work from MMap, it is not well aligned!!");
    }
    let index = BorrowedImplicitSHCIndex {
      n: self.sampling_size,
      n_rows: self.n_tot_rows,
      depth: self.depth,
      sampled_hash: unsafe { &*slice_from_raw_parts(key.as_ptr() as *const T, len) },
    };
    let byte_offsets = unsafe { &*slice_from_raw_parts(val.as_ptr() as *const u64, len) };
    BorrowedExplicitSHCIndex {
      index,
      byte_offsets,
    }
  }
}
*/
