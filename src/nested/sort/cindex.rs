//! Defines a HEALPix Cumulative Index for a HEALPix sorted files.

// Two types of index:
// * naive (starting rcno for given healpix value)
// * constant (index of k*nème source) => 2x binary search for an HEALPix range
//   + index in the array * k => number of rows (implicit)
//   + store (healpix, startig byte)
//   + zuniq indexing (zuniq, starting byte)
//   =>

// AUTRE TYPE D'INDEX:
// stocker à l'index 'i' le cellule HEALPix de la source à la ligne n * i.
// On trouve le range de ligne pour un Range Healpix avec seulement deux binary search!
// * e.g. pour 2.10^9 sources, si on fait prend n = 2_000 lignes, ca fait 1_000_000 ligne d'index,
//   soit moins de 8 MB (idx29 sur u64).

use std::{
  array::TryFromSliceError,
  cmp::Ordering,
  collections::BTreeMap,
  convert::{TryFrom, TryInto},
  fs::File,
  io::{BufWriter, Error as IoError, Write},
  marker::PhantomData,
  mem::align_of,
  ops::{AddAssign, Range, Sub},
  path::Path,
  ptr::slice_from_raw_parts,
  slice::{from_raw_parts, ChunksMut},
  str,
  time::SystemTime,
};

use chrono::{DateTime, SecondsFormat, Utc};
use log::debug;
use memmap2::{Mmap, MmapOptions};
use num_traits::{FromBytes, ToBytes, Zero};

use crate::{
  depth_from_n_hash_unsafe, n_hash,
  nested::map::{
    fits::{
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
    skymap::SkyMap,
    HHash,
  },
};

/// Trait defining the type of (cumulative) value a Healpix Cumulative Index may contains.
pub trait HCIndexValue:
  Sized
  + Zero
  + Sub<Output = Self>
  + AddAssign
  + ToBytes
  + FromBytes<Bytes: for<'a> TryFrom<&'a [u8], Error = TryFromSliceError>>
  + Clone
  + Copy
  + PartialEq
  + 'static
{
  const FITS_DATATYPE: &'static str;
}

impl HCIndexValue for u32 {
  const FITS_DATATYPE: &'static str = "u32";
  /*fn binary_search_u32(data: &[u8], hash: u32) -> Self {
    const SIZE_OF_H: usize = size_of::<u32>();

    let ptr = data.as_ptr() as *const ([u8; size_of::<u32>()], [u8; size_of::<Self>()]);
    let raw_entries = unsafe {
      std::slice::from_raw_parts(ptr, data.len() / (size_of::<u32>() + size_of::<Self>()))
    };
    /*let raw_entries =
    unsafe { transmute::<&[u8], &[([u8; size_of::<u32>()], [u8; size_of::<Self>()])]>(data) };*/
    match raw_entries.binary_search_by_key(&hash, |slice| u32::from_le_bytes(slice.0)) {
      Ok(i) => Self::from_le_bytes(raw_entries[i].1),
      Err(i) => {
        if i == 0 {
          Self::zero()
        } else {
          Self::from_le_bytes(raw_entries[i - 1].1)
        }
      }
    }
  }*/
}

impl HCIndexValue for u64 {
  const FITS_DATATYPE: &'static str = "u64";
}

impl HCIndexValue for f32 {
  const FITS_DATATYPE: &'static str = "f32";
}

impl HCIndexValue for f64 {
  const FITS_DATATYPE: &'static str = "f64";
}

pub enum HCIndexShape {
  Implicit,
  Explicit,
}

/// Defines a Healpix Cumulative Index, i.e. an index in which:
/// * keys are HEALPix cells
/// * all the cells of a given index have the same depth
/// * values are cumulative values, i.e. the value necessarily increase as the HEALPix cell number (hash value) increases.
///
/// # Use cases
/// The current main use case is indexing HEALPix ordered relational tables stored in CSV, VOTable, FITS, ... files.
/// The cumulative index allow to quickly locate physically in a file (from a row numbers or row starting bytes)
/// all rows spatially located in any given HEALPix cell of depth lower or equal to the index depth.
/// Such simple indexes have been use in RCF files, the CDS XMatch Service or the ExXmatch prototype
/// since 2010.
/// We now want to better formalize it in a single generic Rust code in the core HEALPix library,
/// and to provide with a FITS serialization allowing to attach some metadata to those index while providing
/// direct mapping to the structure on the disk without byte re-ordering manipulations.
///
/// It although offer a very quick way to provide the number of rows (or or bytes) covered by a given
/// MOC, without having to sum over all values in a given range.
///
/// It could also be used with cumulative probabilities to quickly compute the taol probability covered
///by an input MOC, on possibly large proba skymap, without having to load the full file in memory.
pub trait HCIndex {
  /// Type of value the index contains.
  type V: HCIndexValue;
  const SIZE_OF_V: usize = size_of::<Self::V>();

  fn shape(&self) -> HCIndexShape;

  /// Returns the HEALPix depth of the index.
  fn depth(&self) -> u8;

  /// Returns the cumulative value associated to the given index.
  /// The starting cumulative value associated to a HEALPix cell of index `i` is obtained by `get(i)`
  /// while the ending cumulative value is obtained by `get(i+1)`.
  /// It means that the largest index equals the number of HEALPix cells at the index depth, and the
  /// number of stored values equals this number plus one.
  ///
  /// # Panics
  /// If `index` larger than `n_hash(depth)`.
  fn get(&self, hash: u64) -> Self::V;

  /// Returns both the starting and the ending values associated with the given input `hash` at the index `depth`.
  fn get_with_hash_at_index_depth(&self, hash: u64) -> Range<Self::V> {
    self.get_with_range_at_index_depth(hash..hash + 1)
  }
  /// Returns both the starting and the ending values associated with the given input `range` at the index `depth`.
  fn get_with_range_at_index_depth(&self, range: Range<u64>) -> Range<Self::V> {
    self.get(range.start)..self.get(range.end)
  }

  /// Returns the noncumulative value associated with the given input `hash` at the index `depth`.
  fn get_noncumulative_with_hash_at_index_depth(&self, hash: u64) -> Self::V {
    self.get_noncumulative_with_range_at_index_depth(hash..hash + 1)
  }
  /// Returns the noncumulative value associated with the given input `range` at the index `depth`.
  fn get_noncumulative_with_range_at_index_depth(&self, range: Range<u64>) -> Self::V {
    self.get(range.end) - self.get(range.start)
  }

  fn get_cell(&self, depth: u8, hash: u64) -> Range<Self::V> {
    match depth.cmp(&self.depth()) {
      Ordering::Equal => self.get_with_hash_at_index_depth(hash),
      Ordering::Greater => {
        let twice_dd = (depth - self.depth()) << 1;
        self.get_with_hash_at_index_depth(hash >> twice_dd)
      }
      Ordering::Less => {
        let twice_dd = (self.depth() - depth) << 1;
        self.get_with_range_at_index_depth((hash << twice_dd)..((hash + 1) << twice_dd))
      }
    }
  }

  fn get_cell_noncumulative(&self, depth: u8, hash: u64) -> Self::V {
    match depth.cmp(&self.depth()) {
      Ordering::Equal => self.get_noncumulative_with_hash_at_index_depth(hash),
      Ordering::Greater => {
        let twice_dd = (depth - self.depth()) << 1;
        self.get_noncumulative_with_hash_at_index_depth(hash >> twice_dd)
      }
      Ordering::Less => {
        let twice_dd = (self.depth() - depth) << 1;
        let start = hash << twice_dd;
        let end = (hash + 1) << twice_dd;
        self.get_noncumulative_with_range_at_index_depth(start..end)
      }
    }
  }

  /// If the given input `depth` is larger than the index `depth`, then the result is the one
  /// obtain from degrading first the given input range to the index depth.
  fn get_range(&self, depth: u8, range: Range<u64>) -> Range<Self::V> {
    match depth.cmp(&self.depth()) {
      Ordering::Equal => self.get_with_range_at_index_depth(range),
      Ordering::Greater => {
        let twice_dd = (depth - self.depth()) << 1;
        let start = range.start >> twice_dd;
        let end = (range.end + !(u64::MAX << twice_dd)) >> twice_dd;
        self.get_with_range_at_index_depth(start..end)
      }
      Ordering::Less => {
        let twice_dd = (self.depth() - depth) << 1;
        let start = range.start << twice_dd;
        let end = range.end << twice_dd;
        self.get_with_range_at_index_depth(start..end)
      }
    }
  }

  /// If the given input `depth` is larger than the index `depth`, then the result is the one
  /// obtain from degrading first the given input range to the index depth.
  fn get_range_noncumulative(&self, depth: u8, range: Range<u64>) -> Self::V {
    match depth.cmp(&self.depth()) {
      Ordering::Equal => self.get_noncumulative_with_range_at_index_depth(range),
      Ordering::Greater => {
        let twice_dd = (depth - self.depth()) << 1;
        let start = range.start >> twice_dd;
        let end = (range.end + !(u64::MAX << twice_dd)) >> twice_dd;
        self.get_noncumulative_with_range_at_index_depth(start..end)
      }
      Ordering::Less => {
        let twice_dd = (self.depth() - depth) << 1;
        let start = range.start << twice_dd;
        let end = range.end << twice_dd;
        self.get_noncumulative_with_range_at_index_depth(start..end)
      }
    }
  }

  /// For the fits serialization. Must write the `n_hash(depth) + 1` values, in ascending order,
  /// and in little-endian, in the given writer.
  /// Exactly `(n_hash(depth) + 1) * size_of::<Self::V>` bytes must be written.
  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError>;

  /// For the fits serialization. Must write all tuples, in ascending hpx order,
  /// and in little-endian, in the given writer.
  /// # WARNING
  /// If the depth is lower or equal to 13, cell hash value **must** be written using `u32`,
  /// else using `u64`.
  /// Exactly TBD!! bytes must be written.
  fn write_all_values_explicit<W: Write>(&self, writer: W) -> Result<usize, IoError>;

  /// Returns the size, in bytes, of the implicit array of values.
  fn byte_size_implicit(&self) -> u64 {
    size_of::<Self::V>() as u64 * n_hash(self.depth())
  }

  /// Returns the size, in bytes, of the list of explicit `(hash, value)` tuples.
  fn byte_size_explicit(&self) -> u64;

  /// Returns the "best" representation ot be used to serialize (in FITS) this index.
  /// The "best" is determined from the size, in byte, of both representation.
  /// # Params
  /// * `impl_over_expl_limit_ratio` limit on the ratio of the implicit byte size
  /// over the explicit byte size.
  ///     + above the limit, the explicit representation is chosen.
  ///     + below the limit, the implicit representation is chosen.
  /// (This allows to put a size limit on the explicit representation: looking to a value in 'explicit'
  /// requires a binary search, which is longer than the direct access allowed by `implicit`).
  fn best_representation(&self, impl_over_expl_limit_ratio: f64) -> HCIndexShape {
    if self.byte_size_implicit() as f64
      <= self.byte_size_explicit() as f64 * impl_over_expl_limit_ratio
    {
      HCIndexShape::Implicit
    } else {
      HCIndexShape::Explicit
    }
  }

  #[allow(clippy::too_many_arguments)]
  fn to_fits<W: Write>(
    &self,
    writer: W,
    shape: HCIndexShape,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    match shape {
      HCIndexShape::Implicit => self.to_fits_implicit(
        writer,
        indexed_file_name,
        indexed_file_len,
        indexed_file_md5,
        indexed_file_last_modif_date,
        indexed_colname_lon,
        indexed_colname_lat,
      ),
      HCIndexShape::Explicit => self.to_fits_explicit(
        writer,
        indexed_file_name,
        indexed_file_len,
        indexed_file_md5,
        indexed_file_last_modif_date,
        indexed_colname_lon,
        indexed_colname_lat,
      ),
    }
  }

  #[allow(clippy::too_many_arguments)]
  fn to_fits_file<P: AsRef<Path>>(
    &self,
    path: P,
    shape: HCIndexShape,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    File::create(path).map_err(FitsError::Io).and_then(|file| {
      self.to_fits(
        BufWriter::new(file),
        shape,
        indexed_file_name,
        indexed_file_len,
        indexed_file_md5,
        indexed_file_last_modif_date,
        indexed_colname_lon,
        indexed_colname_lat,
      )
    })
  }

  /// # Design choices
  /// For compactness, we decided to write the cumulative index into the primary HDU (an empty
  /// primary HDU = 2880 bytes, i.e. almost 3 kB of useless data).
  /// In addition, we want to be able to directly map the array of values in memory without having to
  /// re-order bytes on modern systems (FITS values are encoded in big-endian).
  /// We thus de-compose each value in an array of little-endian ordered bytes.
  /// Regular softwares will thus not load "wrong" integer values but array of bytes, and
  /// cumulative-index-aware softwares (recognizing so from FITS cards) will know how to interpret
  /// those arrays.
  /// * if `NAXIS1=4`, the stored type may be either a `u32` or a `f32`, in little-endian.
  /// * if `NAXIS1=8`, the stored type may be either a `u64` or a `f64`, in little-endian.
  ///   Here the list of card htat must appear, in this order:
  /// * `PRODTYPE`, mandatory: **must** equal `HEALPIX CUMUL INDEX`
  /// * `ORDERING`, mandatory: **must** equal `NESTED`
  /// * `INDXSCHM`, mandatory: so far only `IMPLICIT` is supported. An `EXPLICIT` would rely on
  ///   an ordered ZINDEX followed by the CVALUE.
  /// * `INDXCOV`, mandatory if `INDXSCHM = IMPLICIT`: so far only `FULLSKY` is supported, but we
  ///   may support `PARTIAL` with 2 additional FITS cards providing the starting and ending HEALPix indices.
  ///   It may be usefull, e.g. to index a file containing data located in a given (large) HEALPix cell,
  ///   e.g. individual HATS parquet files.
  /// * `HPXORDER`, mandatory: provide the depth (order) of the cumulative MAP.
  /// * `DATATYPE`, mandatory: **must** be one of `u32`, `u64`, `f32`, `f64`.
  /// * `DTENDIAN`, mandatory: **must** equal `LITTLE`, made to provide information so human readers can understand/guess the data structure.
  /// * `IDXF_NAM`, optional: name of the indexed file if the HCI is a file row index.
  /// * `IDXF_LEN`, optional: length, in bytes, of the indexed file if the HCI is a file row index.
  /// * `IDXF_MD5`, optional: md5 of the indexed file if the HCI is a file row index.
  /// * `IDXF_LMD`, optional: last modification date of the indexed file if the HCI is a file row index.
  /// * `IDXC_LON`, optional: name, in the indexed file, of the column containing the longitudes used to compute HEALPix index.
  /// * `IDXC_LAT`, optional: name, in the indexed file, of the column containing the latitudes used to compute HEALPix index.
  #[allow(clippy::too_many_arguments)]
  fn to_fits_implicit<W: Write>(
    &self,
    mut writer: W,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    let n_values = n_hash(self.depth()) + 1;
    // Perpare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    write_uint_mandatory_keyword_record(
      it.next().unwrap(),
      b"NAXIS1  ",
      size_of::<Self::V>() as u64,
    ); // Len of data axis 1 = number of bytes in dataype
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", n_values); // Len of data axis 2 = n_hash(depth)+1
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..31].copy_from_slice(b"PRODTYPE= 'HEALPIX CUMUL INDEX'"); // Product type
    it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
    it.next().unwrap()[0..20].copy_from_slice(b"INDXSCHM= 'IMPLICIT'");
    it.next().unwrap()[0..20].copy_from_slice(b"INDXCOV = 'FULLSKY '");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth() as u64); // Must be NSIDE in case of RING
    write_str_keyword_record(it.next().unwrap(), b"DATATYPE", Self::V::FITS_DATATYPE);
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    append_optional_card(
      &mut it,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..]).map_err(FitsError::Io)?;
    // Write the data part
    self
      .write_all_values_implicit(&mut writer)
      .map_err(FitsError::Io)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }

  /// # Design choices
  /// For compactness, we decided to write the cumulative index into the primary HDU (an empty
  /// primary HDU = 2880 bytes, i.e. almost 3 kB of useless data).
  /// Both elements of the `(INDEX, VALUES)` tuples are stored in little-endian ordered bytes.
  /// Regular softwares will thus not load "wrong" integer values but array of bytes, and
  /// cumulative-index-aware softwares (recognizing so from FITS cards) will know how to interpret
  /// those arrays.
  /// * if `NAXIS1=4`, the stored type may be either a `u32` or a `f32`, in little-endian.
  /// * if `NAXIS1=8`, the stored type may be either a `u64` or a `f64`, in little-endian.
  ///   Here the list of card htat must appear, in this order:
  /// * `PRODTYPE`, mandatory: **must** equal `HEALPIX CUMUL INDEX`
  /// * `ORDERING`, mandatory: **must** equal `NESTED`
  /// * `INDXSCHM`, mandatory: `EXPLICIT`
  /// * `HPXORDER`, mandatory: provide the depth (order) of the cumulative MAP.
  /// * `INDXTYPE`, mandatory: **must** be `u32` if `HPXORDER <= 13`, `u64` else. (It is the type of stored HEALPix index).
  /// * `DATATYPE`, mandatory: **must** be one of `u32`, `u64`, `f32`, `f64`.
  /// * `DTENDIAN`, mandatory: **must** equal `LITTLE`, made to provide information so human readers can understand/guess the data structure.
  /// * `IDXF_NAM`, optional: name of the indexed file if the HCI is a file row index.
  /// * `IDXF_LEN`, optional: length, in bytes, of the indexed file if the HCI is a file row index.
  /// * `IDXF_MD5`, optional: md5 of the indexed file if the HCI is a file row index.
  /// * `IDXF_LMD`, optional: last modification date of the indexed file if the HCI is a file row index.
  /// * `IDXC_LON`, optional: name, in the indexed file, of the column containing the longitudes used to compute HEALPix index.
  /// * `IDXC_LAT`, optional: name, in the indexed file, of the column containing the latitudes used to compute HEALPix index.
  #[allow(clippy::too_many_arguments)]
  fn to_fits_explicit<W: Write>(
    &self,
    mut writer: W,
    indexed_file_name: Option<&str>,
    indexed_file_len: Option<u64>,
    indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
    indexed_file_last_modif_date: Option<SystemTime>,
    indexed_colname_lon: Option<&str>,
    indexed_colname_lat: Option<&str>,
  ) -> Result<(), FitsError> {
    let use_u32 = self.depth() <= 13;
    let (size_of_h, h_type) = if use_u32 {
      (size_of::<u32>(), u32::FITS_DATATYPE)
    } else {
      (size_of::<u64>(), u64::FITS_DATATYPE)
    };
    let n_values = n_hash(self.depth()) + 1;
    // Perpare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    write_uint_mandatory_keyword_record(
      it.next().unwrap(),
      b"NAXIS1  ",
      (size_of_h + size_of::<Self::V>()) as u64,
    ); // Len of data axis 1 = number of bytes in dataype
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", n_values); // Len of data axis 2 = n_hash(depth)+1
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..31].copy_from_slice(b"PRODTYPE= 'HEALPIX CUMUL INDEX'"); // Product type
    it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
    it.next().unwrap()[0..20].copy_from_slice(b"INDXSCHM= 'EXPLICIT'");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth() as u64); // Must be NSIDE in case of RING
    write_str_keyword_record(it.next().unwrap(), b"INDXTYPE", h_type);
    write_str_keyword_record(it.next().unwrap(), b"DATATYPE", Self::V::FITS_DATATYPE);
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    append_optional_card(
      &mut it,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..]).map_err(FitsError::Io)?;
    // Write the data part
    self
      .write_all_values_implicit(&mut writer)
      .map_err(FitsError::Io)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }
}

fn append_optional_card(
  it: &mut ChunksMut<u8>,
  indexed_file_name: Option<&str>,
  indexed_file_len: Option<u64>,
  indexed_file_md5: Option<[u8; 32]>, // 32 hex characters (128 bits)
  indexed_file_last_modif_date: Option<SystemTime>,
  indexed_colname_lon: Option<&str>,
  indexed_colname_lat: Option<&str>,
) {
  let indexed_file_last_modif_date = indexed_file_last_modif_date.map(DateTime::<Utc>::from);
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
}

#[derive(Debug, PartialEq)]
pub struct OwnedCIndex<T: HCIndexValue> {
  depth: u8,
  values: Box<[T]>,
}
impl<T: HCIndexValue> OwnedCIndex<T> {
  pub fn new_unchecked(depth: u8, values: Box<[T]>) -> Self {
    let len = n_hash(depth) + 1; // + 1 to get the total number of rows in the file at index `n_hash`
    assert_eq!(len as usize, values.len());
    Self { depth, values }
  }
}
impl<T: HCIndexValue> HCIndex for OwnedCIndex<T> {
  type V = T;

  fn shape(&self) -> HCIndexShape {
    HCIndexShape::Implicit
  }

  fn depth(&self) -> u8 {
    self.depth
  }

  fn get(&self, hash: u64) -> Self::V {
    self.values[hash as usize]
  }
  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    /* Replaced by unsafe by arguably faster code
    for v in self.values {
      writer.write_all(v.to_le_bytes().as_ref())?;
    }*/
    write_all_values_implicit(self.values.as_ref(), writer)
  }

  fn write_all_values_explicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_explicit_from_implicit(self.depth, self.values.as_ref(), writer)
  }

  fn byte_size_explicit(&self) -> u64 {
    compute_byte_size_explicit(self.values.as_ref())
  }
}

/// WARNING: **you** have to ensure that the cumulated value still holds in the SkyMap type!!
/// (E.g. it is not the case for u32 if the total count is larger that u32::MAX!)
impl<'a, T: 'a + HCIndexValue, S: SkyMap<'a, ValueType = T>> From<&'a S> for OwnedCIndex<T> {
  fn from(skymap: &'a S) -> Self {
    let depth = skymap.depth();
    let n = n_hash(depth) + 1;
    let mut res = Vec::<T>::with_capacity(n as usize);
    let mut acc = T::zero();
    for v in skymap.values() {
      res.push(acc);
      acc.add_assign(*v);
    }
    res.push(acc);
    Self::new_unchecked(depth, res.into_boxed_slice())
  }
}

/// Implementation based on a borrowed array, possibly directly mapped from the FITS serialization.
#[derive(Debug, PartialEq)]
pub struct BorrowedCIndex<'a, T: HCIndexValue> {
  depth: u8,
  values: &'a [T],
}
impl<T: HCIndexValue> HCIndex for BorrowedCIndex<'_, T> {
  type V = T;
  fn shape(&self) -> HCIndexShape {
    HCIndexShape::Implicit
  }
  fn depth(&self) -> u8 {
    self.depth
  }
  fn get(&self, hash: u64) -> T {
    self.values[hash as usize]
  }
  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    /* Replaced by unsafe by arguably faster code
    for v in self.values {
      writer.write_all(v.to_le_bytes().as_ref())?;
    }*/
    write_all_values_implicit(self.values, writer)
  }

  fn write_all_values_explicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_explicit_from_implicit(self.depth, self.values, writer)
  }

  fn byte_size_explicit(&self) -> u64 {
    compute_byte_size_explicit(self.values)
  }
}

impl<'a, T: HCIndexValue> From<&'a OwnedCIndex<T>> for BorrowedCIndex<'a, T> {
  fn from(cindex: &'a OwnedCIndex<T>) -> Self {
    Self {
      depth: cindex.depth,
      values: cindex.values.as_ref(),
    }
  }
}

fn compute_byte_size_explicit<T: HCIndexValue>(values: &[T]) -> u64 {
  // we start at 1 because 'curr' starts at index 1 and the first element must be included
  values.iter().zip(values.iter().skip(1)).fold(
    1,
    |count, (prev, curr)| if curr != prev { count + 1 } else { count },
  )
}

fn write_all_values_implicit<T: HCIndexValue, W: Write>(
  values: &[T],
  mut writer: W,
) -> Result<usize, IoError> {
  let len = size_of_val(values);
  let ptr = values.as_ptr();
  let offset = ptr.align_offset(align_of::<u8>());
  // I suppose that align with u8 is never a problem, but we test the assumption just in case...
  assert_eq!(offset, 0);
  let bytes: &[u8] = unsafe { std::slice::from_raw_parts(ptr as *const u8, len) };
  writer.write_all(bytes).map(|()| len)
}

fn write_all_values_explicit_from_implicit<T: HCIndexValue, W: Write>(
  depth: u8,
  values: &[T],
  mut writer: W,
) -> Result<usize, IoError> {
  let mut len = 0;
  if depth <= 13 {
    let mut it = values.iter().enumerate();
    if let Some((i, v)) = it.next() {
      let mut prev = *v;
      writer
        .write_all((i as u32).to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
      len += 1;
      for (i, v) in it {
        if *v != prev {
          writer
            .write_all((i as u32).to_le_bytes().as_ref())
            .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
          len += 1;
          prev = *v;
        }
      }
    }
    Ok((size_of::<u32>() + size_of::<T>()) * len)
  } else {
    /*for (i, v) in values.iter().enumerate() {
      writer
        .write_all((i as u64).to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
    }*/
    let mut it = values.iter().enumerate();
    if let Some((i, v)) = it.next() {
      let mut prev = *v;
      writer
        .write_all((i as u64).to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
      len += 1;
      for (i, v) in it {
        if *v != prev {
          writer
            .write_all((i as u64).to_le_bytes().as_ref())
            .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
          len += 1;
          prev = *v;
        }
      }
    }
    Ok((size_of::<u64>() + size_of::<T>()) * len)
  }
}

fn write_all_values_explicit<'a, H: HHash, T: HCIndexValue, I, W: Write>(
  depth: u8,
  it: I,
  mut writer: W,
) -> Result<usize, IoError>
where
  I: Iterator<Item = (&'a H, &'a T)>,
{
  let mut len = 0;
  if depth <= 13 {
    for (h, v) in it {
      len += 1;
      writer
        .write_all(HHash::to_u32(h).to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
    }
    Ok((size_of::<u32>() + size_of::<T>()) * len)
  } else {
    for (h, v) in it {
      len += 1;
      writer
        .write_all(HHash::to_u64(h).to_le_bytes().as_ref())
        .and_then(|()| writer.write_all(&(v.to_le_bytes().as_ref())))?;
    }
    Ok((size_of::<u64>() + size_of::<T>()) * len)
  }
}

fn write_all_values_implicit_from_explicit<H: HHash, T: HCIndexValue, I, W: Write>(
  depth: u8,
  it: I,
  mut writer: W,
) -> Result<usize, IoError>
where
  I: Iterator<Item = (H, T)>,
{
  let mut curr_h = 0;
  let mut curr_v = T::zero();
  let end = n_hash(depth);
  for (h, v) in it {
    let to = HHash::to_u64(&h);
    // Fill missing values by current value.
    while curr_h < to {
      writer.write_all(curr_v.to_le_bytes().as_ref())?;
    }
    // Write current value
    writer.write_all(v.to_le_bytes().as_ref())?;
    curr_v = v;
    curr_h = to;
  }
  while curr_h < end {
    writer.write_all(curr_v.to_le_bytes().as_ref())?;
  }
  Ok(end as usize * size_of::<T>())
}

impl<'a, T: HCIndexValue> BorrowedCIndex<'a, T> {
  /// For each HEALPix cell at the given depth, `values` contains a cumulative quantity.
  /// The value at index `hash` contains the "starting" (inclusive) value for cell of hash value `hash`
  /// while its "ending" (exclusive) value is stored at index `hash + 1`.
  /// Hence, the largest possible index is `n_hash`, the number of cells at the given `depth`,
  /// and `value` is of size `n_hash + 1`.
  /// If the value is a number of rows, `values[0] = 0` and `values[n_hash(depth)] = n_rows_tot`.
  /// # Panics
  /// If the input `values` size does not equal the number of HEALPIx cells at the given depth, plus one.
  pub fn new(depth: u8, values: &'a [T]) -> Self {
    let n = n_hash(depth) + 1;
    assert_eq!(n as usize, values.len());
    Self { depth, values }
  }
}

/// Enum we get when reading a FITS file.
#[derive(Debug)]
pub enum FITSCIndex {
  ImplicitU32(FitsMMappedCIndexImplicit<u32>),
  ImplicitU64(FitsMMappedCIndexImplicit<u64>),
  ImplicitF32(FitsMMappedCIndexImplicit<f32>),
  ImplicitF64(FitsMMappedCIndexImplicit<f64>),
  ExplicitU32U32(FitsMMappedCIndexExplicit<u32, u32>),
  ExplicitU32U64(FitsMMappedCIndexExplicit<u32, u64>),
  ExplicitU32F32(FitsMMappedCIndexExplicit<u32, f32>),
  ExplicitU32F64(FitsMMappedCIndexExplicit<u32, f64>),
  ExplicitU64U32(FitsMMappedCIndexExplicit<u64, u32>),
  ExplicitU64U64(FitsMMappedCIndexExplicit<u64, u64>),
  ExplicitU64F32(FitsMMappedCIndexExplicit<u64, f32>),
  ExplicitU64F64(FitsMMappedCIndexExplicit<u64, f64>),
}
impl FITSCIndex {
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
    let mut end_found = false;
    let mut prodtype_found = false;
    let mut ordering_found = false;
    let mut indxschm_found = false;
    let mut is_ndxschm_explicit = false;
    let mut indxcov_found = false;
    let mut dtendian_found = false;
    let mut hpxoder: Option<u8> = None;
    let mut indxtype: Option<String> = None; // EXPLICIT ONLY
    let mut datatype: Option<String> = None;
    let mut indexed_file_name: Option<String> = None;
    let mut indexed_file_len: Option<u64> = None;
    let mut indexed_file_md5: Option<String> = None;
    let mut indexed_file_last_modif_date: Option<SystemTime> = None;
    let mut indexed_colname_lon: Option<String> = None;
    let mut indexed_colname_lat: Option<String> = None;
    let mut date: Option<SystemTime> = None;

    for kw_record in &mut raw_cards_it {
      match &kw_record[0..8] {
        b"EXTEND  " => check_keyword_and_val(kw_record, b"EXTEND  ", b"F"),
        b"PRODTYPE" => check_keyword_and_str_val(kw_record, b"PRODTYPE", b"HEALPIX CUMUL INDEX")
          .map(|()| prodtype_found = true),
        b"ORDERING" => check_keyword_and_str_val(kw_record, b"ORDERING", b"NESTED")
          .map(|()| ordering_found = true),
        b"INDXSCHM" => check_keyword_and_str_val(kw_record, b"INDXSCHM", b"IMPLICIT")
          .map(|()| indxschm_found = true)
          .or_else(|_err| {
            check_keyword_and_str_val(kw_record, b"INDXSCHM", b"EXPLICIT").map(|()| {
              indxschm_found = true;
              is_ndxschm_explicit = true
            })
          }),
        b"INDXCOV " => check_keyword_and_str_val(kw_record, b"INDXCOV ", b"FULLSKY")
          .map(|()| indxcov_found = true),
        b"HPXORDER" => parse_uint_val::<u8>(kw_record).map(|v| hpxoder = Some(v)),
        b"INDXTYPE" => get_str_val_no_quote(kw_record)
          .map(|v| indxtype = Some(String::from_utf8_lossy(v).to_string())),
        b"DATATYPE" => get_str_val_no_quote(kw_record)
          .map(|v| datatype = Some(String::from_utf8_lossy(v).to_string())),
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
    if !(prodtype_found & ordering_found & indxschm_found & indxcov_found & dtendian_found) {
      return Err(FitsError::new_custom(String::from(
        "One of the HEALPIX CUMULATIVE INDEX mandatory cards is missing on the FITS header!",
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
    if is_ndxschm_explicit {
      if depth <= 13 {
        match datatype.as_deref() {
          Some(u32::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU32U32(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(u64::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU32U64(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(f32::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU32F32(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(f64::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU32F64(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(s) => Err(FitsError::UnexpectedValue {
            keyword: "DATATYPE".to_string(),
            expected: format!("One of: {}, {}.", u32::FITS_DATATYPE, u64::FITS_DATATYPE),
            actual: s.to_string(),
          }),
          None => Err(FitsError::new_custom(String::from(
            "FITS card DATATYPE is missing",
          ))),
        }
      } else {
        match datatype.as_deref() {
          Some(u32::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU64U32(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(u64::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU64U64(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(f32::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU64F32(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(f64::FITS_DATATYPE) => {
            Ok(FITSCIndex::ExplicitU64F64(FitsMMappedCIndexExplicit::new(
              date,
              indexed_file_name,
              indexed_file_len,
              indexed_file_md5,
              indexed_file_last_modif_date,
              indexed_colname_lon,
              indexed_colname_lat,
              depth,
              mmap,
            )))
          }
          Some(s) => Err(FitsError::UnexpectedValue {
            keyword: "DATATYPE".to_string(),
            expected: format!("One of: {}, {}.", u32::FITS_DATATYPE, u64::FITS_DATATYPE),
            actual: s.to_string(),
          }),
          None => Err(FitsError::new_custom(String::from(
            "FITS card DATATYPE is missing",
          ))),
        }
      }
    } else {
      match datatype.as_deref() {
        Some(u32::FITS_DATATYPE) => Ok(FITSCIndex::ImplicitU32(FitsMMappedCIndexImplicit::new(
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
        Some(u64::FITS_DATATYPE) => Ok(FITSCIndex::ImplicitU64(FitsMMappedCIndexImplicit::new(
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
        Some(f32::FITS_DATATYPE) => Ok(FITSCIndex::ImplicitF32(FitsMMappedCIndexImplicit::new(
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
        Some(f64::FITS_DATATYPE) => Ok(FITSCIndex::ImplicitF64(FitsMMappedCIndexImplicit::new(
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
          expected: format!("One of: {}, {}.", u32::FITS_DATATYPE, u64::FITS_DATATYPE),
          actual: s.to_string(),
        }),
        None => Err(FitsError::new_custom(String::from(
          "FITS card DATATYPE is missing",
        ))),
      }
    }
  }
}

/// The result of reading a FITS file, containing a memory map on the data, and from which we can obtain
/// an actual Healpix Cumulative Index object.
#[derive(Debug)]
pub struct FitsMMappedCIndexImplicit<T: HCIndexValue> {
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
impl<T: HCIndexValue> FitsMMappedCIndexImplicit<T> {
  /// Private, only meant to be called from FITS reader.
  #[allow(clippy::too_many_arguments)]
  fn new(
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
  pub fn get_hcindex(&self) -> BorrowedCIndex<'_, T> {
    let offset = self.mmap.as_ptr().align_offset(align_of::<T>());
    if offset != 0 {
      // I assume we never enter here, but the assumption had to be tested!
      panic!("Unable to work from MMap, it is not well aligned!!");
    }
    let len = self.mmap.len() / size_of::<T>();
    BorrowedCIndex::new(self.depth, unsafe {
      &*slice_from_raw_parts(self.mmap.as_ptr() as *const T, len)
    })
  }
}

// EXPLICIT INDEX

// P.S.: to index a HPX sorted CSV file, the best thing to do would be to return an
// iterator on (hpx, starting_byte) and to create a `cds-bstree-file-readonly-rust` structure
// (from an already sorted iterator).
#[derive(Debug, PartialEq)]
pub struct OwnedCIndexExplicit<H: HHash, T: HCIndexValue> {
  depth: u8,
  entries: Vec<(H, T)>,
}
impl<H: HHash, T: HCIndexValue> OwnedCIndexExplicit<H, T> {
  pub fn new_unchecked(depth: u8, entries: Vec<(H, T)>) -> Self {
    Self { depth, entries }
  }
}
impl<H: HHash, T: HCIndexValue> HCIndex for OwnedCIndexExplicit<H, T> {
  type V = T;

  fn shape(&self) -> HCIndexShape {
    HCIndexShape::Explicit
  }

  fn depth(&self) -> u8 {
    self.depth
  }

  fn get(&self, hash: u64) -> Self::V {
    let hash = H::from_u64(hash);
    match self.entries.binary_search_by_key(&hash, |&(h, _)| h) {
      Ok(i) => self.entries[i].1,
      Err(i) => {
        if i == 0 {
          Self::V::zero()
        } else {
          self.entries[i - 1].1
        }
      }
    }
  }
  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_implicit_from_explicit(self.depth, self.entries.iter().cloned(), writer)
  }

  fn write_all_values_explicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_explicit(self.depth, self.entries.iter().map(|(k, v)| (k, v)), writer)
  }

  fn byte_size_explicit(&self) -> u64 {
    let size_of_hash = if self.depth <= 13 {
      size_of::<u32>()
    } else {
      size_of::<u64>()
    } as u64;
    self.entries.len() as u64 * (size_of::<T>() as u64 + size_of_hash)
  }
}

/// `BTreeMap` based structure: slow loading (tree creation), but fast queries.
// P.S.: to index a HPX sorted CSV file, the best thing to do would be to return an
// iterator on (hpx, starting_byte) and to create a `cds-bstree-file-readonly-rust` structur
// (from an already sorted iterator).
#[derive(Debug, PartialEq)]
pub struct OwnedCIndexExplicitBTree<H: HHash, T: HCIndexValue> {
  depth: u8,
  entries: BTreeMap<H, T>,
}
impl<H: HHash, T: HCIndexValue> OwnedCIndexExplicitBTree<H, T> {
  pub fn new(depth: u8, entries: Vec<(H, T)>) -> Self {
    Self {
      depth,
      entries: BTreeMap::from_iter(entries.into_iter()),
    }
  }
}
impl<H: HHash, T: HCIndexValue> HCIndex for OwnedCIndexExplicitBTree<H, T> {
  type V = T;

  fn shape(&self) -> HCIndexShape {
    HCIndexShape::Explicit
  }

  fn depth(&self) -> u8 {
    self.depth
  }

  fn get(&self, hash: u64) -> Self::V {
    let hash = H::from_u64(hash);
    match self.entries.get(&hash) {
      Some(v) => *v,
      None => match self.entries.range(..hash).next_back() {
        Some((_h, v)) => *v,
        None => Self::V::zero(),
      },
    }
  }

  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_implicit_from_explicit(
      self.depth,
      self.entries.iter().map(|(k, v)| (*k, *v)),
      writer,
    )
  }

  fn write_all_values_explicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    write_all_values_explicit(self.depth, self.entries.iter(), writer)
  }

  fn byte_size_explicit(&self) -> u64 {
    let size_of_hash = if self.depth <= 13 {
      size_of::<u32>()
    } else {
      size_of::<u64>()
    } as u64;
    self.entries.len() as u64 * (size_of::<T>() as u64 + size_of_hash)
  }
}

impl<'a, H: HHash, T: HCIndexValue> From<BorrowedCIndexExplicit<'a, H, T>>
  for OwnedCIndexExplicitBTree<H, T>
{
  fn from(value: BorrowedCIndexExplicit<'a, H, T>) -> Self {
    if value.depth <= 13 {
      const SIZE_OF_H: usize = size_of::<u32>();
      Self::new(
        value.depth,
        value
          .bytes
          .chunks(SIZE_OF_H + size_of::<T>())
          .map(|slice| {
            (
              H::from_u32(u32::from_le_bytes(slice[..SIZE_OF_H].try_into().unwrap())),
              T::from_le_bytes(&slice[SIZE_OF_H..].try_into().unwrap()),
            )
          })
          .collect(),
      )
    } else {
      const SIZE_OF_H: usize = size_of::<u64>();
      Self::new(
        value.depth,
        value
          .bytes
          .chunks(SIZE_OF_H + size_of::<T>())
          .map(|slice| {
            (
              H::from_u64(u64::from_le_bytes(slice[..SIZE_OF_H].try_into().unwrap())),
              T::from_le_bytes(&slice[SIZE_OF_H..].try_into().unwrap()),
            )
          })
          .collect(),
      )
    }
  }
}

#[derive(Debug)]
pub struct FitsMMappedCIndexExplicit<H: HHash, T: HCIndexValue> {
  fits_creation_date: Option<SystemTime>,
  indexed_file_name: Option<String>,
  indexed_file_len: Option<u64>,
  indexed_file_md5: Option<String>,
  indexed_file_last_modif_date: Option<SystemTime>,
  indexed_colname_lon: Option<String>,
  indexed_colname_lat: Option<String>,
  depth: u8,
  mmap: Mmap,
  _phantom_h: PhantomData<H>,
  _phantom_v: PhantomData<T>,
}
impl<H: HHash, T: HCIndexValue> FitsMMappedCIndexExplicit<H, T> {
  /// Private, only meant to be called from FITS reader.
  #[allow(clippy::too_many_arguments)]
  fn new(
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
      fits_creation_date,
      indexed_file_name,
      indexed_file_len,
      indexed_file_md5,
      indexed_file_last_modif_date,
      indexed_colname_lon,
      indexed_colname_lat,
      depth,
      mmap,
      _phantom_h: PhantomData,
      _phantom_v: PhantomData,
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
  pub fn get_hcindex(&self) -> BorrowedCIndexExplicit<'_, H, T> {
    BorrowedCIndexExplicit::new(self.depth, self.mmap.as_ref())
  }
}

#[derive(Debug, PartialEq)]
pub struct BorrowedCIndexExplicit<'a, H: HHash, T: HCIndexValue> {
  depth: u8,
  bytes: &'a [u8], // FITS content
  _type_h: PhantomData<H>,
  _type_t: PhantomData<T>,
}
impl<'a, H: HHash, T: HCIndexValue> BorrowedCIndexExplicit<'a, H, T> {
  fn new(depth: u8, bytes: &'a [u8]) -> Self {
    Self {
      depth,
      bytes,
      _type_h: PhantomData,
      _type_t: PhantomData,
    }
  }
  pub fn len(&self) -> usize {
    self.bytes.len()
      / if self.depth <= 13 {
        size_of::<u32>() + size_of::<T>()
      } else {
        size_of::<u64>() + size_of::<T>()
      }
  }
  fn get_h<S: HHash>(&self, index: usize) -> S {
    let h_size = size_of::<S>();
    let hv_size = h_size + size_of::<T>();
    let ptr = self.bytes.as_ptr();
    S::from_le_bytes(
      &unsafe { from_raw_parts(ptr.add(index * hv_size), h_size) }
        .try_into()
        .unwrap(),
    )
  }

  fn get_entry<S: HHash>(&self, index: usize) -> (S, T) {
    let h_size = size_of::<S>();
    let hv_size = h_size + size_of::<T>();
    let ptr = self.bytes.as_ptr();
    let sub_slice = unsafe { from_raw_parts(ptr.add(index * hv_size), hv_size) };
    let (z, v) = sub_slice.split_at(h_size);
    (
      S::from_le_bytes(&z.try_into().unwrap()),
      T::from_le_bytes(&v.try_into().unwrap()),
    )
  }
  /*fn binary_search_get(&self, hash: u64) -> Result<usize, usize> {
    if self.depth <= 13 {
      self.binary_search_get_u32(hash as u32)
    } else {
      self.binary_search_get_u64(hash)
    }
  }*/
  fn binary_search_get_u32(&self, hash: u32) -> Result<usize, usize> {
    self.binary_search_get_gen::<u32>(hash)
  }
  fn binary_search_get_u64(&self, hash: u64) -> Result<usize, usize> {
    self.binary_search_get_gen::<u64>(hash)
  }
  /// Here`S` is the `HHhash` type used in the storage
  fn binary_search_get_gen<S: HHash>(&self, hash: S) -> Result<usize, usize> {
    let h_size = size_of::<S>();
    let hv_size = h_size + size_of::<T>();
    let mut size = self.bytes.len() / hv_size;
    if size == 0 {
      return Err(0);
    }
    let mut base = 0usize;
    while size > 1 {
      let half = size >> 1;
      let mid = base + half;
      let cmp = self.get_h::<S>(mid).cmp(&hash);
      base = if cmp == Ordering::Greater { base } else { mid };
      size -= half;
    }
    let cmp = self.get_h::<S>(base).cmp(&hash);
    if cmp == Ordering::Equal {
      Ok(base)
    } else {
      let result = base + (cmp == Ordering::Less) as usize;
      Err(result)
    }
  }
}

impl<'a, H: HHash, T: HCIndexValue> HCIndex for BorrowedCIndexExplicit<'a, H, T> {
  type V = T;

  fn shape(&self) -> HCIndexShape {
    HCIndexShape::Explicit
  }

  fn depth(&self) -> u8 {
    self.depth
  }

  fn get(&self, hash: u64) -> Self::V {
    if self.depth <= 13 {
      match self.binary_search_get_u32(hash as u32) {
        Ok(i) => self.get_entry::<u32>(i).1,
        Err(i) => {
          if i == 0 {
            Self::V::zero()
          } else {
            self.get_entry::<u32>(i - 1).1
          }
        }
      }
    } else {
      match self.binary_search_get_u64(hash) {
        Ok(i) => self.get_entry::<u64>(i).1,
        Err(i) => {
          if i == 0 {
            Self::V::zero()
          } else {
            self.get_entry::<u64>(i - 1).1
          }
        }
      }
    }
  }

  fn write_all_values_implicit<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    if self.depth <= 13 {
      const SIZE_OF_H: usize = size_of::<u32>();
      write_all_values_implicit_from_explicit(
        self.depth,
        self
          .bytes
          .chunks(SIZE_OF_H + size_of::<Self::V>())
          .map(|slice| {
            (
              u32::from_le_bytes(slice[..SIZE_OF_H].try_into().unwrap()),
              T::from_le_bytes(&slice[SIZE_OF_H..].try_into().unwrap()),
            )
          }),
        writer,
      )
    } else {
      const SIZE_OF_H: usize = size_of::<u64>();
      write_all_values_implicit_from_explicit(
        self.depth,
        self
          .bytes
          .chunks(SIZE_OF_H + size_of::<Self::V>())
          .map(|slice| {
            (
              u64::from_le_bytes(slice[..SIZE_OF_H].try_into().unwrap()),
              T::from_le_bytes(&slice[SIZE_OF_H..].try_into().unwrap()),
            )
          }),
        writer,
      )
    }
  }

  fn write_all_values_explicit<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    writer.write_all(self.bytes).map(|()| self.bytes.len())
  }

  fn byte_size_explicit(&self) -> u64 {
    self.bytes.len() as u64
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[cfg(not(target_arch = "wasm32"))]
  #[test]
  fn testok_cindex_tofrom_fits() {
    let path = "cindex.fits";
    let depth = 10;
    let n = n_hash(depth);
    let values: Vec<u32> = (0..=n as u32).collect();
    let cindex = OwnedCIndex::new_unchecked(depth, values.into_boxed_slice());
    let mut indexex_file_md5 = [0_u8; 32];
    indexex_file_md5.copy_from_slice("0123456789ABCDEF0123456789ABCDEF".as_bytes());
    cindex
      .to_fits_file(
        path,
        HCIndexShape::Implicit,
        Some("filename"),
        Some(42),
        Some(indexex_file_md5),
        None,
        Some("RA"),
        Some("Dec"),
      )
      .unwrap();
    match FITSCIndex::from_fits_file(path) {
      Ok(FITSCIndex::ImplicitU32(hciprovider)) => {
        assert_eq!(
          hciprovider.get_indexed_file_name().map(|s| s.as_str()),
          Some("filename")
        );
        assert_eq!(hciprovider.get_indexed_file_len(), Some(42));
        assert_eq!(
          hciprovider.get_indexed_file_md5().map(|s| s.as_str()),
          Some("0123456789ABCDEF0123456789ABCDEF")
        );
        assert_eq!(
          hciprovider.get_indexed_colname_lon().map(|s| s.as_str()),
          Some("RA")
        );
        assert_eq!(
          hciprovider.get_indexed_colname_lat().map(|s| s.as_str()),
          Some("Dec")
        );

        let cindex2 = hciprovider.get_hcindex();
        // println!("{:?}", &cindex2);
        assert_eq!(cindex2, (&cindex).into());
        for h in 0..n {
          assert_eq!(1, cindex2.get_noncumulative_with_hash_at_index_depth(h));
        }
        for h in 0..n - 10 {
          assert_eq!(
            10,
            cindex2.get_noncumulative_with_range_at_index_depth(h..h + 10)
          );
        }
      }
      Ok(a) => {
        println!("{:?}", &a);
        panic!()
      }
      Err(e) => {
        println!("Err: {:?}", &e);
        panic!()
      }
    }
  }
}
