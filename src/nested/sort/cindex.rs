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
  cmp::Ordering,
  fs::File,
  io::{BufWriter, Error as IoError, Write},
  marker::PhantomData,
  mem::align_of,
  ops::{AddAssign, Range, Sub},
  path::Path,
  ptr::slice_from_raw_parts,
  str,
  time::SystemTime,
};

use chrono::{DateTime, SecondsFormat, Utc};
use log::debug;
#[cfg(feature = "memmap")]
use memmap2::{Mmap, MmapOptions};
use num_traits::{ToBytes, Zero};

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
  },
};

/// Trait defining the type of (cumulative) value a Healpix Cumulative Index may contains.
pub trait HCIndexValue:
  Sized + Zero + Sub<Output = Self> + AddAssign + ToBytes + Clone + Copy
{
  const FITS_DATATYPE: &'static str;
}

impl HCIndexValue for u32 {
  const FITS_DATATYPE: &'static str = "u32";
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

  /// Returns the HEALPix depth of the index.
  fn depth(&self) -> u8;

  /// Returns the cumulative value associated to the given index.
  /// The starting cumulative value associated to a HEALPix cell of index `i` is obtained by `get(i)`
  /// while the ending cumulative value is obtained by `get(i+1)`.
  /// It means that the largest index equals the number of HEALPix cells at the index depth, and the
  /// number of stored values equals thi number plus one.
  ///
  /// # Panics
  /// If `index` larger than `n_hash(depth)`.
  fn get(&self, index: usize) -> Self::V;

  /// For the firs serialization. Must write the `n_hash(depth) + 1` values, in ascnding order,
  /// and in little-endian, in the given writer.
  /// Exactly `(n_hash(depth) + 1) * size_of::<Self::V>` bytes must be written.
  fn write_all_values<W: Write>(&self, writer: W) -> Result<usize, IoError>;

  /// Returns both the starting and the ending values associated with the given input `hash` at the index `depth`.
  fn get_with_hash_at_index_depth(&self, hash: u64) -> Range<Self::V> {
    self.get_with_range_at_index_depth(hash..hash + 1)
  }
  /// Returns both the starting and the ending values associated with the given input `range` at the index `depth`.
  fn get_with_range_at_index_depth(&self, range: Range<u64>) -> Range<Self::V> {
    self.get(range.start as usize)..self.get(range.end as usize)
  }

  /// Returns the noncumulative value associated with the given input `hash` at the index `depth`.
  fn get_noncumulative_with_hash_at_index_depth(&self, hash: u64) -> Self::V {
    self.get_noncumulative_with_range_at_index_depth(hash..hash + 1)
  }
  /// Returns the noncumulative value associated with the given input `range` at the index `depth`.
  fn get_noncumulative_with_range_at_index_depth(&self, range: Range<u64>) -> Self::V {
    self.get(range.end as usize) - self.get(range.start as usize)
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
    let indexed_file_last_modif_date =
      indexed_file_last_modif_date.map(DateTime::<Utc>::from);
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
      .write_all_values(&mut writer)
      .map_err(FitsError::Io)
      .and_then(|n_bytes_written| write_final_padding(writer, n_bytes_written))
  }

  #[allow(clippy::too_many_arguments)]
  fn to_fits_file<P: AsRef<Path>>(
    &self,
    path: P,
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
        indexed_file_name,
        indexed_file_len,
        indexed_file_md5,
        indexed_file_last_modif_date,
        indexed_colname_lon,
        indexed_colname_lat,
      )
    })
  }
}

// Make an 'incomplete' HCI for which the first (or the last) value is provided so that on disk
// the binary size will be 12 * a power of two, which is buffer, cache, ... friendly.

#[derive(Debug, PartialEq)]
pub struct OwnedCIndex<T: HCIndexValue> {
  depth: u8,
  values: Box<[T]>,
}
impl<T: HCIndexValue> OwnedCIndex<T> {
  pub fn new_unchecked(depth: u8, values: Box<[T]>) -> Self {
    Self { depth, values }
  }
}
impl<T: HCIndexValue> HCIndex for OwnedCIndex<T> {
  type V = T;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn get(&self, index: usize) -> Self::V {
    self.values[index]
  }

  fn write_all_values<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    /* Replaced by unsafe by arguably faster code
    for v in self.values {
      writer.write_all(v.to_le_bytes().as_ref())?;
    }*/
    write_all_values(self.values.as_ref(), writer)
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
  fn depth(&self) -> u8 {
    self.depth
  }
  fn get(&self, index: usize) -> T {
    self.values[index]
  }
  fn write_all_values<W: Write>(&self, writer: W) -> Result<usize, IoError> {
    /* Replaced by unsafe by arguably faster code
    for v in self.values {
      writer.write_all(v.to_le_bytes().as_ref())?;
    }*/
    write_all_values(self.values, writer)
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

fn write_all_values<T: HCIndexValue, W: Write>(
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

impl<'a, T: HCIndexValue> BorrowedCIndex<'a, T> {
  /// For each HEALPix cell at the given depth, `values` contains a cumutative quantity.
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
  U32(FitsMMappedCIndex<u32>),
  U64(FitsMMappedCIndex<u64>),
  F32(FitsMMappedCIndex<f32>),
  F64(FitsMMappedCIndex<f64>),
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
    let mut indxcov_found = false;
    let mut dtendian_found = false;
    let mut hpxoder: Option<u8> = None;
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
          .map(|()| indxschm_found = true),
        b"INDXCOV " => check_keyword_and_str_val(kw_record, b"INDXCOV ", b"FULLSKY")
          .map(|()| indxcov_found = true),
        b"HPXORDER" => parse_uint_val::<u8>(kw_record).map(|v| hpxoder = Some(v)),
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
    let n_bytes_data = n_bytes_per_val * n_rows;
    let mmap = unsafe {
      MmapOptions::new()
        .offset(2880)
        .len(n_bytes_data as usize)
        .map(&file)
        .map_err(FitsError::Io)?
    };
    match datatype.as_deref() {
      Some(u32::FITS_DATATYPE) => Ok(FITSCIndex::U32(FitsMMappedCIndex::new(
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
      Some(u64::FITS_DATATYPE) => Ok(FITSCIndex::U64(FitsMMappedCIndex::new(
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
      Some(f32::FITS_DATATYPE) => Ok(FITSCIndex::F32(FitsMMappedCIndex::new(
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
      Some(f64::FITS_DATATYPE) => Ok(FITSCIndex::F64(FitsMMappedCIndex::new(
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

/// The result of reading a FITS file, containing a memory map on the data, and from which we can obtain
/// an actual Healpix Cumulative Index object.
#[derive(Debug)]
pub struct FitsMMappedCIndex<T: HCIndexValue> {
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
impl<T: HCIndexValue> FitsMMappedCIndex<T> {
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
      panic!("Unablle to work from MMap, it is noot well aligned!!");
    }
    let len = self.mmap.len() / size_of::<T>();
    BorrowedCIndex::new(self.depth, unsafe {
      &*slice_from_raw_parts(self.mmap.as_ptr() as *const T, len)
    })
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
        Some("filename"),
        Some(42),
        Some(indexex_file_md5),
        None,
        Some("RA"),
        Some("Dec"),
      )
      .unwrap();
    match FITSCIndex::from_fits_file(path) {
      Ok(FITSCIndex::U32(hciprovider)) => {
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
