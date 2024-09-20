use std::{
  io,
  io::{BufRead, BufReader, Read, Seek},
  num::ParseIntError,
  slice::ChunksExact,
  str::{self, FromStr},
};

use byteorder::{BigEndian, ReadBytesExt};
use log::{debug, warn};

use super::{
  super::skymap::{ImplicitSkyMapArray, SkyMapEnum},
  error::FitsError,
  gz::{is_gz, uncompress},
  keywords::{
    CoordSys, FitsCard, IndexSchema, Nside, Order, Ordering, SkymapKeywords, SkymapKeywordsMap,
    TForm1, TType1,
  },
};
use crate::{depth, is_nside, n_hash};

/// We expect the FITS file to be a BINTABLE containing a map.
/// [Here](https://gamma-astro-data-formats.readthedocs.io/en/latest/skymaps/healpix/index.html)
/// a description of the format.
/// We so far implemented a subset of the format only:
/// * `INDXSCHM= 'IMPLICIT'`
/// * `ORDERING= 'NESTED  '`
/// To be fast (in execution and development), we start by a non-flexible approach in which we
/// expect the BINTABLE extension to contains:
/// ```bash
/// XTENSION= 'BINTABLE'           / binary table extension                         
/// BITPIX  =                    8 / array data type                                
/// NAXIS   =                    2 / number of array dimensions
/// NAXIS1  =                    ?? / length of dimension 1                          
/// NAXIS2  =                   ?? / length of dimension 2                          
/// PCOUNT  =                    0 / number of group parameters                     
/// GCOUNT  =                    1 / number of groups                               
/// TFIELDS =                   ?? / number of table fields
/// TTYPE1  = 'XXX'       // SHOULD STARS WITH 'PROB', else WARNING                                                            
/// TFORM1  = 'XXX'       // MUST CONTAINS D (f64) or E (f32)                                                            
/// TUNIT1  = 'pix-1    '
/// TTYPE2  = ???                                                         
/// TFORM2  = ???                                                            
/// ...
/// MOC     =                    T                                                  
/// PIXTYPE = 'HEALPIX '           / HEALPIX pixelisation                           
/// ORDERING= 'NESTED  '           / Pixel ordering scheme: RING, NESTED, or NUNIQ  
/// COORDSYS= 'C       '  // WARNING if not found
/// NSIDE    =                  ?? / MOC resolution (best nside)
///  or
/// ORDER    =                  ?? / MOC resolution (best order), superseded by NSIDE
///                                / (because NSIDE which are not power of 2 are possible in RING)
/// INDXSCHM= 'IMPLICIT'           / Indexing: IMPLICIT or EXPLICIT
/// ...
/// END
/// ```
///
/// # Params
/// * `reader`: the reader over the FITS content
///
/// # Info
///   Supports gz input stream
///
pub fn from_fits_skymap<R: Read + Seek>(mut reader: BufReader<R>) -> Result<SkyMapEnum, FitsError> {
  if is_gz(&mut reader)? {
    // Probably need to build an explicit map:
    //   highly compressed skymaps are IMPLICIT with a lot of zeros.
    from_fits_skymap_internal(uncompress(reader))
  } else {
    from_fits_skymap_internal(reader)
  }
}

pub fn from_fits_skymap_internal<R: BufRead>(mut reader: R) -> Result<SkyMapEnum, FitsError> {
  let mut header_block = [b' '; 2880];
  consume_primary_hdu(&mut reader, &mut header_block)?;
  // Read the extension HDU
  let mut it80 = next_36_chunks_of_80_bytes(&mut reader, &mut header_block)?;
  // See Table 10 and 17 in https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
  check_keyword_and_val(it80.next().unwrap(), b"XTENSION", b"'BINTABLE'")?;
  check_keyword_and_val(it80.next().unwrap(), b"BITPIX  ", b"8")?;
  check_keyword_and_val(it80.next().unwrap(), b"NAXIS  ", b"2")?;
  let n_bytes_per_row = check_keyword_and_parse_uint_val::<u64>(it80.next().unwrap(), b"NAXIS1  ")?;
  let n_rows = check_keyword_and_parse_uint_val::<u64>(it80.next().unwrap(), b"NAXIS2 ")?;
  check_keyword_and_val(it80.next().unwrap(), b"PCOUNT  ", b"0")?;
  check_keyword_and_val(it80.next().unwrap(), b"GCOUNT  ", b"1")?;
  let n_cols = check_keyword_and_parse_uint_val::<u64>(it80.next().unwrap(), b"TFIELDS ")?;
  let colname_1 = TType1::parse_value(it80.next().unwrap())?;
  debug!("Skymap column name: {}", colname_1.get());
  let coltype_1 = TForm1::parse_value(it80.next().unwrap())?;
  debug!("Skymap column type: {}", coltype_1.to_fits_value());

  // nbits = |BITPIX|xGCOUNTx(PCOUNT+NAXIS1xNAXIS2x...xNAXISn)
  // In our case (bitpix = Depends on TForm, GCOUNT = 1, PCOUNT = 0) => nbytes = n_cells * size_of(T)
  // let data_size n_bytes as usize * n_cells as usize; // N_BYTES ok since BITPIX = 8
  // Read MOC keywords
  let mut skymap_kws = SkymapKeywordsMap::new();
  'hr: loop {
    for kw_record in &mut it80 {
      match SkymapKeywords::is_skymap_kw(kw_record)? {
        Some(kw) => {
          if let Some(previous_mkw) = skymap_kws.insert(kw) {
            // A FITS keyword MUST BE uniq (I may be more relax here, taking the last one and not complaining)
            warn!(
              "Keyword '{}' found more than once in a same HDU! We use the first occurrence.",
              previous_mkw.keyword_str()
            );
            skymap_kws.insert(previous_mkw);
          }
        }
        None => {
          if &kw_record[0..4] == b"END " {
            break 'hr;
          } else {
            debug!("Ignored FITS card: {}", unsafe {
              str::from_utf8_unchecked(kw_record)
            })
          }
        }
      }
    }
    // Read next 2880 bytes
    it80 = next_36_chunks_of_80_bytes(&mut reader, &mut header_block)?;
  }
  // Check constant header params
  skymap_kws.check_pixtype()?; // = HEALPIX
  skymap_kws.check_ordering(Ordering::Nested)?; // So far we support only 'NESTED'
  skymap_kws.check_coordsys(CoordSys::Cel)?; // So far we support only Celestial coordinates (not Galactic)
  skymap_kws.check_index_schema(IndexSchema::Implicit)?; // So far we support only 'IMLPLICIT'
  skymap_kws.check_firstpix(0)?;
  // Check number of columns
  // - we so far support only map having a single column of values
  if n_cols != 1 {
    return Err(FitsError::UnexpectedValue {
      keyword: String::from("TFIELDS"),
      expected: 1.to_string(),
      actual: n_cols.to_string(),
    });
  }
  // Check whether lastpix
  let depth = match skymap_kws.get::<Order>() {
    Some(SkymapKeywords::Order(order)) => Ok(order.get()),
    None => match skymap_kws.get::<Nside>() {
      Some(SkymapKeywords::Nside(nside)) => {
        let nside = nside.get();
        if is_nside(nside) {
          Ok(depth(nside))
        } else {
          Err(FitsError::new_custom(format!(
            "Nside is not valid (to be used in nested mode at least): {}",
            nside
          )))
        }
      }
      None => Err(FitsError::new_custom(String::from(
        "Both keywords 'ORDER' and 'NSIDE' are missing!",
      ))),
      _ => unreachable!(),
    },
    _ => unreachable!(),
  }?;
  let n_hash = n_hash(depth);
  skymap_kws.check_lastpix(n_hash)?;
  // Check whether TForm compatible with N bytes per row
  let n_hash_2 = coltype_1.n_pack() as u64 * n_rows;
  if n_hash != n_hash_2 {
    return Err(FitsError::new_custom(format!(
      "Number of elements {} do not match number of HEALPix cells {}",
      n_hash_2, n_hash
    )));
  }
  // Check n_bytes_per_row
  if n_bytes_per_row != coltype_1.n_bytes() as u64 {
    return Err(FitsError::new_custom(format!(
      "Number of bytes per row {} do not match TFORM1 = {}",
      n_bytes_per_row,
      coltype_1.to_fits_value()
    )));
  }

  // Read data
  match coltype_1 {
    TForm1::B(_) => (0..n_hash)
      .map(|_| reader.read_u8())
      .collect::<Result<Vec<u8>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64U8(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
    TForm1::I(_) => (0..n_hash)
      .map(|_| reader.read_i16::<BigEndian>())
      .collect::<Result<Vec<i16>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64I16(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
    TForm1::J(_) => (0..n_hash)
      .map(|_| reader.read_i32::<BigEndian>())
      .collect::<Result<Vec<i32>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64I32(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
    TForm1::K(_) => (0..n_hash)
      .map(|_| reader.read_i64::<BigEndian>())
      .collect::<Result<Vec<i64>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64I64(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
    TForm1::E(_) => (0..n_hash)
      .map(|_| reader.read_f32::<BigEndian>())
      .collect::<Result<Vec<f32>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64F32(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
    TForm1::D(_) => (0..n_hash)
      .map(|_| reader.read_f64::<BigEndian>())
      .collect::<Result<Vec<f64>, io::Error>>()
      .map(|v| SkyMapEnum::ImplicitU64F64(ImplicitSkyMapArray::new(depth, v.into_boxed_slice()))),
  }
  .map_err(FitsError::Io)
}

const VALUE_INDICATOR: &[u8; 2] = b"= ";

/// # Params
/// - `header_block`: re-usable header block used to avoid multiple allocations
fn consume_primary_hdu<R: BufRead>(
  reader: &mut R,
  header_block: &mut [u8; 2880],
) -> Result<(), FitsError> {
  let mut chunks_of_80 = next_36_chunks_of_80_bytes(reader, header_block)?;
  // SIMPLE = 'T' => file compliant with the FITS standard
  check_keyword_and_val(chunks_of_80.next().unwrap(), b"SIMPLE ", b"T")?;
  chunks_of_80.next().unwrap(); // Do not check for BITPIX (we expect an empty header)
                                // NAXIS = 0 => we only support FITS files with no data in the primary HDU
  check_keyword_and_val(chunks_of_80.next().unwrap(), b"NAXIS ", b"0")?;
  // Ignore possible additional keywords
  while !contains_end(&mut chunks_of_80) {
    // Few chances to enter here (except if someone had a lot of things to say in the header)
    chunks_of_80 = next_36_chunks_of_80_bytes(reader, header_block)?;
  }
  Ok(())
}

fn next_36_chunks_of_80_bytes<'a, R: BufRead>(
  reader: &'a mut R,
  header_block: &'a mut [u8; 2880],
) -> Result<ChunksExact<'a, u8>, FitsError> {
  reader.read_exact(header_block)?;
  Ok(header_block.chunks_exact(80))
}

fn contains_end<'a, I: Iterator<Item = &'a [u8]>>(chunks_of_80: &'a mut I) -> bool {
  for kw_rc in chunks_of_80 {
    debug_assert_eq!(kw_rc.len(), 80);
    if &kw_rc[0..4] == b"END " {
      return true;
    }
  }
  false
}

pub(super) fn check_keyword_and_val(
  keyword_record: &[u8],
  expected_kw: &[u8],
  expected_val: &[u8],
) -> Result<(), FitsError> {
  check_expected_keyword(keyword_record, expected_kw)?;
  check_for_value_indicator(keyword_record)?;
  check_expected_value(keyword_record, expected_val)
}

fn check_keyword_and_parse_uint_val<T>(
  keyword_record: &[u8],
  expected_kw: &[u8],
) -> Result<T, FitsError>
where
  T: Into<u64> + FromStr<Err = ParseIntError>,
{
  check_expected_keyword(keyword_record, expected_kw)?;
  check_for_value_indicator(keyword_record)?;
  parse_uint_val::<T>(keyword_record)
}

#[allow(dead_code)]
pub(super) fn check_keyword_and_get_str_val<'a>(
  keyword_record: &'a [u8],
  expected_kw: &[u8],
) -> Result<&'a str, FitsError> {
  check_expected_keyword(keyword_record, expected_kw)?;
  check_for_value_indicator(keyword_record)?;
  // We go unsafe because FITS headers are not supposed to contain non-ASCII chars
  get_str_val_no_quote(keyword_record).map(|bytes| unsafe { str::from_utf8_unchecked(bytes) })
}

pub(super) fn check_expected_keyword(
  keyword_record: &[u8],
  expected: &[u8],
) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80); // length of a FITS keyword-record
  debug_assert!(expected.len() <= 8); // length of a FITS keyword
  if &keyword_record[..expected.len()] == expected {
    Ok(())
  } else {
    // We know what we put in it, so unsafe is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) }.trim_end());
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&keyword_record[..expected.len()])
      .trim_end()
      .to_string();
    // panic!("Ecpected: {}, Actual: {}", expected, String::from_utf8_lossy(&src[..]));
    Err(FitsError::UnexpectedKeyword { expected, actual })
  }
}

pub(super) fn check_for_value_indicator(keyword_record: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80); // length of a FITS keyword-record
  if get_value_indicator(keyword_record) == VALUE_INDICATOR {
    Ok(())
  } else {
    let keyword_record = String::from_utf8_lossy(keyword_record)
      .trim_end()
      .to_string();
    Err(FitsError::ValueIndicatorNotFound { keyword_record })
  }
}

pub(super) fn get_keyword(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[..8]
}
pub(super) fn get_value_indicator(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[8..10]
}
pub(super) fn get_value(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[10..]
}
pub(super) fn get_left_trimmed_value(keyword_record: &[u8]) -> &[u8] {
  get_value(keyword_record).trim_ascii_start()
}

pub(super) fn check_expected_value(
  keyword_record: &[u8],
  expected: &[u8],
) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80); // length of a FITS keyword-record
  let src = get_value(keyword_record);
  let lt_src = src.trim_ascii_start();
  if lt_src.len() >= expected.len() && &lt_src[..expected.len()] == expected {
    Ok(())
  } else {
    let keyword = String::from_utf8_lossy(&src[0..8]).trim_end().to_string();
    // We know what we put in it, so unsae is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) });
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&lt_src[..expected.len()]).to_string();
    Err(FitsError::UnexpectedValue {
      keyword,
      expected,
      actual,
    })
  }
}

/// We know that the expected value does not contains a simple quote.
/// A trim_end is applied so that the result does not contain leading or trailing spaces.
pub(super) fn get_str_val_no_quote(keyword_record: &[u8]) -> Result<&[u8], FitsError> {
  let mut it = get_left_trimmed_value(keyword_record).split_inclusive(|c| *c == b'\'');
  if let Some([b'\'']) = it.next() {
    if let Some([subslice @ .., b'\'']) = it.next() {
      return Ok(subslice.trim_ascii_end());
    }
  }
  let keyword_record = String::from_utf8_lossy(keyword_record)
    .trim_end()
    .to_string();
  Err(FitsError::StringValueNotFound { keyword_record })
}

pub(super) fn parse_uint_val<T>(keyword_record: &[u8]) -> Result<T, FitsError>
where
  T: Into<u64> + FromStr<Err = ParseIntError>,
{
  let src = get_left_trimmed_value(keyword_record);
  let to = index_of_last_digit(src);
  if to == 0 {
    let keyword_record = String::from_utf8_lossy(keyword_record)
      .trim_end()
      .to_string();
    Err(FitsError::UintValueNotFound { keyword_record })
  } else {
    // we go unsafe and unwrap because we already tested for regular digits
    let str_val = unsafe { str::from_utf8_unchecked(&src[..to]) };
    str_val.parse::<T>().map_err(|e| FitsError::WrongUintValue {
      context: str_val.to_string(),
      err: e,
    })
  }
}

pub(super) fn index_of_last_digit(src: &[u8]) -> usize {
  for (i, c) in src.iter().enumerate() {
    if !c.is_ascii_digit() {
      return i;
    }
  }
  src.len()
}
