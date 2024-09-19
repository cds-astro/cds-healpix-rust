use std::io::Write;

use num_traits::ToBytes;

use crate::{n_hash, nside};

use super::error::FitsError;

/// Trait marking the type of the values written in the FITS map.
pub trait FitsCompatibleValue: ToBytes {
  /// Size, in bytes, of a value.
  fn fits_naxis1() -> u8;
  fn fits_tform() -> &'static str;
}

impl FitsCompatibleValue for u8 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "B"
  }
}

impl FitsCompatibleValue for i16 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "I"
  }
}

impl FitsCompatibleValue for i32 {
  fn fits_naxis1() -> u8 {
    4
  }

  fn fits_tform() -> &'static str {
    "J"
  }
}

impl FitsCompatibleValue for i64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "K"
  }
}

impl FitsCompatibleValue for u32 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "J"
  }
}

impl FitsCompatibleValue for u64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "K"
  }
}

impl FitsCompatibleValue for f32 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "E"
  }
}

impl FitsCompatibleValue for f64 {
  fn fits_naxis1() -> u8 {
    size_of::<Self>() as u8
  }

  fn fits_tform() -> &'static str {
    "D"
  }
}

// Add read map from fits!!

/// Write the values in a FITS HEALPix map.
///
/// # Params
/// * `values`: array of values, the index of the value correspond to the HEALPix cell the value is
/// associated with.
pub fn write_implicit_skymap_fits<R: Write, T: FitsCompatibleValue>(
  mut writer: R,
  values: &[T],
) -> Result<(), FitsError> {
  let n_cells = values.len() as u64;
  let depth = (values.len() / 12).trailing_zeros() as u8 >> 1;
  if n_cells != n_hash(depth) {
    return Err(FitsError::new_custom(format!(
      "Number of cell {} not compatible with an HEALPix depth of {}. Expected: {}.",
      n_cells,
      depth,
      n_hash(depth)
    )));
  }
  write_skymap_fits_header::<_, T>(&mut writer, depth)?;
  for value in values {
    writer.write_all(value.to_be_bytes().as_ref())?;
  }
  // Complete FITS block of 2880 bytes
  let mod2880 = (n_cells as usize * T::fits_naxis1() as usize) % 2880;
  if mod2880 != 0 {
    writer.write_all(&vec![0_u8; 2880 - mod2880])?;
  }
  Ok(())
}

#[allow(dead_code)]
pub(crate) fn write_uint_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: u64) {
  write_keyword_record(dest, keyword, &val.to_string())
}

#[allow(dead_code)]
pub(crate) fn write_str_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: &str) {
  write_keyword_record(dest, keyword, &format!("'{}'", val))
}

/// # Params
/// - `dest` destination, must contains 80 bytes
/// - `value_part` is string, must be already quote: `'str_value'`
pub(crate) fn write_keyword_record(dest: &mut [u8], keyword: &[u8; 8], value_part: &str) {
  const VALUE_INDICATOR: &[u8; 2] = b"= ";
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val_bytes = value_part.as_bytes();
  dest[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
}

fn write_uint_mandatory_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: u64) {
  const VALUE_INDICATOR: &[u8; 2] = b"= ";
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val = val.to_string();
  let val_bytes = val.as_bytes();
  dest[30 - val_bytes.len()..30].copy_from_slice(val_bytes);
}

fn write_str_mandatory_keyword_record(dest: &mut [u8], keyword: &[u8; 8], s: &str) {
  const VALUE_INDICATOR: &[u8; 2] = b"= ";
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val = format!("'{:>8}'", s);
  let val_bytes = val.as_bytes();
  dest[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
}

fn write_primary_hdu<R: Write>(writer: &mut R) -> Result<(), FitsError> {
  let mut header_block = [b' '; 2880];
  header_block[0..30].copy_from_slice(b"SIMPLE  =                    T");
  header_block[80..110].copy_from_slice(b"BITPIX  =                    8");
  header_block[160..190].copy_from_slice(b"NAXIS   =                    0");
  header_block[240..270].copy_from_slice(b"EXTEND  =                    T");
  header_block[320..323].copy_from_slice(b"END");
  writer.write_all(&header_block[..]).map_err(FitsError::Io)
}

fn write_skymap_fits_header<R: Write, T: FitsCompatibleValue>(
  mut writer: R,
  depth: u8,
) -> Result<(), FitsError> {
  let nside = nside(depth);
  let n_cells = n_hash(depth);

  write_primary_hdu(&mut writer)?;
  let mut header_block = [b' '; 2880];
  let mut it = header_block.chunks_mut(80);
  // Write BINTABLE specific keywords in the buffer
  it.next().unwrap()[0..20].copy_from_slice(b"XTENSION= 'BINTABLE'");
  it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8");
  it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2");
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS1  ", T::fits_naxis1() as u64);
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", n_cells);
  it.next().unwrap()[0..30].copy_from_slice(b"PCOUNT  =                    0");
  it.next().unwrap()[0..30].copy_from_slice(b"GCOUNT  =                    1");
  it.next().unwrap()[0..30].copy_from_slice(b"TFIELDS =                    1");
  it.next().unwrap()[0..20].copy_from_slice(b"TTYPE1  = 'T       '");
  write_str_mandatory_keyword_record(it.next().unwrap(), b"TFORM1  ", T::fits_tform());
  it.next().unwrap()[0..20].copy_from_slice(b"PIXTYPE = 'HEALPIX '");
  it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
  it.next().unwrap()[0..20].copy_from_slice(b"EXTNAME = 'xtension'");
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NSIDE   ", nside as u64);
  it.next().unwrap()[0..30].copy_from_slice(b"FIRSTPIX=                    0");
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"LASTPIX ", n_cells);
  it.next().unwrap()[0..20].copy_from_slice(b"INDXSCHM= 'IMPLICIT'");
  it.next().unwrap()[0..20].copy_from_slice(b"OBJECT  = 'FULLSKY '");
  it.next().unwrap()[0..26].copy_from_slice(b"CREATOR = 'CDS HEALPix Rust'");
  it.next().unwrap()[0..3].copy_from_slice(b"END");
  // Do write the header
  writer.write_all(&header_block[..]).map_err(FitsError::Io)
}
