use std::io::{Error as IoError, Write};

use num_traits::ToBytes;

use crate::{depth_from_n_hash_unsafe, n_hash, nested::map::skymap::SkyMapValue, nside};

use super::error::FitsError;

// Add read map from fits!!

/// Write the values in a FITS HEALPix map.
///
/// # Params
/// * `values`: array of values, the index of the value correspond to the HEALPix cell the value is
///   associated with.
pub fn write_implicit_skymap_fits<R: Write, T: SkyMapValue>(
  mut writer: R,
  values: &[T],
) -> Result<(), FitsError> {
  // We could have called 'write_implicit_skymap_fits_from_it' but we
  // wanted to avoid counting the number of elements in the iterator
  // since we already know it.
  // In the end, using monomorphisation, the binary code is duplicated so...
  let n_cells = values.len() as u64;
  let depth = depth_from_n_hash_unsafe(values.len() as u64);
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
  write_final_padding(writer, n_cells as usize * T::FITS_NAXIS1 as usize)
}

pub fn write_implicit_skymap_fits_from_it<R, I>(
  mut writer: R,
  depth: u8,
  it: I,
) -> Result<(), FitsError>
where
  R: Write,
  I: Iterator,
  I::Item: SkyMapValue,
{
  let n_cells = n_hash(depth);
  write_skymap_fits_header::<_, I::Item>(&mut writer, depth)?;
  let mut n = 0;
  for value in it {
    writer.write_all(value.to_be_bytes().as_ref())?;
    n += 1;
  }
  if n != n_cells {
    return Err(FitsError::new_custom(format!(
      "Wrong number of HEALPix cells writen at depth {}. Expected: {}. Actual: {}.",
      depth, n_cells, n
    )));
  }
  write_final_padding(writer, n_cells as usize * I::Item::FITS_NAXIS1 as usize)
}

/// Possible add blanks at the end of the FITS file to complete the last
/// 2880 bytes block.
pub(crate) fn write_final_padding<R: Write>(
  writer: R,
  n_bytes_already_written: usize,
) -> Result<(), FitsError> {
  write_final_padding_ioerr(writer, n_bytes_already_written).map_err(FitsError::Io)
}

pub(crate) fn write_final_padding_ioerr<R: Write>(
  mut writer: R,
  n_bytes_already_written: usize,
) -> Result<(), IoError> {
  let mod2880 = n_bytes_already_written % 2880;
  if mod2880 != 0 {
    writer.write_all(&vec![0_u8; 2880 - mod2880])
  } else {
    Ok(())
  }
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

pub(crate) fn write_uint_mandatory_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: u64) {
  const VALUE_INDICATOR: &[u8; 2] = b"= ";
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val = val.to_string();
  let val_bytes = val.as_bytes();
  dest[30 - val_bytes.len()..30].copy_from_slice(val_bytes);
}

pub(crate) fn write_str_mandatory_keyword_record(dest: &mut [u8], keyword: &[u8; 8], s: &str) {
  const VALUE_INDICATOR: &[u8; 2] = b"= ";
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val = format!("'{:<8}'", s);
  let val_bytes = val.as_bytes();
  dest[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
}

pub(crate) fn write_primary_hdu<R: Write>(writer: &mut R) -> Result<(), FitsError> {
  write_primary_hdu_ioerr(writer).map_err(FitsError::Io)
}
pub(crate) fn write_primary_hdu_ioerr<R: Write>(writer: &mut R) -> Result<(), IoError> {
  let mut header_block = [b' '; 2880];
  header_block[0..30].copy_from_slice(b"SIMPLE  =                    T");
  header_block[80..110].copy_from_slice(b"BITPIX  =                    8");
  header_block[160..190].copy_from_slice(b"NAXIS   =                    0");
  header_block[240..270].copy_from_slice(b"EXTEND  =                    T");
  header_block[320..323].copy_from_slice(b"END");
  writer.write_all(&header_block[..])
}

fn write_skymap_fits_header<R: Write, T: SkyMapValue>(
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
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS1  ", T::FITS_NAXIS1 as u64);
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", n_cells);
  it.next().unwrap()[0..30].copy_from_slice(b"PCOUNT  =                    0");
  it.next().unwrap()[0..30].copy_from_slice(b"GCOUNT  =                    1");
  it.next().unwrap()[0..30].copy_from_slice(b"TFIELDS =                    1");
  it.next().unwrap()[0..20].copy_from_slice(b"TTYPE1  = 'T       '");
  write_str_mandatory_keyword_record(it.next().unwrap(), b"TFORM1  ", T::FITS_TFORM);
  it.next().unwrap()[0..20].copy_from_slice(b"PIXTYPE = 'HEALPIX '");
  it.next().unwrap()[0..20].copy_from_slice(b"ORDERING= 'NESTED  '");
  it.next().unwrap()[0..20].copy_from_slice(b"COORDSYS= 'C       '");
  it.next().unwrap()[0..20].copy_from_slice(b"EXTNAME = 'xtension'");
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"NSIDE   ", nside as u64);
  it.next().unwrap()[0..30].copy_from_slice(b"FIRSTPIX=                    0");
  write_uint_mandatory_keyword_record(it.next().unwrap(), b"LASTPIX ", n_cells - 1);
  it.next().unwrap()[0..20].copy_from_slice(b"INDXSCHM= 'IMPLICIT'");
  it.next().unwrap()[0..20].copy_from_slice(b"OBJECT  = 'FULLSKY '");
  // it.next().unwrap()[0..28].copy_from_slice(b"CREATOR = 'CDS HEALPix Rust'");
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
