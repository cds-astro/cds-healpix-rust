use std::io::{BufRead, BufReader, Read, Seek};

use flate2::read::GzDecoder;

const GZ_MAGIC_NUM: [u8; 2] = [0x1F, 0x8B];
const GZ_MAGIC_NUM_LEN: usize = GZ_MAGIC_NUM.len();

pub fn is_gz<R: Read + Seek>(reader: &mut BufReader<R>) -> Result<bool, std::io::Error> {
  let mut gz_magic_bytes = [0u8; 2];
  reader.read_exact(&mut gz_magic_bytes)?;
  reader.seek_relative(-(GZ_MAGIC_NUM_LEN as i64))?;
  Ok(gz_magic_bytes == GZ_MAGIC_NUM)
}

/// Returns an object implementing `BufRead` and decompressing on-th-fly the input `BufRead`.
pub fn uncompress<R: BufRead>(reader: R) -> BufReader<GzDecoder<R>> {
  BufReader::new(GzDecoder::new(reader))
}
