use std::{
  error::Error,
  fs::File,
  io::{BufRead, BufReader, BufWriter},
  path::PathBuf,
};

use clap::Args;
use hpxlib::nested::{
  get, n_hash,
  sort::{
    cindex::{HCIndex, OwnedCIndex},
    get_hpx,
  },
};

use crate::input::csv::Csv;

/// Create an index on an HEALPix NESTED sorted CSV file to then quickly retrieve rows in a given HEALPix cell.
#[derive(Debug, Args)]
pub struct HealpixCumulIndex {
  #[command(flatten)]
  csv: Csv,
  /// Path of the output FITS file containing the HEALPix Cumulative Index.
  #[clap(short = 'o', long = "out", value_name = "FILE")]
  output: PathBuf,
  /// Field number of the longitude(in degrees) used to compute the HEALPix number, starting from 1.
  #[clap(short = 'l', long, value_name = "FIELD", default_value_t = 1)]
  lon: usize,
  /// Field number of the latitude (in degrees) used to compute the HEALPix number, starting from 1.
  #[clap(short = 'b', long, value_name = "FIELD", default_value_t = 2)]
  lat: usize,
  /// Depth of the HEALPix cumulative index (around 6 to 10, then output file will be large).
  #[arg(long, default_value_t = 8_u8)]
  depth: u8,
}

impl HealpixCumulIndex {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    if self.csv.input == PathBuf::from(r"-") {
      return Err(String::from("Command 'hcidx' not compatible with reading from stdin.").into());
    }
    // Open the file and get line iterator
    let mut reader = File::open(&self.csv.input).map(BufReader::new)?;
    let mut line = String::new();
    let mut n_bytes_read = 0;
    let mut n_bytes = reader.read_line(&mut line)?;
    // Handle comments lines
    while n_bytes > 0 && line.starts_with('#') {
      n_bytes_read += n_bytes;
      line.clear();
      n_bytes = reader.read_line(&mut line)?;
    }
    // Handle header line (if any)
    if self.csv.header && n_bytes > 0 {
      n_bytes_read += n_bytes;
      line.clear();
      n_bytes = reader.read_line(&mut line)?;
    }
    // Prepare the method computing HEALPix hash values from a CSV line
    let hpx = get_hpx(
      self.lon - 1,
      self.lat - 1,
      self.csv.delimiter,
      get(self.depth),
    );
    // Prepare output
    let fits_file = File::create(self.output)?;
    let out_fits_write = BufWriter::new(fits_file);
    // Prepare building the map
    let len = n_hash(self.depth) + 1;
    let mut map: Vec<u64> = Vec::with_capacity(len as usize);
    // Read line by line
    let mut irow = 0;
    while n_bytes > 0 {
      line.pop(); // removes the ending `'\n'
      let icell = hpx(&line);
      if icell + 1 < map.len() as u64 {
        return Err(
          format!(
            "HEALPix error at row {}: the file seems not ot be sorted!",
            irow
          )
          .into(),
        );
      }
      // Push only the starting byte of the first row having a given cell number.
      // Copy the value for all empty cells between two non-empty cells.
      for _ in map.len() as u64..=icell {
        //info!("Push row: {}; bytes: {:?}", irow, &byte_range);
        map.push(n_bytes_read as u64);
      }
      n_bytes_read += n_bytes;
      irow += 1;
      line.clear();
      n_bytes = reader.read_line(&mut line)?;
    }
    // Complete the map if necessary
    for _ in map.len() as u64..len {
      map.push(n_bytes_read as u64);
    }
    // Write the cumulative map
    let file_metadata = self.csv.input.metadata().ok();
    OwnedCIndex::new_unchecked(self.depth, map.into_boxed_slice())
      .to_fits(
        out_fits_write,
        self.csv.input.file_name().and_then(|name| name.to_str()),
        file_metadata.as_ref().map(|meta| meta.len()),
        None, // So far we do not compute the md5 of the VOTable!
        file_metadata.as_ref().and_then(|meta| meta.modified().ok()),
        Some(format!("#{}", self.lon).as_str()),
        Some(format!("#{}", self.lat).as_str()),
      )
      .map_err(|e| e.into())
  }
}
