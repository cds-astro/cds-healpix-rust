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
    cindex::{HCIndex, HCIndexShape, OwnedCIndex, OwnedCIndexExplicit},
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
  #[arg(long, default_value_t = 8_u8)] // '-d' used for delimiter in CSV
  depth: u8,
  #[arg(short = 'e', long)]
  /// Use in-memory `explicit` representation instead of `implicit`.
  explicit: bool,
  #[arg(short = 'r', long = "ratio")]
  /// Limit on the ratio of the implicit over the explicit byte sizes for the FITS serialisation.
  /// Above the limit, the explicit representation is chosen.
  /// Below the limit, the implicit representation is chosen.
  /// If unset, the in-memory representation is chosen.
  implicit_over_explicit_ratio: Option<f64>,
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

    if self.explicit {
      let len = n_hash(self.depth) + 1;
      let mut entries: Vec<(u64, u64)> = Vec::with_capacity(len.min(1_000_000) as usize);
      // Read line by line
      let mut prev_icell = 0;
      let mut irow = 0;
      while n_bytes > 0 {
        line.pop(); // removes the ending '\n'
        let icell = hpx(&line);
        if icell + 1 < prev_icell {
          return Err(
            format!(
              "HEALPix error at row {}: the file seems not to be sorted!",
              irow
            )
            .into(),
          );
        }
        // Push only the starting byte of the first row having a given cell number.
        if icell > prev_icell || (icell == 0 && entries.is_empty()) {
          entries.push((icell, n_bytes_read as u64));
        }
        n_bytes_read += n_bytes;
        irow += 1;
        line.clear();
        n_bytes = reader.read_line(&mut line)?;
        prev_icell = icell;
      }
      entries.push((prev_icell + 1, n_bytes_read as u64));
      // Write the cumulative map
      let explicit_index = OwnedCIndexExplicit::new_unchecked(self.depth, entries);
      self.write_index(explicit_index)
    } else {
      // Prepare building the map
      let len = n_hash(self.depth) + 1;
      let mut map: Vec<u64> = Vec::with_capacity(len as usize);
      // Read line by line
      let mut irow = 0;
      while n_bytes > 0 {
        line.pop(); // removes the ending '\n'
        let icell = hpx(&line);
        if icell + 1 < map.len() as u64 {
          return Err(
            format!(
              "HEALPix error at row {}: the file seems not to be sorted!",
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
      let implicit_index = OwnedCIndex::new_unchecked(self.depth, map.into_boxed_slice());
      self.write_index(implicit_index)
    }
  }

  fn write_index<H: HCIndex>(self, cindex: H) -> Result<(), Box<dyn Error>> {
    // Prepare output
    let fits_file = File::create(self.output)?;
    let out_fits_write = BufWriter::new(fits_file);
    let file_metadata = self.csv.input.metadata().ok();
    let best_repr = self
      .implicit_over_explicit_ratio
      .map(|ratio| cindex.best_representation(ratio))
      .unwrap_or_else(|| {
        if self.explicit {
          HCIndexShape::Explicit
        } else {
          HCIndexShape::Implicit
        }
      });
    match best_repr {
      HCIndexShape::Implicit => cindex.to_fits_implicit(
        out_fits_write,
        self.csv.input.file_name().and_then(|name| name.to_str()),
        file_metadata.as_ref().map(|meta| meta.len()),
        None, // So far we do not compute the md5 of the VOTable!
        file_metadata.as_ref().and_then(|meta| meta.modified().ok()),
        Some(format!("#{}", self.lon).as_str()),
        Some(format!("#{}", self.lat).as_str()),
      ),
      HCIndexShape::Explicit => cindex.to_fits_explicit(
        out_fits_write,
        self.csv.input.file_name().and_then(|name| name.to_str()),
        file_metadata.as_ref().map(|meta| meta.len()),
        None, // So far we do not compute the md5 of the VOTable!
        file_metadata.as_ref().and_then(|meta| meta.modified().ok()),
        Some(format!("#{}", self.lon).as_str()),
        Some(format!("#{}", self.lat).as_str()),
      ),
    }
    .map_err(|e| e.into())
  }
}
