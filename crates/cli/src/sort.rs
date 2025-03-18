use std::{error::Error, path::PathBuf};

use clap::Args;

use hpxlib::nested::sort::{
  hpx_external_sort_csv_file, hpx_external_sort_csv_file_stdout, hpx_external_sort_csv_stdin,
  hpx_external_sort_csv_stdin_stdout, SimpleExtSortParams,
};

use crate::input::csv::Csv;

/// Sorts a CSV file by order 29 HEALPix NESTED indices, uses external sort to support huge files.
#[derive(Debug, Clone, Args)]
pub struct Sort {
  #[command(flatten)]
  csv: Csv,
  /// Field number of the longitude(in degrees) used to compute the HEALPix number, starting from 1.
  #[clap(short = 'l', long, value_name = "FIELD", default_value_t = 1)]
  lon: usize,
  /// Field number of the latitude (in degrees) used to compute the HEALPix number, starting from 1.
  #[clap(short = 'b', long, value_name = "FIELD", default_value_t = 2)]
  lat: usize,
  /// Path of the output file, else write to stdout
  #[clap(short = 'o', long = "out", value_name = "FILE")]
  output: Option<PathBuf>,
  /// Set the number of threads used [default: use all available threads]
  #[arg(long, value_name = "N")]
  parallel: Option<usize>,
  /// Do not use external sort. Faster, but use only with table holding in RAM.
  #[arg(short = 'f', long = "full-in-mem")]
  fully_in_memory: bool,
  /// Directory containing the temporary directories/files for external sort.
  #[arg(long, default_value = ".sort_tmp/")]
  tmp_dir: PathBuf,
  /// Number of rows per external sort chunk (2 chunks are simultaneously loaded in memory).
  #[arg(long, default_value_t = 50_000_usize)]
  chunk_size: usize,
  /// Depth of the computed HEALPix count map for the external sort. Should be deep enough so that
  /// the largest count map value is smaller than `chunk-size`.
  #[arg(long, default_value_t = 8_u8)]
  depth: u8,
  /// Save the computed count map in the given FITS file path.
  #[arg(long)]
  count_map_path: Option<PathBuf>,
}

impl Sort {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    self.choose_input_and_exec()
  }

  pub fn choose_input_and_exec(self) -> Result<(), Box<dyn Error>> {
    let sort_params =
      SimpleExtSortParams::new(self.tmp_dir, self.chunk_size as u32, self.parallel, true);
    if self.csv.input == PathBuf::from(r"-") {
      match self.output {
        None => hpx_external_sort_csv_stdin_stdout(
          self.lon - 1,
          self.lat - 1,
          self.csv.header,
          Some(self.csv.delimiter),
          self.depth,
          self.count_map_path,
          Some(sort_params),
        ),
        Some(path) => hpx_external_sort_csv_stdin(
          path,
          true,
          self.lon - 1,
          self.lat - 1,
          self.csv.header,
          Some(self.csv.delimiter),
          self.depth,
          self.count_map_path,
          Some(sort_params),
        ),
      }
    } else {
      match self.output {
        None => hpx_external_sort_csv_file_stdout(
          self.csv.input,
          self.lon - 1,
          self.lat - 1,
          self.csv.header,
          Some(self.csv.delimiter),
          self.depth,
          self.count_map_path,
          Some(sort_params),
        ),
        Some(path) => hpx_external_sort_csv_file(
          self.csv.input,
          path,
          true,
          self.lon - 1,
          self.lat - 1,
          self.csv.header,
          Some(self.csv.delimiter),
          self.depth,
          self.count_map_path,
          Some(sort_params),
        ),
      }
    }
  }
}
