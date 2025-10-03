use std::{
  error::Error,
  fs::{self, File},
  path::PathBuf,
};

use clap::{Args, Subcommand};
use memmap2::MmapOptions;

use hpxlib::nested::sort::cindex::FitsMMappedCIndex;
use hpxlib::nested::{
  bmoc::BMOC,
  sort::cindex::{FITSCIndex, HCIndex},
};

/// Type of the query (HEALPix cell or BMOC).
#[derive(Debug, Subcommand)]
pub enum QueryType {
  /// Query a given HEALPix cell
  Cell {
    /// Depth of the queried HEALPix cell.
    depth: u8,
    /// Index of the queried HEALPix cell.
    ipix: u64,
  },
  /// Query a given BMOC
  BMOC {
    /// Path of the BMOC
    path: PathBuf,
  },
}
impl QueryType {
  // FitsMMappedCIndex
  //pub fn exec<H: HCIndex<V = u64>>(self, csv_name: String, hci: H) -> Result<(), Box<dyn Error>> {
  pub fn exec<'a, H, T>(self, fits_hci: &'a T) -> Result<(), Box<dyn Error>>
  where
    H: HCIndex<V = u64>,
    T: FitsMMappedCIndex<'a, HCIndexType = H> + 'a,
  {
    let csv_name = fits_hci
      .get_indexed_file_name()
      .expect("No file name found in the FITS HCI file.");
    let expected_csv_len = fits_hci
      .get_indexed_file_len()
      .expect("No file length found in the FITS HCI file.");
    check_file_exists_and_check_file_len(csv_name, expected_csv_len)?;

    let hci = fits_hci.get_hcindex();
    let hci_depth = hci.depth();
    match self {
      Self::Cell { depth, ipix } => {
        if hci_depth < depth {
          return Err(
            format!(
              "Query depth ({}) larger than index depth ({}) not supported yet.",
              depth, hci_depth,
            )
            .into(),
          );
        }
        // Performs the query
        let data_start = hci.get(0) as usize;
        let bytes_range = hci.get_cell(depth, ipix);
        let file = File::open(csv_name)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let mut stdout = std::io::stdout();
        std::io::copy(&mut &mmap[0..data_start], &mut stdout)
          .and_then(|_| {
            std::io::copy(
              &mut &mmap[bytes_range.start as usize..bytes_range.end as usize],
              &mut stdout,
            )
          })
          .map(|_| ())
          .map_err(|e| e.into())
      }
      Self::BMOC { path } => {
        let bmoc = BMOC::from_fits_file(path)?;
        let bmoc_depth = bmoc.get_depth_max();
        // Check bmoc depth
        if hci_depth < bmoc_depth {
          return Err(
            format!(
              "BMOC depth ({}) larger than index depth ({}) not supported yet.",
              bmoc_depth, hci_depth,
            )
            .into(),
          );
        }
        // Performs the query
        let data_start = hci.get(0) as usize;
        let file = File::open(csv_name)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let mut stdout = std::io::stdout();
        // Copy header
        std::io::copy(&mut &mmap[0..data_start], &mut stdout)?;
        for range in bmoc.to_ranges() {
          let bytes_range = hci.get_range(bmoc_depth, range);
          std::io::copy(
            &mut &mmap[bytes_range.start as usize..bytes_range.end as usize],
            &mut stdout,
          )?;
        }
        Ok(())
      }
    }
  }
}

/// Query an HEALPix sorted and indexed CSV file (see the 'hcidx' command).
#[derive(Debug, Args)]
pub struct QueryHCIndex {
  /// Path of the input HEALPix Cumulative Index FITS file (it contains the name of the CSV to query).
  #[clap(value_name = "NAME.hci.fits")]
  hcindex: PathBuf,
  #[clap(subcommand)]
  query: QueryType,
}

impl QueryHCIndex {
  pub fn exec(self) -> Result<(), Box<dyn Error>> {
    match FITSCIndex::from_fits_file(self.hcindex)? {
      FITSCIndex::ImplicitU64(fits_hci) => self.query.exec(&fits_hci),
      FITSCIndex::ExplicitU32U64(fits_hci) => self.query.exec(&fits_hci),
      FITSCIndex::ExplicitU64U64(fits_hci) => self.query.exec(&fits_hci),
      _ => Err(
        String::from("Wrong data type in the FITS Healpix Cumulative Index type. Expected: u64.")
          .into(),
      ),
    }
  }
}

fn check_file_exists_and_check_file_len(
  csv_name: &String,
  expected_csv_len: u64,
) -> Result<(), Box<dyn Error>> {
  // Check if file exists
  fs::exists(csv_name)
    .map_err(|e| e.into())
    .and_then(|exists| {
      if exists {
        Ok::<(), Box<dyn Error>>(())
      } else {
        Err(format!("File `{}` not found in the current directory.", csv_name).into())
      }
    })?;
  // Check file len
  let actual_csv_len = fs::metadata(csv_name).map(|metadata| metadata.len())?;
  if actual_csv_len != expected_csv_len {
    Err(
      format!(
        "Local CSV file `{}` len does not match index info. Expected: {}. Actual: {}.",
        csv_name, expected_csv_len, actual_csv_len
      )
      .into(),
    )
  } else {
    Ok(())
  }
}
