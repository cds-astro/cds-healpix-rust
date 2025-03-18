//! Here deals with sky positions, more precisely with Equatorial coordiantes `(lon, lat)`.

use std::{
  error::Error,
  f64::consts::FRAC_PI_2,
  fs::File,
  io::{stdout, BufRead, BufReader, StdoutLock, Write},
  path::PathBuf,
};

use clap::{Args, Subcommand};

use hpxlib::TWICE_PI;

/// Perform an operation consuming an iterator of results of (longitude, latitude) tuples, in radians.
pub trait PosItOperation {
  type R;
  /// Perform an operation an an iterator of (longitude, latitude) tuples, in radians.
  fn exec<I>(self, pos_it: I) -> Result<Self::R, Box<dyn Error>>
  where
    I: Iterator<Item = Result<(f64, f64), Box<dyn Error>>>;
}

#[derive(Debug, Clone, Subcommand)]
pub enum PosListInput {
  List(PosList),
  Csv(PosCsvConsumed),
}
impl PosListInput {
  /// # Params
  /// * `op`: operation consuming the read (longitude, latitude) tuples, in radians, to produce a result.
  pub fn exec<F>(self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: PosItOperation,
  {
    match self {
      Self::List(e) => e.exec_op(op),
      Self::Csv(e) => e.exec_op(op),
    }
  }
}
#[derive(Debug, Clone, Subcommand)]
pub enum PosInput {
  Value(PosVal),
  List(PosList),
  Csv(PosCsv),
}
impl PosInput {
  /// # Params
  /// * `f`: input parameters are
  ///     + longitude in decimal degrees
  ///     + latitude in decimal degrees
  ///     + writer
  ///     + seprator
  pub fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    match self {
      Self::Value(e) => e.exec(f),
      Self::List(e) => e.exec(f),
      Self::Csv(e) => e.exec(f),
    }
  }
}

/// Single position provided in the command line
#[derive(Debug, Clone, Args)]
pub struct PosVal {
  /// Longitude, in degrees
  #[clap(value_name = "LON_DEG")]
  pub lon: f64,
  /// Latitude, in degrees
  #[clap(allow_negative_numbers = true, value_name = "LAT_DEG")]
  pub lat: f64,
}

impl PosVal {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    let lon_rad = lon_deg2rad(self.lon)?;
    let lat_rad = lat_deg2rad(self.lat)?;
    f(lon_rad, lat_rad, &mut lock, ' ').and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))
  }
}

/// List of DELIM separated longitude and latitudes, in decimal degrees; one position per line.
#[derive(Debug, Clone, Args)]
pub struct PosList {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Use DELIM delimiter
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
}

impl PosList {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    if self.input == PathBuf::from(r"-") {
      let stdin = std::io::stdin();
      self.exec_on_reader(stdin.lock(), f)
    } else {
      let file = File::open(&self.input)?;
      let reader = BufReader::new(file);
      self.exec_on_reader(reader, f)
    }
  }

  fn exec_on_reader<R, F>(&self, read: R, f: F) -> Result<(), Box<dyn Error>>
  where
    R: BufRead,
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    for line in read.lines() {
      listline2lonlat(self.delimiter, line)
        .and_then(|(lon_rad, lat_rad)| f(lon_rad, lat_rad, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
    }
    Ok(())
  }

  pub fn exec_op<F>(&self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: PosItOperation,
  {
    if self.input == PathBuf::from(r"-") {
      let stdin = std::io::stdin();
      self.exec_op_on_reader(stdin.lock(), op)
    } else {
      let file = File::open(&self.input)?;
      let reader = BufReader::new(file);
      self.exec_op_on_reader(reader, op)
    }
  }

  fn exec_op_on_reader<R, F>(&self, read: R, op: F) -> Result<F::R, Box<dyn Error>>
  where
    R: BufRead,
    F: PosItOperation,
  {
    let it = read
      .lines()
      .map(|line| listline2lonlat(self.delimiter, line));
    op.exec(it)
  }
}

/// Longitude and latitude, in decimal degrees, among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct PosCsv {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the longitude, starting from 1
  #[clap(short = 'l', long, value_name = "FIELD", default_value_t = 1)]
  pub lon: usize,
  /// Field number of the latitude, starting from 1
  #[clap(short = 'b', long, value_name = "FIELD", default_value_t = 2)]
  pub lat: usize,
  /// Use DELIM delimiter
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
  /// Treat the first line as field headers, print in output
  #[clap(long)]
  pub header: bool,
  /// Rewrite the full file, adding the computed column(s)
  #[clap(long)]
  pub paste: bool,
  /// Add the given header when options --paste and --header are used
  #[clap(long, value_name = "ADD", default_value = "ADDITIONAL_HEADER")]
  pub paste_header: String,
}

impl PosCsv {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    if self.input == PathBuf::from(r"-") {
      let stdin = std::io::stdin();
      self.exec_on_reader(stdin.lock(), f)
    } else {
      let file = File::open(&self.input)?;
      let reader = BufReader::new(file);
      self.exec_on_reader(reader, f)
    }
  }

  fn exec_on_reader<R, F>(&self, read: R, f: F) -> Result<(), Box<dyn Error>>
  where
    R: BufRead,
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let (index_first_col, offset_to_second_col, longitude_first) = if self.lon < self.lat {
      (self.lon - 1, self.lat - self.lon - 1, true)
    } else {
      (self.lat - 1, self.lon - self.lat - 1, false)
    };

    let mut lock = stdout().lock();
    let mut lines = read.lines();
    if self.paste {
      if self.header {
        if let Some(first_line) = lines.next() {
          first_line.and_then(|line| {
            writeln!(
              &mut lock,
              "{}{}{}",
              &line, self.delimiter, &self.paste_header
            )
          })?;
        }
      }
      for line in lines {
        line
          .and_then(|line| write!(&mut lock, "{}{}", &line, self.delimiter).map(|()| line))
          .map_err(|e| e.into())
          .and_then(|line| {
            csvline2lonlat(
              self.delimiter,
              line,
              index_first_col,
              offset_to_second_col,
              longitude_first,
            )
          })
          .and_then(|(lon_rad, lat_rad)| f(lon_rad, lat_rad, &mut lock, self.delimiter))
          .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
      }
    } else {
      if self.header {
        // Should we write the header in paster_header?
        if let Some(first_line) = lines.next() {
          first_line?;
        }
      }
      for line in lines {
        csvline2lonlat_res(
          self.delimiter,
          line,
          index_first_col,
          offset_to_second_col,
          longitude_first,
        )
        .and_then(|(lon_rad, lat_rad)| f(lon_rad, lat_rad, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
      }
    }
    Ok(())
  }
}

/// Longitude and latitude, in decimal degrees, among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct PosCsvConsumed {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the longitude, starting from 1
  #[clap(short = 'l', long, value_name = "FIELD", default_value_t = 1)]
  pub lon: usize,
  /// Field number of the latitude, starting from 1
  #[clap(short = 'b', long, value_name = "FIELD", default_value_t = 2)]
  pub lat: usize,
  /// Use DELIM delimiter
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
  /// Treat the first line as field headers, print in output
  #[clap(long)]
  pub header: bool,
  /// Set the number of threads used [default: use all available threads]
  #[arg(long, value_name = "N")]
  pub parallel: Option<usize>,
  /// Number of rows to be processed in parallel (2 chunks are simultaneously loaded in memory)
  #[arg(long, default_value_t = 50_000_usize)]
  pub chunk_size: usize,
}

impl PosCsvConsumed {
  /// WARNING: this method does not takes into account `parallel` and `chunk_size`.
  fn exec_op<F>(&self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: PosItOperation,
  {
    if self.input == PathBuf::from(r"-") {
      let stdin = std::io::stdin();
      self.exec_op_on_reader(stdin.lock(), op)
    } else {
      let file = File::open(&self.input)?;
      let reader = BufReader::new(file);
      self.exec_op_on_reader(reader, op)
    }
  }

  /// WARNING: this method does not takes into account `parallel` and `chunk_size`.
  fn exec_op_on_reader<R, F>(&self, read: R, op: F) -> Result<F::R, Box<dyn Error>>
  where
    R: BufRead,
    F: PosItOperation,
  {
    let (index_first_col, offset_to_second_col, longitude_first) = if self.lon < self.lat {
      (self.lon - 1, self.lat - self.lon - 1, true)
    } else {
      (self.lat - 1, self.lon - self.lat - 1, false)
    };

    let mut lines = read.lines();
    if self.header {
      if let Some(first_line) = lines.next() {
        first_line?;
      }
    }

    op.exec(lines.map(|line| {
      csvline2lonlat_res(
        self.delimiter,
        line,
        index_first_col,
        offset_to_second_col,
        longitude_first,
      )
    }))
  }
}

fn listline2lonlat(
  separator: char,
  line: std::io::Result<String>,
) -> Result<(f64, f64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    line
      .trim()
      .split_once(separator)
      .ok_or_else(|| format!("Split failed on '{}' with delimiter '{}'.", line, separator).into())
      .and_then(|(lon_deg_str, lat_deg_str)| {
        let lon_rad = lon_deg_str
          .parse::<f64>()
          .map_err(|e| format!("Error parsing: '{}': {:?}", lon_deg_str, e).into())
          .and_then(lon_deg2rad);
        let lat_rad = lat_deg_str
          .parse::<f64>()
          .map_err(|e| format!("Error parsing: '{}': {:?}", lat_deg_str, e).into())
          .and_then(lat_deg2rad);
        match (lon_rad, lat_rad) {
          (Ok(lon), Ok(lat)) => Ok((lon, lat)),
          (Err(e), _) => Err(e),
          (_, Err(e)) => Err(e),
        }
      })
  })
}

fn csvline2lonlat_res(
  separator: char,
  line: std::io::Result<String>,
  index_first_col: usize,
  offset_to_second_col: usize,
  longitude_first: bool,
) -> Result<(f64, f64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    csvline2lonlat(
      separator,
      line,
      index_first_col,
      offset_to_second_col,
      longitude_first,
    )
  })
}

fn csvline2lonlat(
  separator: char,
  line: String,
  index_first_col: usize,
  offset_to_second_col: usize,
  longitude_first: bool,
) -> Result<(f64, f64), Box<dyn Error>> {
  let mut field_it = line.as_str().split(separator);
  let field1 = field_it
    .nth(index_first_col)
    .ok_or_else(|| {
      format!(
        "Pb retrieving 1st field: less than {} field delimited by '{}' in row '{}'",
        index_first_col + 1,
        separator,
        line
      )
      .into()
    })
    .and_then(|s| {
      s.parse::<f64>()
        .map_err(|e| format!("Error parsing: '{}': {:?}", s, e).into())
    });
  let field2 = field_it
    .nth(offset_to_second_col)
    .ok_or_else(|| {
      format!(
        "Pb retrieving 2nd field: less than {} field delimited by '{}' in row '{}'",
        index_first_col + offset_to_second_col + 1,
        separator,
        line
      )
      .into()
    })
    .and_then(|s| {
      s.parse::<f64>()
        .map_err(|e| format!("Error parsing: '{}': {:?}", s, e).into())
    });
  let (lon, lat) = if longitude_first {
    (field1, field2)
  } else {
    (field2, field1)
  };
  match (lon.and_then(lon_deg2rad), lat.and_then(lat_deg2rad)) {
    (Ok(lon), Ok(lat)) => Ok((lon, lat)),
    (Err(e), _) => Err(e),
    (_, Err(e)) => Err(e),
  }
}

fn lon_deg2rad(lon_deg: f64) -> Result<f64, Box<dyn Error>> {
  let mut lon = lon_deg.to_radians();
  if lon == TWICE_PI {
    lon = 0.0;
  }
  if !(0.0..TWICE_PI).contains(&lon) {
    Err(String::from("Longitude must be in [0, 2pi[").into())
  } else {
    Ok(lon)
  }
}

fn lat_deg2rad(lat_deg: f64) -> Result<f64, Box<dyn Error>> {
  let lat = lat_deg.to_radians();
  if !(-FRAC_PI_2..=FRAC_PI_2).contains(&lat) {
    Err(String::from("Latitude must be in [-pi/2, pi/2]").into())
  } else {
    Ok(lat)
  }
}
