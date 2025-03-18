//! Here, deal with `(x, y)` coordinates in the HEALPix (2D) projection plane.

use std::{
  error::Error,
  fs::File,
  io::{stdout, BufRead, BufReader, StdoutLock, Write},
  path::PathBuf,
};

use clap::{Args, Subcommand};

#[derive(Debug, Clone, Subcommand)]
pub enum CooInput {
  Value(CooVal),
  List(CooList),
  Csv(CooCsv),
}
impl CooInput {
  /// # Params
  /// * `f`: input parameters are
  ///     + x ∊ ]-8, 8[
  ///     + y ∊ [-2, 2]
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
pub struct CooVal {
  /// x ∊ [0, 8[ (x can be negative, possibly resulting in a negative longitude)
  #[clap(value_name = "X")]
  pub x: f64,
  /// y ∊ [-2, 2]
  #[clap(allow_negative_numbers = true, value_name = "Y")]
  pub y: f64,
}

impl CooVal {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(f64, f64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    let x = check_x(self.x)?;
    let y = check_y(self.y)?;
    f(x, y, &mut lock, ' ').and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))
  }
}

/// List of DELIM separated projected (x, y) coordinates; one position per line.
#[derive(Debug, Clone, Args)]
pub struct CooList {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Use DELIM delimiter
  #[clap(value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
}

impl CooList {
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
      listline2xy(self.delimiter, line)
        .and_then(|(x, y)| f(x, y, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
    }
    Ok(())
  }
}

/// Projected coordinates (x, y), among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct CooCsv {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the x coordinate, starting from 1
  #[clap(short = 'x', long, value_name = "FIELD", default_value_t = 1)]
  pub x: usize,
  /// Field number of the y coordinate, starting from 1
  #[clap(short = 'y', long, value_name = "FIELD", default_value_t = 2)]
  pub y: usize,
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

impl CooCsv {
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
    let (index_first_col, offset_to_second_col, x_first) = if self.x < self.y {
      (self.x - 1, self.y - self.x - 1, true)
    } else {
      (self.y - 1, self.x - self.y - 1, false)
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
            csvline2xy(
              self.delimiter,
              line,
              index_first_col,
              offset_to_second_col,
              x_first,
            )
          })
          .and_then(|(x, y)| f(x, y, &mut lock, self.delimiter))
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
        csvline2xy_res(
          self.delimiter,
          line,
          index_first_col,
          offset_to_second_col,
          x_first,
        )
        .and_then(|(x, y)| f(x, y, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
      }
    }
    Ok(())
  }
}

fn listline2xy(
  separator: char,
  line: std::io::Result<String>,
) -> Result<(f64, f64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    line
      .trim()
      .split_once(separator)
      .ok_or_else(|| format!("Split failed on '{}' with delimiter '{}'.", line, separator).into())
      .and_then(|(x_str, y_str)| {
        let x = x_str.parse::<f64>().map_err(|e| e.into()).and_then(check_x);
        let y = y_str.parse::<f64>().map_err(|e| e.into()).and_then(check_y);
        match (x, y) {
          (Ok(x), Ok(y)) => Ok((x, y)),
          (Err(e), _) => Err(e),
          (_, Err(e)) => Err(e),
        }
      })
  })
}

fn csvline2xy_res(
  separator: char,
  line: std::io::Result<String>,
  index_first_col: usize,
  offset_to_second_col: usize,
  x_first: bool,
) -> Result<(f64, f64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    csvline2xy(
      separator,
      line,
      index_first_col,
      offset_to_second_col,
      x_first,
    )
  })
}

fn csvline2xy(
  separator: char,
  line: String,
  index_first_col: usize,
  offset_to_second_col: usize,
  x_first: bool,
) -> Result<(f64, f64), Box<dyn Error>> {
  let mut field_it = line.as_str().split(separator);
  let field1 = field_it
    .nth(index_first_col)
    .ok_or_else(|| {
      format!(
        "Less than {} field delimited by '{}' in row '{}'",
        index_first_col + 1,
        separator,
        line
      )
      .into()
    })
    .and_then(|s| s.parse::<f64>().map_err(|e| e.into()));
  let field2 = field_it
    .nth(offset_to_second_col)
    .ok_or_else(|| {
      format!(
        "Less than {} field delimited by '{}' in row '{}'",
        index_first_col + offset_to_second_col + 1,
        separator,
        line
      )
      .into()
    })
    .and_then(|s| s.parse::<f64>().map_err(|e| e.into()));
  let (x, y) = if x_first {
    (field1, field2)
  } else {
    (field2, field1)
  };
  match (x.and_then(check_x), y.and_then(check_y)) {
    (Ok(x), Ok(y)) => Ok((x, y)),
    (Err(e), _) => Err(e),
    (_, Err(e)) => Err(e),
  }
}

fn check_x(x: f64) -> Result<f64, Box<dyn Error>> {
  if !(-8.0..8.0).contains(&x) {
    Err(String::from("x must be in [-8, 8]").into())
  } else {
    Ok(x)
  }
}

fn check_y(y: f64) -> Result<f64, Box<dyn Error>> {
  if !(-2.0..=2.0).contains(&y) {
    Err(String::from("y must be in [-2, 2]").into())
  } else {
    Ok(y)
  }
}
