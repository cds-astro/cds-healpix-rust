//! Deals with HEALPix cells numbers at a fixed depth.

use std::{
  error::Error,
  fs::File,
  io::{stdout, BufRead, BufReader, StdoutLock, Write},
  path::PathBuf,
};

use clap::{Args, Subcommand};

/// Perform an operation consuming an iterator of results of (longitude, latitude) tuples, in radians.
pub trait HashItOperation {
  type R;
  /// Perform an operation an an iterator of (longitude, latitude) tuples, in radians.
  fn exec<I>(self, has_it: I) -> Result<Self::R, Box<dyn Error>>
  where
    I: Iterator<Item = Result<u64, Box<dyn Error>>>;
}

#[derive(Debug, Clone, Subcommand)]
pub enum HashListInput {
  List(HashList),
  Csv(HashCsvConsumed),
}
impl HashListInput {
  /// # Params
  /// * `op`: operation consuming the read (longitude, latitude) tuples, in radians, to produce a result.
  pub fn exec<F>(self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: HashItOperation,
  {
    match self {
      Self::List(e) => e.exec_op(op),
      Self::Csv(e) => e.exec_op(op),
    }
  }
}

#[derive(Debug, Clone, Subcommand)]
pub enum HashInput {
  Value(HashVal),
  List(HashList),
  Csv(HashCsv),
}
impl HashInput {
  /// # Params
  /// * `f`: input parameters are
  ///     + hash value
  ///     + writer
  pub fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    match self {
      Self::Value(e) => e.exec(f),
      Self::List(e) => e.exec(f),
      Self::Csv(e) => e.exec(f),
    }
  }
}

/// Single hash value (aka ipix or cell number) provided in the command line
#[derive(Debug, Clone, Args)]
pub struct HashVal {
  /// Healpix hash value (aka ipix or cell number)
  #[clap(value_name = "HASH/IPIX")]
  pub hash: u64,
}

impl HashVal {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    f(self.hash, &mut lock, ' ').and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))
  }

  /*pub fn exec_op<F>(&self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: HashItOperation,
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
    F: HashItOperation,
  {
    let it = read.lines().map(listline2hash);
    op.exec(it)
  }*/
}

/// List of hash value (aka ipix or cell number); one value per line.
#[derive(Debug, Clone, Args)]
pub struct HashList {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Use DELIM as output delimiter
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
}

impl HashList {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
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
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    for line in read.lines() {
      listline2hash(line)
        .and_then(|hash| f(hash, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
    }
    Ok(())
  }

  pub fn exec_op<F>(&self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: HashItOperation,
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
    F: HashItOperation,
  {
    let it = read.lines().map(listline2hash);
    op.exec(it)
  }
}

/// Hash value (aka ipix or cell number) among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct HashCsv {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the hash value, starting from 1
  #[clap(short, long, value_name = "FIELD", default_value_t = 1)]
  pub field: usize,
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

impl HashCsv {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
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
    F: Fn(u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let ifield = self.field - 1;
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
          .and_then(|line| csvline2hash(self.delimiter, line, ifield))
          .and_then(|hash| f(hash, &mut lock, self.delimiter))
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
        csvline2hash_res(self.delimiter, line, ifield)
          .and_then(|hash| f(hash, &mut lock, self.delimiter))
          .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
      }
    }
    Ok(())
  }
}

/// Hash value (aka ipix or cell number) among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct HashCsvConsumed {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the hash value, starting from 1
  #[clap(short, long, value_name = "FIELD", default_value_t = 1)]
  pub field: usize,
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
impl HashCsvConsumed {
  /// WARNING: this method does not takes into account `parallel` and `chunk_size`.
  fn exec_op<F>(&self, op: F) -> Result<F::R, Box<dyn Error>>
  where
    F: HashItOperation,
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
    F: HashItOperation,
  {
    let ifield = self.field - 1;

    let mut lines = read.lines();
    if self.header {
      if let Some(first_line) = lines.next() {
        first_line?;
      }
    }

    op.exec(lines.map(|line| csvline2hash_res(self.delimiter, line, ifield)))
  }
}

fn listline2hash(line: std::io::Result<String>) -> Result<u64, Box<dyn Error>> {
  line
    .map_err(|e| e.into())
    .and_then(|line| line.trim().parse::<u64>().map_err(|e| e.into()))
}

fn csvline2hash_res(
  separator: char,
  line: std::io::Result<String>,
  index_col: usize,
) -> Result<u64, Box<dyn Error>> {
  line
    .map_err(|e| e.into())
    .and_then(|line| csvline2hash(separator, line, index_col))
}

fn csvline2hash(separator: char, line: String, index_col: usize) -> Result<u64, Box<dyn Error>> {
  line
    .as_str()
    .split(separator)
    .nth(index_col)
    .ok_or_else(|| {
      format!(
        "Less than {} field delimited by '{}' in row '{}'",
        index_col + 1,
        separator,
        line
      )
      .into()
    })
    .and_then(|s| {
      s.trim()
        .parse::<u64>()
        .map_err(|e| format!("Error parsing '{}': {:?}", s, e).into())
    })
}
