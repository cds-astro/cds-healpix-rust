use std::{
  error::Error,
  fs::File,
  io::{stdout, BufRead, BufReader, StdoutLock, Write},
  path::PathBuf,
};

use clap::{Args, Subcommand};

use hpxlib::{n_hash, DEPTH_MAX};

#[derive(Debug, Clone, Subcommand)]
pub enum DepthHashInput {
  Value(DepthHashVal),
  List(DepthHashList),
  Csv(DepthHashCsv),
}
impl DepthHashInput {
  /// # Params
  /// * `f`: input parameters are
  ///     + depth ∊ [0, 29]
  ///     + hash ∊ [0, n_hash(depth)]
  ///     + writer
  ///     + seprator
  pub fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    match self {
      Self::Value(e) => e.exec(f),
      Self::List(e) => e.exec(f),
      Self::Csv(e) => e.exec(f),
    }
  }
}

/// Single depth/hash value provided in the command line
#[derive(Debug, Clone, Args)]
pub struct DepthHashVal {
  /// Healpix depth (aka order) ∊ [0, 29]
  #[clap(value_name = "DEPTH")]
  pub depth: u8,
  /// Healpix hash (aka ipix or cell number)
  #[clap(allow_negative_numbers = true, value_name = "HASH")]
  pub hash: u64,
}

impl DepthHashVal {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    let depth = check_depth(self.depth)?;
    let hash = check_hash(self.hash, n_hash(depth))?;
    f(depth, hash, &mut lock, ' ').and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))
  }
}

/// List of DELIM separated depth and hash value; one tuple per line.
#[derive(Debug, Clone, Args)]
pub struct DepthHashList {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Use DELIM delimiter
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
}

impl DepthHashList {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
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
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let mut lock = stdout().lock();
    for line in read.lines() {
      listline2dh(self.delimiter, line)
        .and_then(|(depth, hash)| f(depth, hash, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
    }
    Ok(())
  }
}

/// Depth (aka order) and hash value (aka ipix), among the fields of a DELIM separated file.
#[derive(Debug, Clone, Args)]
pub struct DepthHashCsv {
  /// Path of the input file ('-' for stdin)
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Field number of the depth, starting from 1
  #[clap(short = '1', long = "depth", value_name = "FIELD", default_value_t = 1)]
  pub depth_fieldnum: usize,
  /// Field number of the hash value, starting from 1
  #[clap(short = '2', long = "hash", value_name = "FIELD", default_value_t = 2)]
  pub hash_fieldnum: usize,
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

impl DepthHashCsv {
  fn exec<F>(self, f: F) -> Result<(), Box<dyn Error>>
  where
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
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
    F: Fn(u8, u64, &mut StdoutLock<'static>, char) -> Result<(), Box<dyn Error>>,
  {
    let (index_first_col, offset_to_second_col, id, ih) =
      if self.depth_fieldnum < self.hash_fieldnum {
        (
          self.depth_fieldnum - 1,
          self.hash_fieldnum - self.depth_fieldnum - 1,
          0,
          1,
        )
      } else {
        (
          self.hash_fieldnum - 1,
          self.depth_fieldnum - self.hash_fieldnum - 1,
          1,
          0,
        )
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
            csvline2dh(
              self.delimiter,
              line,
              index_first_col,
              offset_to_second_col,
              id,
              ih,
            )
          })
          .and_then(|(d, h)| f(d, h, &mut lock, self.delimiter))
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
        csvline2dh_res(
          self.delimiter,
          line,
          index_first_col,
          offset_to_second_col,
          id,
          ih,
        )
        .and_then(|(d, h)| f(d, h, &mut lock, self.delimiter))
        .and_then(|()| writeln!(&mut lock).map_err(|e| e.into()))?;
      }
    }
    Ok(())
  }
}

fn listline2dh(
  separator: char,
  line: std::io::Result<String>,
) -> Result<(u8, u64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    line
      .trim()
      .split_once(separator)
      .ok_or_else(|| format!("Split failed on '{}' with delimiter '{}'.", line, separator).into())
      .and_then(|(d_str, h_str)| dhstr2dh(d_str, h_str))
  })
}

fn csvline2dh_res(
  separator: char,
  line: std::io::Result<String>,
  index_first_col: usize,
  offset_to_second_col: usize,
  id: usize,
  ih: usize,
) -> Result<(u8, u64), Box<dyn Error>> {
  line.map_err(|e| e.into()).and_then(|line| {
    csvline2dh(
      separator,
      line,
      index_first_col,
      offset_to_second_col,
      id,
      ih,
    )
  })
}

fn csvline2dh(
  separator: char,
  line: String,
  index_first_col: usize,
  offset_to_second_col: usize,
  id: usize,
  ih: usize,
) -> Result<(u8, u64), Box<dyn Error>> {
  let mut field_it = line.as_str().split(separator);
  let field1 = field_it.nth(index_first_col).ok_or_else(|| {
    format!(
      "Less than {} field delimited by '{}' in row '{}'",
      index_first_col + 1,
      separator,
      line
    )
  })?;
  let field2 = field_it.nth(offset_to_second_col).ok_or_else(|| {
    format!(
      "Less than {} field delimited by '{}' in row '{}'",
      index_first_col + offset_to_second_col + 1,
      separator,
      line
    )
  })?;
  let elems = [field1, field2];
  dhstr2dh(elems[id], elems[ih])
}

fn dhstr2dh(d_str: &str, h_str: &str) -> Result<(u8, u64), Box<dyn Error>> {
  d_str
    .parse::<u8>()
    .map_err(|e| e.into())
    .and_then(check_depth)
    .and_then(|d| {
      h_str
        .parse::<u64>()
        .map_err(|e| e.into())
        .and_then(|h| check_hash(h, n_hash(d)))
        .map(|h| (d, h))
    })
}

fn check_depth(depth: u8) -> Result<u8, Box<dyn Error>> {
  if depth > DEPTH_MAX {
    Err(format!("Wrong depth. Expected: <{}. Actual: {}.", DEPTH_MAX, depth).into())
  } else {
    Ok(depth)
  }
}

fn check_hash(hash: u64, n_hash: u64) -> Result<u64, Box<dyn Error>> {
  if hash >= n_hash {
    Err(format!("Wrong hash value. Expected: <{}. Actual: {}", hash, n_hash).into())
  } else {
    Ok(hash)
  }
}
