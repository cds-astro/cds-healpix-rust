use std::path::PathBuf;

use clap::Args;

#[derive(Debug, Clone, Args)]
pub struct Csv {
  /// Path of the input file ('-' for stdin).
  #[clap(value_name = "FILE", default_value = "-")]
  pub input: PathBuf,
  /// Use DELIM delimiter.
  #[clap(short, long, value_name = "DELIM", default_value_t = '\t')]
  pub delimiter: char,
  /// Treat the first line as field headers.
  #[clap(long)]
  pub header: bool,
}
