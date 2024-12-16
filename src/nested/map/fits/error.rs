use std::{io, num::ParseIntError};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum FitsError {
  /// IO error
  #[error("I/O error.")]
  Io(#[from] io::Error),
  #[error("I/O error. Path: {path:}. Error: {err:?}.")]
  IoWithPath { path: String, err: io::Error },
  #[error("Wrong FITS keyword. Expected: {expected:?}. Actual: {actual:?}).")]
  UnexpectedKeyword { expected: String, actual: String },
  #[error("Value indicator not found in keyword record '{keyword_record:?}'.")]
  ValueIndicatorNotFound { keyword_record: String },
  #[error("Wrong value for keyword '{keyword:}'. Expected: '{expected:}'. Actual: '{actual:}'.")]
  UnexpectedValue {
    keyword: String,
    expected: String,
    actual: String,
  },
  #[error("Unsigned int value not found in keyword record '{keyword_record:}'.")]
  UintValueNotFound { keyword_record: String },
  #[error("String value no found in keyword record '{keyword_record:}'.")]
  StringValueNotFound { keyword_record: String },
  #[error("Parse {context:}. Error: {err:?}")]
  WrongUintValue { context: String, err: ParseIntError },
  #[error("FITS not valid. Multiple Keyword '{keyword:}'.")]
  MultipleKeyword { keyword: String },
  #[error("Missing keyword '{keyword:}'.")]
  MissingKeyword { keyword: String },
  #[error("Incompatible keyword values for {keyword1:} and {keyword2:}.")]
  IncompatibleKeywordContent { keyword1: String, keyword2: String },
  #[error("More data than the expected!")]
  RemainingData,
  #[error("Less data than expected!")]
  PrematureEndOfData,
  #[error("Unexpected number of data written!")]
  UnexpectedWrittenSize,
  #[error("unexpected depth. Max expected: {depth_max:}. Actual: {depth:}")]
  UnexpectedDepth { depth: u8, depth_max: u8 },
  #[error("FITS not valid: '{msg:}'.")]
  Custom { msg: String },
}

impl FitsError {
  pub fn new_custom(msg: String) -> Self {
    Self::Custom { msg }
  }
}
