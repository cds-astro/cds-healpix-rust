//! This module contains the MOC specific FITS keywords

use std::{fmt, slice::ChunksMut, str};

use log::warn;

use super::{
  error::FitsError,
  read::{get_keyword, get_str_val_no_quote, parse_uint_val},
  write::write_keyword_record,
};

pub trait FitsCard: Sized {
  const KEYWORD: &'static [u8; 8];

  fn keyword_str() -> &'static str {
    unsafe { str::from_utf8_unchecked(Self::KEYWORD) }
  }

  fn keyword_string() -> String {
    unsafe { String::from_utf8_unchecked(Self::KEYWORD.to_vec()) }
  }

  fn parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    Self::specific_parse_value(keyword_record)
  }

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError>;

  fn write_keyword_record(&self, keyword_record: &mut [u8]) -> Result<(), FitsError> {
    write_keyword_record(keyword_record, Self::KEYWORD, &self.to_fits_value());
    Ok(())
  }

  /// Must be in quotes `'val'` is value type is string
  fn to_fits_value(&self) -> String;

  /// Generate an error in case the parsed value does not match a pre-define list of possible values
  /// To be called in `specific_parse_value`.
  /// Essentially, it converts &str in String (because once the error is raised, the str in the
  /// read buffer are out-of-scope.
  fn predefine_val_err(parsed_value: &[u8], expected_values: &[&[u8]]) -> FitsError {
    FitsError::UnexpectedValue {
      keyword: Self::keyword_string(),
      expected: format!(
        "{:?}",
        expected_values
          .iter()
          .map(|v| unsafe { String::from_utf8_unchecked(v.to_vec()) })
          .collect::<Vec<String>>()
      ),
      actual: String::from_utf8_lossy(parsed_value).to_string(),
    }
  }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Ordering {
  Nested,
  Ring,
}
impl FitsCard for Ordering {
  const KEYWORD: &'static [u8; 8] = b"ORDERING";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"NESTED" => Ok(Ordering::Nested),
      b"RING" => Ok(Ordering::Ring),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"NESTED", b"RING"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      Ordering::Nested => "'NESTED'",
      Ordering::Ring => "'RING'",
    })
  }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum CoordSys {
  Cel,
  Gal,
}
impl FitsCard for CoordSys {
  const KEYWORD: &'static [u8; 8] = b"COORDSYS";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"C" => Ok(CoordSys::Cel),
      b"G" => Ok(CoordSys::Gal),
      b"CEL" => {
        warn!("COORDSYS value should be 'C', not 'CEL'");
        Ok(CoordSys::Cel)
      }
      b"GAL" => {
        warn!("COORDSYS value should be 'G', not 'CAL'");
        Ok(CoordSys::Gal)
      }
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"C", b"G"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      CoordSys::Cel => "'C'",
      CoordSys::Gal => "'G'",
    })
  }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum PixType {
  Healpix,
}
impl FitsCard for PixType {
  const KEYWORD: &'static [u8; 8] = b"PIXTYPE ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"HEALPIX" => Ok(PixType::Healpix),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"TCB"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from("'HEALPIX'")
  }
}

#[derive(Debug)]
pub struct Order(u8);
impl Order {
  pub fn get(&self) -> u8 {
    self.0
  }
}
impl FitsCard for Order {
  const KEYWORD: &'static [u8; 8] = b"ORDER   ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u8>(keyword_record).map(Self)
  }
  fn to_fits_value(&self) -> String {
    format!("{}", &self.0)
  }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum TForm1 {
  B(Option<u32>), // u8
  I(Option<u32>), // u16
  J(Option<u32>), // u32
  K(Option<u32>), // u64
  E(Option<u32>), // f32
  D(Option<u32>), // f64
}
impl TForm1 {
  pub fn type_n_bytes(&self) -> usize {
    match self {
      TForm1::B(_) => size_of::<u8>(),
      TForm1::I(_) => size_of::<u16>(),
      TForm1::J(_) => size_of::<u32>(),
      TForm1::K(_) => size_of::<u64>(),
      TForm1::E(_) => size_of::<f32>(),
      TForm1::D(_) => size_of::<f64>(),
    }
  }
  pub fn n_bytes(&self) -> usize {
    match self {
      TForm1::B(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<u8>(),
      TForm1::I(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<u16>(),
      TForm1::J(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<u32>(),
      TForm1::K(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<u64>(),
      TForm1::E(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<f32>(),
      TForm1::D(opt_n) => opt_n.unwrap_or(1) as usize * size_of::<f64>(),
    }
  }
  pub fn n_pack(&self) -> u32 {
    match self {
      TForm1::B(opt_n) => opt_n.unwrap_or(1),
      TForm1::I(opt_n) => opt_n.unwrap_or(1),
      TForm1::J(opt_n) => opt_n.unwrap_or(1),
      TForm1::K(opt_n) => opt_n.unwrap_or(1),
      TForm1::E(opt_n) => opt_n.unwrap_or(1),
      TForm1::D(opt_n) => opt_n.unwrap_or(1),
    }
  }
}
impl fmt::Display for TForm1 {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "{} = {}", Self::keyword_str(), self.to_fits_value())
  }
}
impl FitsCard for TForm1 {
  const KEYWORD: &'static [u8; 8] = b"TFORM1  ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    let val = get_str_val_no_quote(keyword_record)?;
    if val.len() == 1 {
      match val {
        b"B" => Ok(TForm1::B(None)),
        b"I" => Ok(TForm1::I(None)),
        b"J" => Ok(TForm1::J(None)),
        b"K" => Ok(TForm1::K(None)),
        b"E" => Ok(TForm1::E(None)),
        b"D" => Ok(TForm1::D(None)),
        parsed_val => Err(Self::predefine_val_err(
          parsed_val,
          &[b"nB", b"nI", b"nJ", b"nK", b"nE", b"nD"],
        )),
      }
    } else {
      let (n, k) = val.split_at(val.len() - 1);
      let n_str = unsafe { str::from_utf8_unchecked(n) };
      let n = n_str
        .parse::<u32>()
        .map_err(|err| FitsError::WrongUintValue {
          context: n_str.to_string(),
          err,
        })?;
      match k {
        b"B" => Ok(TForm1::B(Some(n))),
        b"I" => Ok(TForm1::I(Some(n))),
        b"J" => Ok(TForm1::J(Some(n))),
        b"K" => Ok(TForm1::K(Some(n))),
        b"E" => Ok(TForm1::E(Some(n))),
        b"D" => Ok(TForm1::D(Some(n))),
        parsed_val => Err(Self::predefine_val_err(
          parsed_val,
          &[b"nB", b"nI", b"nJ", b"nK", b"nE", b"nD"],
        )),
      }
    }
  }

  fn to_fits_value(&self) -> String {
    match self {
      TForm1::B(opt_n) => format!(
        "'{}B'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
      TForm1::I(opt_n) => format!(
        "'{}I'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
      TForm1::J(opt_n) => format!(
        "'{}J'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
      TForm1::K(opt_n) => format!(
        "'{}K'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
      TForm1::E(opt_n) => format!(
        "'{}E'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
      TForm1::D(opt_n) => format!(
        "'{}D'",
        opt_n.map(|n| n.to_string()).unwrap_or(String::from(""))
      ),
    }
  }
}

#[derive(Debug)]
pub struct TType1(String);
impl TType1 {
  pub fn get(&self) -> &str {
    self.0.as_str()
  }
}
impl FitsCard for TType1 {
  const KEYWORD: &'static [u8; 8] = b"TTYPE1  ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    get_str_val_no_quote(keyword_record)
      .map(|s| String::from_utf8_lossy(s).to_string())
      .map(Self)
  }
  fn to_fits_value(&self) -> String {
    format!("'{}'", &self.0)
  }
}

#[derive(Debug)]
pub struct Nside(u32);
impl Nside {
  pub fn get(&self) -> u32 {
    self.0
  }
}
impl FitsCard for Nside {
  const KEYWORD: &'static [u8; 8] = b"NSIDE   ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u32>(keyword_record).map(Self)
  }
  fn to_fits_value(&self) -> String {
    format!("{}", &self.0)
  }
}

#[derive(Debug)]
pub struct FirstPix(u64);
impl FirstPix {
  pub fn get(&self) -> u64 {
    self.0
  }
}
impl Default for FirstPix {
  fn default() -> Self {
    Self(0)
  }
}
impl FitsCard for FirstPix {
  const KEYWORD: &'static [u8; 8] = b"FIRSTPIX";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u64>(keyword_record).map(Self)
  }
  fn to_fits_value(&self) -> String {
    format!("{}", &self.0)
  }
}

#[derive(Debug)]
pub struct LastPix(u64);
impl LastPix {
  pub fn get(&self) -> u64 {
    self.0
  }
}
impl FitsCard for LastPix {
  const KEYWORD: &'static [u8; 8] = b"LASTPIX ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u64>(keyword_record).map(Self)
  }
  fn to_fits_value(&self) -> String {
    format!("{}", &self.0)
  }
}

// Healpix map specific
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum IndexSchema {
  Implicit,
  Explicit,
  Sparse,
}
impl FitsCard for IndexSchema {
  const KEYWORD: &'static [u8; 8] = b"INDXSCHM";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"IMPLICIT" => Ok(IndexSchema::Implicit),
      b"EXPLICIT" => Ok(IndexSchema::Explicit),
      b"SPARSE" => Ok(IndexSchema::Sparse),
      parsed_val => Err(Self::predefine_val_err(
        parsed_val,
        &[b"IMPLICIT", b"EXPLICIT", b"SPARSE"],
      )),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      IndexSchema::Implicit => "'IMPLICIT'",
      IndexSchema::Explicit => "'EXPLICIT'",
      IndexSchema::Sparse => "'SPARSE'",
    })
  }
}

// Uses the index in an array of Option(SkymapKeywords) for fast retrieving of the Card :)
pub trait SkymapCard: FitsCard {
  const INDEX: u8;
}
impl SkymapCard for Ordering {
  const INDEX: u8 = 0;
}
impl SkymapCard for CoordSys {
  const INDEX: u8 = 1;
}
impl SkymapCard for Order {
  const INDEX: u8 = 2;
}
impl SkymapCard for PixType {
  const INDEX: u8 = 3;
}
impl SkymapCard for TForm1 {
  const INDEX: u8 = 4;
}
impl SkymapCard for TType1 {
  const INDEX: u8 = 5;
}
impl SkymapCard for Nside {
  const INDEX: u8 = 6;
}
impl SkymapCard for FirstPix {
  const INDEX: u8 = 7;
}
impl SkymapCard for LastPix {
  const INDEX: u8 = 8;
}
impl SkymapCard for IndexSchema {
  const INDEX: u8 = 9;
}

#[derive(Debug)]
pub(super) struct SkymapKeywordsMap {
  entries: [Option<SkymapKeywords>; 10],
}
impl SkymapKeywordsMap {
  pub(super) fn new() -> SkymapKeywordsMap {
    Self {
      entries: [None, None, None, None, None, None, None, None, None, None],
    }
  }

  pub(super) fn insert(&mut self, entry: SkymapKeywords) -> Option<SkymapKeywords> {
    self.entries[entry.index()].replace(entry)
  }

  pub(super) fn get<T: SkymapCard>(&self) -> Option<&SkymapKeywords> {
    self.entries[T::INDEX as usize].as_ref()
  }

  #[allow(dead_code)]
  pub(super) fn write_all(&self, keyword_records: &mut ChunksMut<u8>) -> Result<(), FitsError> {
    for kw in self.entries.iter().filter_map(|v| v.as_ref()) {
      kw.write_keyword_record(keyword_records.next().unwrap())?;
    }
    Ok(())
  }

  pub(super) fn check_pixtype(&self) -> Result<(), FitsError> {
    match self.get::<PixType>() {
      Some(SkymapKeywords::PixType(PixType::Healpix)) => Ok(()),
      None => Err(FitsError::MissingKeyword {
        keyword: PixType::keyword_string(),
      }),
      _ => unreachable!(), // since there is only one elem in PixType enum
    }
  }

  pub(super) fn check_coordsys(
    &self,
    expected: CoordSys,
    accept_not_found: bool,
  ) -> Result<(), FitsError> {
    match self.get::<CoordSys>() {
      Some(SkymapKeywords::CoordSys(actual)) => {
        if *actual == expected {
          Ok(())
        } else {
          Err(FitsError::UnexpectedValue {
            keyword: CoordSys::keyword_string(),
            expected: CoordSys::Cel.to_fits_value(),
            actual: actual.to_fits_value(),
          })
        }
      }
      None => {
        if accept_not_found {
          warn!(
            "Missing keyword '{}'; Value '{}' is assumed!",
            CoordSys::keyword_str(),
            expected.to_fits_value()
          );
          Ok(())
        } else {
          Err(FitsError::MissingKeyword {
            keyword: CoordSys::keyword_string(),
          })
        }
      }
      _ => unreachable!(), // since there is only one elem in CoorSys enum
    }
  }

  pub(super) fn check_ordering(
    &self,
    expected: Ordering,
    accept_not_found: bool,
  ) -> Result<(), FitsError> {
    match self.get::<Ordering>() {
      Some(SkymapKeywords::Ordering(actual)) => {
        if *actual == expected {
          Ok(())
        } else {
          Err(FitsError::UnexpectedValue {
            keyword: Ordering::keyword_string(),
            expected: expected.to_fits_value(),
            actual: actual.to_fits_value(),
          })
        }
      }
      _ => {
        if accept_not_found {
          warn!(
            "Missing keyword '{}'; Value '{}' is assumed!",
            Ordering::keyword_str(),
            expected.to_fits_value()
          );
          Ok(())
        } else {
          Err(FitsError::MissingKeyword {
            keyword: Ordering::keyword_string(),
          })
        }
      }
    }
  }

  pub(super) fn check_firstpix(
    &self,
    expected: u64,
    accept_not_found: bool,
  ) -> Result<(), FitsError> {
    match self.get::<FirstPix>() {
      Some(SkymapKeywords::FirstPix(FirstPix(actual))) => {
        if *actual == expected {
          Ok(())
        } else {
          Err(FitsError::UnexpectedValue {
            keyword: FirstPix::keyword_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
          })
        }
      }
      None => {
        if accept_not_found {
          warn!(
            "Missing keyword '{}'; Value '{}' is assumed!",
            FirstPix::keyword_str(),
            expected
          );
          Ok(())
        } else {
          Err(FitsError::MissingKeyword {
            keyword: FirstPix::keyword_string(),
          })
        }
      }
      _ => unreachable!(),
    }
  }

  pub(super) fn check_lastpix(
    &self,
    expected: u64,
    accept_not_found: bool,
  ) -> Result<(), FitsError> {
    match self.get::<LastPix>() {
      Some(SkymapKeywords::LastPix(LastPix(actual))) => {
        if *actual == expected {
          Ok(())
        } else {
          Err(FitsError::UnexpectedValue {
            keyword: LastPix::keyword_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
          })
        }
      }
      None => {
        if accept_not_found {
          warn!(
            "Missing keyword '{}'; Value '{}' is assumed!",
            LastPix::keyword_str(),
            &expected
          );
          Ok(())
        } else {
          Err(FitsError::MissingKeyword {
            keyword: LastPix::keyword_string(),
          })
        }
      }
      _ => unreachable!(),
    }
  }

  pub(super) fn check_index_schema(&self, expected: IndexSchema) -> Result<(), FitsError> {
    match self.get::<IndexSchema>() {
      Some(SkymapKeywords::IndexSchema(actual)) => {
        if *actual == expected {
          Ok(())
        } else {
          Err(FitsError::UnexpectedValue {
            keyword: IndexSchema::keyword_string(),
            expected: expected.to_fits_value(),
            actual: actual.to_fits_value(),
          })
        }
      }
      _ => Err(FitsError::MissingKeyword {
        keyword: IndexSchema::keyword_string(),
      }),
    }
  }
}

#[derive(Debug)]
pub enum SkymapKeywords {
  Ordering(Ordering),
  CoordSys(CoordSys),
  Order(Order),
  PixType(PixType),
  TForm1(TForm1),
  TType1(TType1),
  Nside(Nside),
  FirstPix(FirstPix),
  LastPix(LastPix),
  IndexSchema(IndexSchema),
}
impl SkymapKeywords {
  pub(super) fn is_skymap_kw(keyword_record: &[u8]) -> Result<Option<Self>, FitsError> {
    // I have not yet found how to match on the FitsCard::KEYWORD associated constant :o/
    match get_keyword(keyword_record) {
      b"ORDERING" => Ordering::parse_value(keyword_record)
        .map(SkymapKeywords::Ordering)
        .map(Some),
      b"COORDSYS" => CoordSys::parse_value(keyword_record)
        .map(SkymapKeywords::CoordSys)
        .map(Some),
      b"ORDER   " => Order::parse_value(keyword_record)
        .map(SkymapKeywords::Order)
        .map(Some),
      b"PIXTYPE " => PixType::parse_value(keyword_record)
        .map(SkymapKeywords::PixType)
        .map(Some),
      b"TFORM1  " => TForm1::parse_value(keyword_record)
        .map(SkymapKeywords::TForm1)
        .map(Some),
      b"TTYPE1  " => TType1::parse_value(keyword_record)
        .map(SkymapKeywords::TType1)
        .map(Some),
      b"NSIDE   " => Nside::parse_value(keyword_record)
        .map(SkymapKeywords::Nside)
        .map(Some),
      b"FIRSTPIX" => FirstPix::parse_value(keyword_record)
        .map(SkymapKeywords::FirstPix)
        .map(Some),
      b"LASTPIX " => LastPix::parse_value(keyword_record)
        .map(SkymapKeywords::LastPix)
        .map(Some),
      b"INDXSCHM" => IndexSchema::parse_value(keyword_record)
        .map(SkymapKeywords::IndexSchema)
        .map(Some),
      _ => Ok(None),
    }
  }

  fn index(&self) -> usize {
    (match self {
      SkymapKeywords::Ordering(_) => Ordering::INDEX,
      SkymapKeywords::CoordSys(_) => CoordSys::INDEX,
      SkymapKeywords::Order(_) => Order::INDEX,
      SkymapKeywords::PixType(_) => PixType::INDEX,
      SkymapKeywords::TForm1(_) => TForm1::INDEX,
      SkymapKeywords::TType1(_) => TType1::INDEX,
      SkymapKeywords::Nside(_) => Nside::INDEX,
      SkymapKeywords::FirstPix(_) => FirstPix::INDEX,
      SkymapKeywords::LastPix(_) => LastPix::INDEX,
      SkymapKeywords::IndexSchema(_) => IndexSchema::INDEX,
    }) as usize
  }

  pub(super) fn keyword(&self) -> &'static [u8; 8] {
    match self {
      SkymapKeywords::Ordering(_) => Ordering::KEYWORD,
      SkymapKeywords::CoordSys(_) => CoordSys::KEYWORD,
      SkymapKeywords::Order(_) => Order::KEYWORD,
      SkymapKeywords::PixType(_) => PixType::KEYWORD,
      SkymapKeywords::TForm1(_) => TForm1::KEYWORD,
      SkymapKeywords::TType1(_) => TType1::KEYWORD,
      SkymapKeywords::Nside(_) => Nside::KEYWORD,
      SkymapKeywords::FirstPix(_) => FirstPix::KEYWORD,
      SkymapKeywords::LastPix(_) => LastPix::KEYWORD,
      SkymapKeywords::IndexSchema(_) => IndexSchema::KEYWORD,
    }
  }

  pub(super) fn keyword_str(&self) -> &str {
    unsafe { str::from_utf8_unchecked(self.keyword()) }.trim_end()
  }

  fn write_keyword_record(&self, keyword_record: &mut [u8]) -> Result<(), FitsError> {
    match self {
      SkymapKeywords::Ordering(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::CoordSys(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::Order(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::PixType(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::TForm1(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::TType1(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::Nside(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::FirstPix(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::LastPix(kw) => kw.write_keyword_record(keyword_record),
      SkymapKeywords::IndexSchema(kw) => kw.write_keyword_record(keyword_record),
    }
  }
}
