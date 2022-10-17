//! Module defining the commodity enums Cardinal (points) and MainWind (points). 

use std::mem;

/// Cardinal points 
#[derive(Debug, PartialEq, Eq)]
pub enum Cardinal {
  /// South direction
  S,
  /// East direction
  E,
  /// North direction
  N,
  /// West direction
  W
}

impl Cardinal {

  /// Returns a Cardinal point give an index.
  /// We define this method for calls from external languages
  /// - 0 => S
  /// - 1 => E
  /// - 2 => N
  /// - 3 => W
  /// 
  /// # Input
  /// - `i` index in `[0, 3]`
  /// 
  /// # Output
  /// - cardinal point
  /// 
  /// # Panics
  /// If the given index in not in `[0, 3]`
  /// 
  /// # Example
  /// ```rust
  /// use cdshealpix::compass_point::{Cardinal};
  /// 
  /// assert_eq!(Cardinal::from_index(0), Cardinal::S);
  /// assert_eq!(Cardinal::from_index(1), Cardinal::E);
  /// assert_eq!(Cardinal::from_index(2), Cardinal::N);
  /// assert_eq!(Cardinal::from_index(3), Cardinal::W);
  /// ```
  pub fn from_index(i: u8) -> Cardinal {
    match i {
      0 => Cardinal::S,
      1 => Cardinal::E,
      2 => Cardinal::N,
      3 => Cardinal::W,
      _ => panic!("Wrong Cardinal index: Expected value in [0, 3]; Actual value {}.", i),
    }
  }

  pub(crate) fn index(&self) -> u8 {
    match *self {
      Cardinal::S => 0,
      Cardinal::E => 1,
      Cardinal::N => 2,
      Cardinal::W => 3,
    }
  }
  
  /// Returns:
  /// - South: -offset
  /// - North: +offset
  /// - East and West: 0
  pub(super) fn offset_sn(&self, offset: f64) -> f64 {
    match *self {
      Cardinal::S => -offset,
      Cardinal::N => offset,
      _ => 0f64,
    }
  }

  /// Returns:
  /// - West: -offset
  /// - East: +offset
  /// - South and North: 0
  pub(super) fn offset_we(&self, offset: f64) -> f64 {
    match *self {
      Cardinal::W => -offset,
      Cardinal::E => offset,
      _ => 0f64,
    }
  }

  #[allow(dead_code)]
  pub(super) fn next_clockwise(&self) -> Cardinal {
    match *self {
      Cardinal::S => Cardinal::W,
      Cardinal::E => Cardinal::S,
      Cardinal::N => Cardinal::E,
      Cardinal::W => Cardinal::N,
    }
  }

  #[allow(dead_code)]
  pub(super) fn next_counter_clockwise(&self) -> Cardinal {
    match *self {
      Cardinal::S => Cardinal::E,
      Cardinal::E => Cardinal::N,
      Cardinal::N => Cardinal::W,
      Cardinal::W => Cardinal::S,
    }
  }

  pub(super) fn clockwise_cycle(&self) -> (Cardinal, Cardinal, Cardinal, Cardinal) {
    match *self {
      Cardinal::S => (Cardinal::S, Cardinal::W, Cardinal::N, Cardinal::E),
      Cardinal::E => (Cardinal::E, Cardinal::S, Cardinal::W, Cardinal::N),
      Cardinal::N => (Cardinal::N, Cardinal::E, Cardinal::S, Cardinal::W),
      Cardinal::W => (Cardinal::W, Cardinal::N, Cardinal::E, Cardinal::S),
    }
  }

  pub(super) fn counter_clockwise_cycle(&self) -> (Cardinal, Cardinal, Cardinal, Cardinal) {
    match *self {
      Cardinal::S => (Cardinal::S, Cardinal::E, Cardinal::N, Cardinal::W),
      Cardinal::E => (Cardinal::E, Cardinal::N, Cardinal::W, Cardinal::S),
      Cardinal::N => (Cardinal::N, Cardinal::W, Cardinal::S, Cardinal::E),
      Cardinal::W => (Cardinal::W, Cardinal::S, Cardinal::E, Cardinal::N),
    }
  }
}

/// Cardinal set. 
/// Internally the information is stored on 4 bits.
pub struct CardinalSet {
  // Byte storing the 4 boolean values
  byte: u8,
}

impl Default for CardinalSet {
  fn default() -> Self {
    Self::new()
  }
}

impl CardinalSet {
  /// Returns an empty cardinal set
  pub fn new() -> CardinalSet {
    CardinalSet {
      byte: 0_u8,
    }
  }
  
  /// Retusn a cardinal set with all directions set
  pub fn all() -> CardinalSet {
    CardinalSet {
      byte: 0b00001111_u8,
    }
  }
  
  /// Add or remove (or do nothing) the given direction to the set.
  pub fn set(&mut self, key: Cardinal, value: bool) {
    let i = key.index() as u8;
    let mask = 1u8 << i;
    if value {
      self.byte |= mask;  
    } else {
      self.byte &= !mask;
    }
  }
  
  /// Returns `true` if the given direction is in the set.
  pub fn get(&self, key: Cardinal) -> bool {
    let i = key.index();
    self.get_from_index(i)
  }
  
  fn get_from_index(&self, index: u8) -> bool {
    let mask = 1u8 << index;
    self.byte & mask != 0u8
  }
  
  /// Remove all directions from the set
  pub fn clear(&mut self) {
    self.byte = 0u8;
  }
}

impl IntoIterator for CardinalSet {
  type Item = Cardinal;
  type IntoIter = CardinalSetIterator;
  
  fn into_iter(self) -> Self::IntoIter {
    CardinalSetIterator {
      cardinal_set: self,
      index: 0
    }
  }
}

/// Structure used to iterate over a CardinalSet
pub struct CardinalSetIterator {
  cardinal_set: CardinalSet,
  index: u8,
}

impl Iterator for CardinalSetIterator {
  type Item = Cardinal;
  
  fn next(&mut self) -> Option<Cardinal> {
    while self.index < 4 {
      if self.cardinal_set.get_from_index(self.index) {
        let card = Cardinal::from_index(self.index);
        self.index += 1;
        return Some(card)
      } else {
        self.index += 1;
      }
    }
    None
  }
}

/// Equivalent of a Java EnumMap for cardinal directions.
/// We require T to implement the Copy trait since internally we use an array stored on the stack.
pub struct CardinalMap<T: Copy> {
  array: [Option<T>; 4],
}

impl<T: Copy> Default for CardinalMap<T> {
  fn default() -> Self {
    Self::new()
  }
}

/// Equivalent of a Java EnumMap for the cardinal points.
/// We require T to implement the Copy trait since internally we use an array stored on the stack.
impl<V: Copy> CardinalMap<V> {

  /// Creates a new empty map.
  pub fn new() ->  CardinalMap<V> {
    CardinalMap {
      array: [Option::None; 4],
    }
  }

  /// Associate the given value with the given direction
  pub fn put(&mut self, key: Cardinal, value: V) -> Option<V> {
    mem::replace(&mut self.array[key.index() as usize], Some(value))
  }

  /// Get a pointer to the value associated with the given direction
  pub fn get(&self, key: Cardinal) -> Option<&V> {
    self.array[key.index() as usize].as_ref()
  }
  
  /// Replace all values by None
  pub fn clear(&mut self) {
    // fill_with_none(self.array);
    for i in 0..4 {
      self.array[i] = None;
    }
  }
}

/////////////
// Ordinal //
/////////////

/// Cardinal points 
#[derive(Debug, PartialEq, Eq)]
pub enum Ordinal {
  /// Southeast direction
  SE,
  /// Southwest direction
  SW,
  /// Northeast direction
  NE,
  /// Northwest direction
  NW
}

impl Ordinal {
  pub(crate) fn index( & self ) -> u8 {
    match * self {
      Ordinal::SE => 0,
      Ordinal::SW => 1,
      Ordinal::NE => 2,
      Ordinal::NW => 3,
    }
  }
}

/// Equivalent of a Java EnumMap for ordinal directions.
/// We require T to implement the Copy trait since internally we use an array stored on the stack.
pub struct OrdinalMap<T: Copy> {
  array: [Option<T>; 4],
}

impl<T: Copy> Default for OrdinalMap<T> {
  fn default() -> Self {
    Self::new()
  }
}

/// Equivalent of a Java EnumMap for the cardinal points.
/// We require T to implement the Copy trait since internally we use an array stored on the stack.
impl<V: Copy> OrdinalMap<V> {

  /// Creates a new empty map.
  pub fn new() ->  OrdinalMap<V> {
    OrdinalMap {
      array: [Option::None; 4],
    }
  }

  /// Associate the given value with the given direction
  pub fn put(&mut self, key: Ordinal, value: V) -> Option<V> {
    mem::replace(&mut self.array[key.index() as usize], Some(value))
  }

  /// Get a pointer to the value associated with the given direction
  pub fn get(&self, key: Ordinal) -> Option<&V> {
    self.array[key.index() as usize].as_ref()
  }

  /// Replace all values by None
  pub fn clear(&mut self) {
    for i in 0..4 {
      self.array[i] = None;
    }
  }
}

///////////////
// MAIN WIND //
///////////////

/// Main winds directions
#[derive(Debug, PartialEq, Eq)]
pub enum MainWind {
  /// South
  S,
  /// Southeast
  SE,
  /// East
  E,
  /// Southwest
  SW,
  /// Center (not a real main winds)
  C,
  /// Northeast
  NE,
  /// West
  W,
  /// Norhtwest
  NW,
  /// North
  N
}

impl MainWind {

  /// Returns a Main wind direction give an index.
  /// We define this method for calls from external languages
  /// - 0 => S
  /// - 1 => SE
  /// - 2 => E
  /// - 3 => SW
  /// - 4 => C
  /// - 5 => NE
  /// - 6 => W
  /// - 7 => NW
  /// - 8 => N
  /// 
  /// # Input
  /// - `i` index in `[0, 8]`
  /// 
  /// # Output
  /// - main wind direction
  /// 
  /// # Panics
  /// If the given index in not in `[0, 8]`
  /// 
  /// # Example
  /// 
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// 
  /// assert_eq!(MainWind::from_index(0), MainWind::S);
  /// assert_eq!(MainWind::from_index(1), MainWind::SE);
  /// assert_eq!(MainWind::from_index(2), MainWind::E);
  /// assert_eq!(MainWind::from_index(3), MainWind::SW);
  /// assert_eq!(MainWind::from_index(4), MainWind::C);
  /// assert_eq!(MainWind::from_index(5), MainWind::NE);
  /// assert_eq!(MainWind::from_index(6), MainWind::W);
  /// assert_eq!(MainWind::from_index(7), MainWind::NW);
  /// assert_eq!(MainWind::from_index(8), MainWind::N);
  /// ```
  pub fn from_index(i: u8) -> MainWind {
    match i {
      0 => MainWind::S,
      1 => MainWind::SE,
      2 => MainWind::E,
      3 => MainWind::SW,
      4 => MainWind::C,
      5 => MainWind::NE,
      6 => MainWind::W,
      7 => MainWind::NW,
      8 => MainWind::N,
      _ => panic!("Wrong MainWind index: Expected value in [0, 7]; Actual value {}.", i),
    }
  }

  /// Returns `true` is the Main Wind is more particularly a Cardinal point.
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// 
  /// assert_eq!(true, MainWind::S.is_cardinal());
  /// assert_eq!(true, MainWind::E.is_cardinal());
  /// assert_eq!(true, MainWind::N.is_cardinal());
  /// assert_eq!(true, MainWind::W.is_cardinal());
  /// assert_eq!(false, MainWind::C.is_cardinal());
  /// assert_eq!(false, MainWind::SE.is_cardinal());
  /// assert_eq!(false, MainWind::SW.is_cardinal());
  /// assert_eq!(false, MainWind::NE.is_cardinal());
  /// assert_eq!(false, MainWind::NW.is_cardinal());
  /// ```
  pub fn is_cardinal(&self) -> bool {
    matches!(*self, MainWind::S | MainWind::E | MainWind::N | MainWind::W)
  }

  /// Convert this main wind into a Cardinal point.
  /// # Panics
  /// If the maind wind is not a cardinal point
  pub fn to_cardinal(&self) -> Cardinal {
    // use self::Cardinal;
    match *self {
      MainWind::S => Cardinal::S,
      MainWind::E => Cardinal::E,
      MainWind::N => Cardinal::N,
      MainWind::W => Cardinal::W,
      _ => panic!("Main wind '{:?}' can't be converted to cardinal!", &self),
    }
  }

  /// Returns `true` is the Main Wind is more particularly an Ordinal point.
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// 
  /// assert_eq!(false, MainWind::S.is_ordinal());
  /// assert_eq!(false, MainWind::E.is_ordinal());
  /// assert_eq!(false, MainWind::N.is_ordinal());
  /// assert_eq!(false, MainWind::W.is_ordinal());
  /// assert_eq!(false, MainWind::C.is_ordinal());
  /// assert_eq!(true, MainWind::SE.is_ordinal());
  /// assert_eq!(true, MainWind::SW.is_ordinal());
  /// assert_eq!(true, MainWind::NE.is_ordinal());
  /// assert_eq!(true, MainWind::NW.is_ordinal());
  /// ```
  pub fn is_ordinal(&self) -> bool {
    matches!(*self, MainWind::SE | MainWind::SW | MainWind::NE | MainWind::NW)
  }

  /// Convert this main wind into an Ordinal point.
  /// # Panics
  /// If the maind wind is not an orinal point
  pub fn to_ordinal(&self) -> Ordinal {
    match *self {
      MainWind::SE => Ordinal::SE,
      MainWind::SW => Ordinal::SW,
      MainWind::NE => Ordinal::NE,
      MainWind::NW => Ordinal::NW,
      _ => panic!("Main wind '{:?}' can't be converted to ordinal!", &self),
    }
  }
  
  /// Returns the given main wind opposite direction.
  /// 
  /// # Example
  /// 
  /// ```rust
  /// use cdshealpix::compass_point::{MainWind};
  /// 
  /// assert_eq!(MainWind::S.opposite(),  MainWind::N);
  /// assert_eq!(MainWind::SE.opposite(), MainWind::NW);
  /// assert_eq!(MainWind::E.opposite(),  MainWind::W);
  /// assert_eq!(MainWind::SW.opposite(), MainWind::NE);
  /// assert_eq!(MainWind::C.opposite(),  MainWind::C);
  /// assert_eq!(MainWind::NE.opposite(), MainWind::SW);
  /// assert_eq!(MainWind::W.opposite(),  MainWind::E);
  /// assert_eq!(MainWind::NW.opposite(), MainWind::SE);
  /// assert_eq!(MainWind::N.opposite(),  MainWind::S);
  /// ```
  pub fn opposite(&self) -> MainWind {
    match *self {
        MainWind::S  => MainWind::N,
        MainWind::SE => MainWind::NW,
        MainWind::E  => MainWind::W,
        MainWind::SW => MainWind::NE,
        MainWind::C  => MainWind::C,
        MainWind::NE => MainWind::SW,
        MainWind::W  => MainWind::E,
        MainWind::NW => MainWind::SE,
        MainWind::N  => MainWind::S,
    }
  }
  
  /// Returns the given Main Wind direction according to the given offsets.
  /// - `offset_se must` be in `[-1, 1]`
  /// - `offset_sw must` be in `[-1, 1]`
  /// 
  /// # Example
  /// 
  /// ```rust
  /// use cdshealpix::compass_point::MainWind;
  /// use cdshealpix::compass_point::MainWind::{S, SE, E, SW, C, NE, W, NW, N};
  /// 
  /// assert_eq!(MainWind::from_offsets(-1, -1),  S);
  /// assert_eq!(MainWind::from_offsets( 0, -1), SE);
  /// assert_eq!(MainWind::from_offsets( 1, -1),  E);
  /// assert_eq!(MainWind::from_offsets(-1,  0), SW);
  /// assert_eq!(MainWind::from_offsets( 0,  0),  C);
  /// assert_eq!(MainWind::from_offsets( 1,  0), NE);
  /// assert_eq!(MainWind::from_offsets(-1,  1),  W);
  /// assert_eq!(MainWind::from_offsets( 0,  1), NW);
  /// assert_eq!(MainWind::from_offsets( 1,  1),  N);
  /// ```
  pub fn from_offsets(offset_se: i8, offset_sw: i8) -> MainWind {
    debug_assert!((-1_i8..=1_i8).contains(&offset_se));
    debug_assert!((-1_i8..=1_i8).contains(&offset_sw));
    let mut i = (offset_sw + 1_i8) as u8;
    i += (i << 1) + (offset_se + 1_i8) as u8;
    MainWind::from_index(i)
  }
  
  fn index(&self) -> u8 {
    match *self {
      MainWind::S  => 0,
      MainWind::SE => 1,
      MainWind::E  => 2,
      MainWind::SW => 3,
      MainWind::C  => 4,
      MainWind::NE => 5,
      MainWind::W  => 6,
      MainWind::NW => 7,
      MainWind::N  => 8,
    }
  }

  /// Returns:
  /// _W NW _N
  /// SW _C NE
  /// _S SE _E 
  /// ----------> SE
  /// -1  0  1
  pub(super) fn offset_se(&self) -> i8 {
    match *self {
      MainWind::S  => -1,
      MainWind::SE =>  0,
      MainWind::E  =>  1,
      MainWind::SW => -1,
      MainWind::C  =>  0,
      MainWind::NE =>  1,
      MainWind::W  => -1,
      MainWind::NW =>  0,
      MainWind::N  =>  1,
    }
  }

  /// Returns:
  ///    ^
  ///  1 | _W NW _N
  ///  0 | SW _C NE
  /// -1 | _S SE _E
  pub(super) fn offset_sw(&self) -> i8 {
    match *self {
      MainWind::S  => -1,
      MainWind::SE => -1,
      MainWind::E  => -1,
      MainWind::SW => 0,
      MainWind::C  => 0,
      MainWind::NE => 0,
      MainWind::W  => 1,
      MainWind::NW => 1,
      MainWind::N  => 1,
    }
  }
}

/// Equivalent of a Java EnumMap for the main winds.
/// We require T to implement the Copy trait since internally we use an array stored on the stack.
#[derive(Debug)]
pub struct MainWindMap<T: Copy> {
  array: [Option<T>; 9],
}

impl<T: Copy> Default for MainWindMap<T> {
  fn default() -> Self {
    Self::new()
  }
}

impl<V: Copy> MainWindMap<V> {

  /// Creates a new empty map.
  pub fn new() -> MainWindMap<V> {
    MainWindMap {
      array: [Option::None; 9],
    }
  }

  /// Associate the given value with the given direction
  pub fn put(&mut self, key: MainWind, value: V) -> Option<V> {
    mem::replace(&mut self.array[key.index() as usize], Some(value))
  }

  /// Associate None with the given direction
  pub fn put_none(&mut self, key: MainWind) -> Option<V> {
    mem::replace(&mut self.array[key.index() as usize], None)
  }

  /// Associate the given Option with the given direction
  pub fn put_opt(&mut self, key: MainWind, value: Option<V>) -> Option<V> {
    mem::replace(&mut self.array[key.index() as usize], value)
  }

  /// Get a pointer to the value associated with the given direction
  pub fn get(&self, key: MainWind) -> Option<&V> {
    self.get_from_index(key.index() as usize)
  }

  fn get_from_index(&self, index: usize) -> Option<&V> {
    self.array[index].as_ref()
  }
  
  /// Replace all values by None 
  pub fn clear(&self) {
    let mut a = self.array; 
    fill_with_none(&mut a);
  }

  /// Returns a vector of values
  pub fn values_vec(&self) -> Vec<V> {
    self.array.iter().filter_map(|&o| o).collect::<Vec<V>>()
  }
  
  pub fn entries(&self) -> Box<[(MainWind, V)]> {
    self.entries_vec().into_boxed_slice()
  }

  pub fn entries_vec(&self) -> Vec<(MainWind, V)> {
    self.array.iter().enumerate()
      .filter_map(|(i, &o)| o.map(|v| (MainWind::from_index(i as u8), v)))
      .collect::<Vec<(MainWind, V)>>()
  }
  
}

impl<V: Copy + Ord> MainWindMap<V> {
  
  /// Returns the values contained in the map, ordered in their natural order in a fixed length array
  pub fn sorted_values(&self) -> Box<[V]> {
    //let mut values = self.values_vec().into_boxed_slice();
    //values.sort_unstable();
    //values
    self.sorted_values_vec().into_boxed_slice()
  }

  /// Returns the values contained in the map, ordered in their natural order in a growable array
  pub fn sorted_values_vec(&self) -> Vec<V> {
    let mut values: Vec<V> = self.values_vec();
    values.sort_unstable();
    values
  }

  pub fn sorted_entries(&self) -> Box<[(MainWind, V)]> {
    self.sorted_entries_vec().into_boxed_slice()
  }
  
  pub fn sorted_entries_vec(&self) -> Vec<(MainWind, V)> {
    let mut entries = self.entries_vec();
    entries.sort_unstable_by(|(_, v1), (_, v2)| v1.partial_cmp(v2).unwrap());
    entries
  }
}


/// Fill the given (mutable reference on a) slice of Option with the None value
fn fill_with_none<T>(array_of_option: &mut [Option<T>]) {
  array_of_option.iter_mut().for_each(|o| *o = None);
}
