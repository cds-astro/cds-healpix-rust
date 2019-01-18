//! Module defining the commodity enums Cardinal (points) and MainWind (points). 

use std::mem;

/// Cardinal points 
#[derive(Debug, PartialEq)]
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

  fn index(&self) -> u8 {
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
}

/// Cardinal set. 
/// Internally the information is stored on 4 bits.
pub struct CardinalSet {
  // Byte storing the 4 boolean values
  byte: u8,
}

impl CardinalSet {
  
  /// Retuns a new empty cardinal set.
  pub fn new() -> CardinalSet {
    CardinalSet {
      byte: 0u8,
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
    match self.array[key.index() as usize] {
      Some(ref v) => Some(v),
      None => None,
    }
  }
  
  /// Replace all values by None
  pub fn clear(&mut self) {
    // fill_with_none(self.array);
    for i in 0..4 {
      self.array[i] = None;
    }
  }
}


///////////////
// MAIN WIND //
///////////////

/// Main winds directions
#[derive(Debug, PartialEq)]
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
    debug_assert!(-1_i8 <= offset_se && offset_se <= 1_i8);
    debug_assert!(-1_i8 <= offset_sw && offset_sw <= 1_i8);
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
    match self.array[index] {
      Some(ref v) => Some(v),
      None => None,
    }
  }
  
  /// Replace all values by None 
  pub fn clear(&self) {
    let mut a = self.array; 
    fill_with_none(&mut a);
  }

  /// Returns a vector of values
  pub fn values_vec(&self) -> Vec<V> {
    self.array.into_iter().filter_map(|&o| o).collect::<Vec<V>>()
  }
  
}

impl<V: Copy + Ord> MainWindMap<V> {
  
  /// Returns the values contained in the map, ordered in their natural order in a fixed length array
  pub fn sorted_values(&self) -> Box<[V]> {
    let mut values = self.values_vec().into_boxed_slice();
    values.sort_unstable();
    values
  }

  /// Returns the values contained in the map, ordered in their natural order in a growable array
  pub fn sorted_values_vec(&self) -> Vec<V> {
    let mut values: Vec<V> = self.values_vec();
    values.sort_unstable();
    values
  }
}


/// Fill the given (mutable reference on a) slice of Option with the None value
fn fill_with_none<T>(array_of_option: &mut [Option<T>]) {
  array_of_option.iter_mut().for_each(|o| *o = None);
}
