//! Definition of a BMOC, i.e. a MOC storing an additional flag telling if a cell is fully
//! or partially covered by the MOC.

use std::slice::Iter;

use super::super::{nside_square_unsafe};

/// A very basic and simple BMOC Builder: we push elements in it assuming that we provide them
/// in the write order, without duplicates and small cells included in larger cells.
pub(super) struct BMOCBuilderUnsafe {
  depth_max: u8,
  entries: Option<Vec<u64>>,
}

impl BMOCBuilderUnsafe {
  
  pub fn new(depth_max: u8, capacity: usize) -> BMOCBuilderUnsafe {
    BMOCBuilderUnsafe {
      depth_max,
      entries: Some(Vec::with_capacity(capacity)),
    }
  }
  
  /* Commented because not used so far
  /// Clear the content and start a fresh builder with the given initial capacity.
  #[warn(dead_code)]
  pub fn re_init(&mut self, capacity: usize) -> &mut BMOCBuilderUnsafe {
    self.entries = Some(Vec::with_capacity(capacity));
    self
  }*/
  
  pub fn push(&mut self, depth: u8, hash: u64, is_full: bool) -> &mut BMOCBuilderUnsafe {
    if let Some(ref mut v) = self.entries {
      v.push(build_raw_value(depth, hash, is_full, self.depth_max));
    } else {
      panic!("Empty builder, you have to re-init it before re-using it!");
    }
    self
  }
  
  pub fn push_all(&mut self, depth: u8, from_hash: u64, to_hash: u64, are_full: bool) -> &mut BMOCBuilderUnsafe {
    if let Some(ref mut v) = self.entries {
      for h in from_hash..to_hash {
        v.push(build_raw_value(depth, h, are_full, self.depth_max));
      }
    } else {
      panic!("Empty builder, you have to re-init it before re-using it!");
    }
    self
  }
  
  pub fn to_bmoc(&mut self) -> BMOC {
    BMOC::create_unsafe(self.depth_max, self.entries.take().expect("Empty builder!").into_boxed_slice())
  }

  fn pack(&mut self) -> Vec<u64> {
    let mut entries = self.entries.take().expect("Empty builder!");
    // On-place pack
    let mut prev_to_index = 0_usize;
    let mut curr_to_index = entries.len();
    while prev_to_index != curr_to_index { // changes occurs
      prev_to_index = curr_to_index;
      let mut i_prev_moc = 0_usize;
      let mut i_curr_moc = 0_usize;
      while i_prev_moc < prev_to_index {
        let mut curr_cell = entries[i_prev_moc];
        i_prev_moc += 1;
        let mut curr_cell_depth = get_depth(curr_cell, self.depth_max);
        let mut curr_cell_hash = get_hash_from_delta_depth(curr_cell, self.depth_max - curr_cell_depth);
        // Look for the first cell of the larger cell (depth - 1)  (=> 2 last bits = 00), the cell must be FULL
        while i_prev_moc < prev_to_index &&
          (curr_cell_depth == 0 || is_partial(curr_cell) || is_not_first_cell_of_larger_cell(curr_cell_hash)) {
          if i_curr_moc != i_prev_moc {
            entries[i_curr_moc] = curr_cell;
            i_curr_moc += 1;
          }
          curr_cell = entries[i_prev_moc];
          i_prev_moc += 1;
          curr_cell_depth = get_depth(curr_cell, self.depth_max);
          curr_cell_hash = get_hash_from_delta_depth(curr_cell, self.depth_max - curr_cell_depth);
        }
        // Look at the 3 siblings
        if i_prev_moc + 2 < prev_to_index
          && entries[i_prev_moc + 0] == build_raw_value(curr_cell_depth, curr_cell_hash | 1, true, self.depth_max)
          && entries[i_prev_moc + 1] == build_raw_value(curr_cell_depth, curr_cell_hash | 2, true, self.depth_max)
          && entries[i_prev_moc + 2] == build_raw_value(curr_cell_depth, curr_cell_hash | 3, true, self.depth_max) {
          entries[i_curr_moc] = build_raw_value(curr_cell_depth - 1, curr_cell_hash >> 2, true, self.depth_max);
          i_curr_moc += 1;
          i_prev_moc += 3;
        } else if i_curr_moc != i_prev_moc {
          entries[i_curr_moc] = curr_cell;
          i_curr_moc += 1;
        }
      }
      curr_to_index = i_curr_moc;
    }
    // We may find a better algorithm doing a single pass on the input MOC
    // Here the number of passes max = mocDepth - smallestDepthOfACellInOutputMoc
    // YEP: new idea: do it like a buffer with a cursor on the last "unmergeable" element!!
    entries.truncate(curr_to_index);
    entries
  }
  
  // We assume the given entries form a valid BMOC (already packef, ordered, ...)
  fn to_lower_depth(&self, new_depth: u8, mut entries: Vec<u64>) -> Vec<u64> {
    if new_depth >= self.depth_max {
      panic!("The given depth must be lower than the depth max of the BMOC");
    }
    let mut i_new = 0_usize;
    let mut prev_hash_at_new_depth = loop {
      if i_new < entries.len() {
        break None;
      }
      let raw_value = entries[i_new];
      let depth = self.get_depth(raw_value);
      if depth <= new_depth {
        let twice_delta_depth = (self.depth_max - new_depth) << 1;
        entries[i_new] = (raw_value >> twice_delta_depth) | (raw_value & 1_u64);
        i_new += 1;
      } else {
        // let twice_delta_depth = (depth - new_depth) << 1;
        break Some(get_hash_from_delta_depth(raw_value, self.depth_max - depth)); // >> twice_delta_depth);
      }
    };
    for i in (i_new + 1)..entries.len() {
      let raw_value = entries[i];
      let depth = self.get_depth(raw_value);
      if depth <= new_depth {
        if prev_hash_at_new_depth.is_some() {
          entries[i_new] = (prev_hash_at_new_depth.take().unwrap() << 2) | 2_u64;
          i_new += 1;
        }
        let twice_delta_depth = (self.depth_max - new_depth) << 1;
        entries[i_new] = (raw_value >> twice_delta_depth) | (raw_value & 1_u64);
        i_new += 1;
      } else {
        let curr_hash_at_new_depth = get_hash_from_delta_depth(raw_value, self.depth_max - new_depth);
        if let Some(prev_val_at_new_depth) = prev_hash_at_new_depth {
          if prev_val_at_new_depth != curr_hash_at_new_depth {
            entries[i_new] = (prev_val_at_new_depth << 2) | 2_u64; // sentinel bit + flag = 0
            i_new += 1;
            prev_hash_at_new_depth.replace(curr_hash_at_new_depth);
          }
        } else {
          prev_hash_at_new_depth.replace(curr_hash_at_new_depth);
        }
      }
    }
    if prev_hash_at_new_depth.is_some() {
      entries[i_new] = (prev_hash_at_new_depth.take().unwrap() << 2) | 2_u64;
      i_new += 1;
    }
    entries.truncate(i_new);
    entries
  }
  
  pub fn to_bmoc_packing(&mut self) -> BMOC {
    let entries = self.pack();
    BMOC::create_unsafe(self.depth_max, entries.into_boxed_slice())
  }

  pub fn to_lower_depth_bmoc_packing(&mut self, new_depth: u8) -> BMOC {
    let entries = self.pack();
    let entries = self.to_lower_depth(new_depth, entries);
    BMOC::create_unsafe(new_depth, entries.into_boxed_slice())
  }

  fn get_depth(&self, raw_value: u64) -> u8 {
    self.get_depth_no_flag(rm_flag(raw_value))
  }
  /// Works both with no flag or with flag set to 0
  fn get_depth_no_flag(&self, raw_value_no_flag: u64) -> u8 {
    self.depth_max - (raw_value_no_flag.trailing_zeros() >> 1) as u8
  }
}

/// Structure defining a simple BMOC.
pub struct BMOC {
  pub depth_max: u8,
  pub entries: Box<[u64]>,
}

#[derive(Debug)]
pub struct Cell {
  pub raw_value: u64,
  pub depth: u8,
  pub hash: u64,
  pub is_full: bool,
}

impl Cell {
  fn new(raw_value: u64, depth_max: u8) -> Cell {
    // Extract the flag
    let is_full = (raw_value | 1_u64) == 1_u64;
    // Remove the flag bit, then divide by 2 (2 bits per level)
    let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
    // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
    let hash = raw_value >> (2 + (delta_depth << 1));
    let depth = depth_max - delta_depth;
    Cell { raw_value, depth, hash, is_full }
  }
}


impl BMOC {
  
  /* Use this for a BMOC builder!
  pub(super) fn new(depth_max: u8, capacity: usize) -> BMOC {
    BMOC { depth_max, entries: Vec::with_capacity(capacity) }
  }*/

  /// We suppose here that the entries are already sorted (ASC natural ordering) with
  /// no duplicates and no small cells included into larger one's.
  pub(super) fn create_unsafe(depth_max: u8, entries: Box<[u64]>) -> BMOC {
    BMOC { depth_max, entries}
  }
  
  /*pub(super) fn add(&mut self, depth: u8, hash: u64, is_full: u8) {
    self.entries.push(build_raw_value(depth, hash, is_full, self.depth_max));
  }*/
  
  pub fn from_raw_value(&self, raw_value: u64) -> Cell {
    Cell::new(raw_value, self.depth_max)
  }
  
  /// Returns the number of cells at depth `depth_max` the moc contains, i.e.
  /// the sum for each cell of the number of cells at depth `depth_max`.
  pub fn deep_size(&self) -> usize {
    let mut sum = 0_usize;
    for &raw_value in self.entries.iter() {
      let depth = self.get_depth(raw_value);
      sum += nside_square_unsafe(self.depth_max - depth) as usize;
    }
    sum
  }
  
  /// Iterator on the BMOC raw values
  /// See method `` to extract informations from a raw value
  pub fn iter(&self) -> Iter<u64> {
    self.entries.iter()
  }
  
  /// Returns an iterator iterating over all cells at the BMOC maximum depth
  /// (the iteration is made in the natural cell order).
  pub fn flat_iter(&self) -> BMOCFlatIter {
    BMOCFlatIter::new(self.depth_max, self.deep_size(),self.entries.iter())
  }
  
  /// Returns an array containing all the BMOC cells flattened at the maximum depth.
  /// This is an utility methods bascailly calling `deep_size` to initialize an array
  /// and `flat_iter` to retrieve all cells.
  pub fn to_flat_array(&self) -> Box<[u64]> {
    let mut res: Vec<u64> = Vec::with_capacity(self.deep_size());
    for cell in self.flat_iter() {
      res.push(cell);
    }
    res.into_boxed_slice()
  }
  
  fn get_depth(&self, raw_value: u64) -> u8 {
    self.get_depth_no_flag(rm_flag(raw_value))
  }
  
  /// Works both with no flag or with flag set to 0
  fn get_depth_no_flag(&self, raw_value_no_flag: u64) -> u8 {
    self.depth_max - (raw_value_no_flag.trailing_zeros() >> 1) as u8
  }
  
}

#[inline]
fn rm_flag(raw_value: u64) -> u64 {
  raw_value >> 1
}

#[inline]
fn is_partial(raw_value: u64) -> bool {
  (raw_value & 1_u64) == 0_u64
}

#[inline]
fn is_not_first_cell_of_larger_cell(hash: u64) -> bool {
  (hash & 3_u64) != 0_u64
}

#[inline]
fn get_depth(raw_value: u64, depth_max: u8) -> u8 {
  get_depth_no_flag(rm_flag(raw_value), depth_max)
}

#[inline]
fn get_depth_no_flag(raw_value_no_flag: u64, depth_max: u8) -> u8 {
  depth_max - (raw_value_no_flag.trailing_zeros() >> 1) as u8
}

#[inline]
fn get_hash_from_delta_depth(raw_value: u64, delta_depth: u8) -> u64 {
  raw_value >> (2 + (delta_depth << 1))
}



pub struct BMOCFlatIter<'a> {
  depth_max: u8,
  deep_size: usize,
  raw_val_iter: Iter<'a, u64>,
  curr_val: Option<u64>,
  curr_val_max: u64,
  n_returned: usize,
}

impl<'a> BMOCFlatIter<'a> {
  fn new(depth_max: u8, deep_size: usize, raw_val_iter: Iter<'a, u64>) -> BMOCFlatIter<'a> {
    let mut flat_iter = BMOCFlatIter { 
      depth_max, deep_size, raw_val_iter, 
      curr_val: None, curr_val_max: 0_u64, n_returned: 0_usize
    };
    flat_iter.next_cell();
    flat_iter
  }
  
  pub fn deep_size(&self) -> usize {
    self.deep_size
  }
  
  pub fn depth(&self) -> u8 {
    self.depth_max
  }
  
  fn next_cell(&mut self) -> Option<u64> {
    match self.raw_val_iter.next() {
      None => self.curr_val.take(),
      Some(&raw_value) => {
        // Remove the flag bit, then divide by 2 (2 bits per level)
        let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
        let twice_delta_depth = delta_depth << 1;
        // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
        let hash = raw_value >> (2 + twice_delta_depth);
        let val = hash << twice_delta_depth;
        self.curr_val_max = val | ((1_u64 << twice_delta_depth) - 1_u64);
        self.curr_val.replace(val)

        /*// Remove the flag bit, then divide by 2 (2 bits per level)
        let twice_delta_depth = (raw_value >> 1).trailing_zeros() as u8;
        // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
        let mask = 0xFFFFFFFFFFFFFFFC_u64 << twice_delta_depth;
        let min = raw_value & mask;
        self.curr_val_max = min | ((!mask) >> 1);
        self.curr_val.replace(min)*/
      },
    }
  }

}

impl<'a> Iterator for BMOCFlatIter<'a> {
  type Item = u64;

  fn next(&mut self) -> Option<u64> {
    if let Some(val) = self.curr_val {
      self.n_returned += 1;
      if val < self.curr_val_max {
        self.curr_val.replace(val + 1)
      } else {
        self.next_cell() 
      }
    } else {
      None
    }
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let n = self.deep_size - self.n_returned;
    (n, Some(n))
  }
}


pub struct BMOCIter<'a> {
  depth_max: u8,
  iter: Iter<'a, u64>,
}

impl<'a> Iterator for BMOCIter<'a> {
  type Item = Cell;

  fn next(&mut self) -> Option<Cell> {
    match self.iter.next() {
      None => None,
      Some(&raw_value) => Some(Cell::new(raw_value, self.depth_max)),
    }
  }
  
  fn size_hint(&self) -> (usize, Option<usize>) {
    self.iter.size_hint()
  }
}

impl<'a> IntoIterator for &'a BMOC {
  type Item = Cell;
  type IntoIter = BMOCIter<'a>;

  fn into_iter(self) -> Self::IntoIter {
    BMOCIter { depth_max: self.depth_max, iter: self.entries.iter() }
  }
}


/// Create a BMOC raw value coding the depth, the hash and a flag in a way such that
/// the natural ordering follow a z-order curve.
///
/// # Inputs
/// - `depth`: depth of the hash value
/// - `hash`: hash value
/// - `is_full`: must be `false` (not full) or `true` (full)
/// - `depth_max`: the depth of the BMOC (we can use 29 for a unique raw value, but it will work
///   only with languages supporting unsigned 64 bit integers)
/// 
/// # Outputs
/// - the value coded like this:
///   - BBBBxx...xxS00...00F if depth < depth_max
///   - BBBBxx...xxxx...xxSF if depth = depht_max
///   - with in bith cases:
///     -  B: the 4 bits coding the base hash [0- 11]
///     - xx: the 2 bits of level x
///     -  S: the sentinel bit coding the depth
///     - 00: if (depth != depht_max) those bits are unused bits
///     -  F: the flag bit (0: partial, 1: full)
#[inline]
fn build_raw_value(depth: u8, hash: u64, is_full: bool, depth_max: u8) -> u64 {
  // Set the sentinel bit
  let mut hash = (hash << 1) | 1_u64;
  // Shift according to the depth and add space for the flag bit
  hash <<= 1 + ((depth_max - depth) << 1);
  // Set the flag bit if needed
  hash | (is_full as u64) // see https://doc.rust-lang.org/std/primitive.bool.html
}
