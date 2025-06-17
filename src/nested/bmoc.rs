//! Definition of a BMOC, i.e. a MOC storing an additional flag telling if a cell is fully
//! or partially covered by the MOC.
//!
//! So far, all BMOC logical operations (not, and, or, xor) are made from the BMOC representation.
//! It is probably simpler and faster to work on ranges (but we have to handle the flag).

use std::{
  cmp::{max, Ordering},
  fs::File,
  io::{BufRead, BufReader, Error as IoError, Write},
  path::Path,
  slice::Iter,
  str,
  time::SystemTime,
  vec::IntoIter,
};

use base64::{engine::general_purpose::STANDARD, DecodeError, Engine};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use chrono::{DateTime, SecondsFormat, Utc};
use log::debug;

use crate::{
  nested::{
    hash,
    map::fits::{
      error::FitsError,
      read::{
        check_keyword_and_parse_uint_val, check_keyword_and_str_val, check_keyword_and_val,
        get_str_val_no_quote, next_36_chunks_of_80_bytes, parse_uint_val,
      },
      write::{
        write_final_padding_ioerr, write_keyword_record, write_primary_hdu_ioerr,
        write_str_keyword_record, write_uint_mandatory_keyword_record,
      },
    },
    to_range,
  },
  nside_square_unsafe,
};

/// A very basic and simple BMOC Builder: we push elements in it assuming that we provide them
/// in the right order, without duplicates and small cells included in larger cells.
#[derive(Debug)]
pub struct BMOCBuilderUnsafe {
  // (super) // removed because of external test
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
    // println!("push {:?}", Cell::new(*v.last().unwrap(), self.depth_max));
    } else {
      panic!("Empty builder, you have to re-init it before re-using it!");
    }
    self
  }

  fn push_raw_unsafe(&mut self, raw_value: u64) -> &mut BMOCBuilderUnsafe {
    // println!("push {:?}", Cell::new(raw_value, self.depth_max));
    if let Some(ref mut v) = self.entries {
      v.push(raw_value);
    } else {
      panic!("Empty builder, you have to re-init it before re-using it!");
    }
    self
  }

  pub fn push_all(
    &mut self,
    depth: u8,
    from_hash: u64,
    to_hash: u64,
    are_full: bool,
  ) -> &mut BMOCBuilderUnsafe {
    if let Some(ref mut v) = self.entries {
      for h in from_hash..to_hash {
        v.push(build_raw_value(depth, h, are_full, self.depth_max));
      }
    } else {
      panic!("Empty builder, you have to re-init it before re-using it!");
    }
    self
  }

  #[allow(clippy::wrong_self_convention)]
  pub fn to_bmoc(mut self) -> BMOC {
    BMOC::create_unsafe(
      self.depth_max,
      self
        .entries
        .take()
        .expect("Empty builder!")
        .into_boxed_slice(),
    )
  }

  /// We consider that the pushed elements are not ordered, but they come from a valid BMOC (i.e.
  /// no cell included in another cell)
  #[allow(clippy::wrong_self_convention)]
  pub fn to_bmoc_from_unordered(mut self) -> BMOC {
    let mut res = self.entries.take().expect("Empty builder!");
    res.sort_unstable();
    BMOC::create_unsafe(self.depth_max, res.into_boxed_slice())
  }

  fn pack(&mut self) -> Vec<u64> {
    let mut entries = self.entries.take().expect("Empty builder!");
    // On-place pack
    let mut prev_to_index = 0_usize;
    let mut curr_to_index = entries.len();
    while prev_to_index != curr_to_index {
      // changes occurs
      prev_to_index = curr_to_index;
      let mut i_prev_moc = 0_usize;
      let mut i_curr_moc = 0_usize;
      while i_prev_moc < prev_to_index {
        let mut curr_cell = entries[i_prev_moc];
        i_prev_moc += 1;
        let mut curr_cell_depth = get_depth(curr_cell, self.depth_max);
        let mut curr_cell_hash =
          get_hash_from_delta_depth(curr_cell, self.depth_max - curr_cell_depth);
        // Look for the first cell of the larger cell (depth - 1)  (=> 2 last bits = 00), the cell must be FULL
        while i_prev_moc < prev_to_index
          && (curr_cell_depth == 0
            || is_partial(curr_cell)
            || is_not_first_cell_of_larger_cell(curr_cell_hash))
        {
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
          && entries[i_prev_moc]
            == build_raw_value(curr_cell_depth, curr_cell_hash | 1, true, self.depth_max)
          && entries[i_prev_moc + 1]
            == build_raw_value(curr_cell_depth, curr_cell_hash | 2, true, self.depth_max)
          && entries[i_prev_moc + 2]
            == build_raw_value(curr_cell_depth, curr_cell_hash | 3, true, self.depth_max)
        {
          entries[i_curr_moc] = build_raw_value(
            curr_cell_depth - 1,
            curr_cell_hash >> 2,
            true,
            self.depth_max,
          );
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

  fn low_depth_raw_val_at_lower_depth(&self, raw_value: u64, new_depth: u8) -> u64 {
    debug_assert!(self.get_depth(raw_value) <= new_depth);
    debug_assert!(new_depth <= self.depth_max);
    let twice_delta_depth = (self.depth_max - new_depth) << 1;
    (raw_value >> twice_delta_depth) | (raw_value & 1_u64)
  }

  // We assume the given entries form a valid BMOC (already packef, ordered, ...)
  fn to_lower_depth(&self, new_depth: u8, mut entries: Vec<u64>) -> Vec<u64> {
    if new_depth >= self.depth_max {
      panic!("The given depth must be lower than the depth max of the BMOC");
    }
    let mut i_new = 0_usize;
    let mut prev_hash_at_new_depth = loop {
      if i_new == entries.len() {
        // All cells have a depth <= new_depth
        break None;
      }
      let raw_value = entries[i_new];
      let depth = self.get_depth(raw_value);
      if depth <= new_depth {
        entries[i_new] = self.low_depth_raw_val_at_lower_depth(raw_value, new_depth);
        i_new += 1;
      } else {
        break Some(get_hash_from_delta_depth(
          raw_value,
          self.depth_max - new_depth,
        ));
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
        entries[i_new] = self.low_depth_raw_val_at_lower_depth(raw_value, new_depth);
        i_new += 1;
      } else {
        let curr_hash_at_new_depth =
          get_hash_from_delta_depth(raw_value, self.depth_max - new_depth);
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

  #[allow(clippy::wrong_self_convention)]
  pub fn to_bmoc_packing(&mut self) -> BMOC {
    let entries = self.pack();
    BMOC::create_unsafe(self.depth_max, entries.into_boxed_slice())
  }

  #[allow(clippy::wrong_self_convention)]
  pub fn to_lower_depth_bmoc(&mut self, new_depth: u8) -> BMOC {
    let entries = self.entries.take().expect("Empty builder!");
    let entries = self.to_lower_depth(new_depth, entries);
    BMOC::create_unsafe(new_depth, entries.into_boxed_slice())
  }

  #[allow(clippy::wrong_self_convention)]
  pub fn to_lower_depth_bmoc_packing(&mut self, new_depth: u8) -> BMOC {
    let entries = self.pack();
    let entries = self.to_lower_depth(new_depth, entries);
    BMOC::create_unsafe(new_depth, entries.into_boxed_slice())
  }

  #[inline]
  fn get_depth(&self, raw_value: u64) -> u8 {
    self.get_depth_no_flag(rm_flag(raw_value))
  }
  #[inline]
  /// Works both with no flag or with flag set to 0
  fn get_depth_no_flag(&self, raw_value_no_flag: u64) -> u8 {
    self.depth_max - (raw_value_no_flag.trailing_zeros() >> 1) as u8
  }
}

pub enum Status {
  /// The point is in the MOC
  IN,
  /// The point is out of the MOC
  OUT,
  /// The point may be in or out of the MOC
  UNKNOWN,
}

/// Builder taking cell at the MOC maximum depth.
pub struct BMOCBuilderFixedDepth {
  depth: u8,
  bmoc: Option<BMOC>,
  is_full: bool,
  buffer: Vec<u64>,
  sorted: bool,
}

impl BMOCBuilderFixedDepth {
  /// # Inputs
  ///  - `depth`: BMOC depth.
  ///  - `is_full`: the flag to be set for each cell number (I expect`true` to be used for example
  ///    when building catalogues MOC.
  ///
  /// The results of logical operations between BMOC having the flag of each of their cells
  /// set to `true` must equal the results of regular MOC logical operations.
  pub fn new(depth: u8, is_full: bool) -> BMOCBuilderFixedDepth {
    BMOCBuilderFixedDepth::with_capacity(depth, is_full, 10_000_000)
  }

  pub fn with_capacity(depth: u8, is_full: bool, buff_capacity: usize) -> BMOCBuilderFixedDepth {
    BMOCBuilderFixedDepth {
      depth,
      bmoc: None,
      is_full,
      buffer: Vec::with_capacity(buff_capacity),
      sorted: true,
    }
  }

  /// The hash must be at the builder depth
  pub fn push(&mut self, hash: u64) {
    if let Some(h) = self.buffer.last() {
      if *h == hash {
        return;
      } else if self.sorted && *h > hash {
        self.sorted = false;
      }
    }
    self.buffer.push(hash);
    if self.buffer.len() == self.buffer.capacity() {
      self.drain_buffer();
    }
  }

  #[allow(clippy::wrong_self_convention)]
  pub fn to_bmoc(&mut self) -> Option<BMOC> {
    // if self.buffer.len() > 0 {
    self.drain_buffer();
    // }
    self.bmoc.take()
  }

  fn drain_buffer(&mut self) {
    if !self.sorted {
      // Sort and remove duplicates
      self.buffer.sort_unstable();
      self.buffer.dedup();
    }
    let new_bmoc = self.buff_to_bmoc();
    self.clear_buff();
    self.bmoc = Some(match self.bmoc.take() {
      Some(prev_bmoc) => prev_bmoc.or(&new_bmoc),
      None => new_bmoc,
    })
  }

  fn buff_to_bmoc(&mut self) -> BMOC {
    let mut i = 0_usize;
    let mut k = 0_usize;
    while i < self.buffer.len() {
      let h = self.buffer[i];
      let sequence_len = self.largest_lower_cell_sequence_len(h, &self.buffer[i..]);
      /*{
        // Look at the maximum number of cell that could be merge if the hash is the first of a cell
        let delta_depth = (h.trailing_zeros() >> 1).min(self.depth); // low_res_cell_depth = self.depth - delta_depth
        let num_cells = 1_usize << (dd << 1); // number of depth self.depth cells in the low_res_cell = (2^dd)^2 = 2^(2*dd)
        // Look for a sequence
        let mut j = i + 1;
        let mut expected_h = h + 1_u64;

        while j < self.buffer.len() && sequence_len < num_cells && self.buffer[j] == expected_h {
          j += 1;
          sequence_len = 1;
          expected_h += 1;
        }
      }*/
      // Look at the actual low_res_cell the sequence correspond to
      let delta_depth = sequence_len.next_power_of_two();
      let delta_depth = if delta_depth > sequence_len {
        delta_depth.trailing_zeros() >> 2 // take previous value and divide by 2
      } else {
        debug_assert_eq!(delta_depth, sequence_len);
        delta_depth.trailing_zeros() >> 1 // divide by 2
      } as u8;
      let twice_dd = delta_depth << 1;
      let sequence_len = 1_usize << twice_dd;
      // Write the value
      self.buffer[k] = build_raw_value(
        self.depth - delta_depth,
        h >> twice_dd,
        self.is_full,
        self.depth,
      );
      k += 1;
      i += sequence_len;
    }
    // self.buffer.truncate(k);
    BMOC::create_unsafe_copying(self.depth, &self.buffer[0..k])
  }

  #[inline]
  fn largest_lower_cell_sequence_len(&self, mut h: u64, entries: &[u64]) -> usize {
    // Look for the maximum number of cells that could be merged if the hash is the first of a cell
    let dd = ((h.trailing_zeros() >> 1) as u8).min(self.depth); // low_res_cell_depth = self.depth - delta_depth
    let n = 1_usize << (dd << 1); // number of depth self.depth cells in the low_res_cell = (2^dd)^2 = 2^(2*dd)
                                  // Look for a sequence
    let n = n.min(entries.len());
    /*for i in 1..n {
      h += 1;
      if entries[i] != h {
        return i;
      }
    }*/
    for (i, e) in entries.iter().enumerate().take(n).skip(1) {
      h += 1;
      if *e != h {
        return i;
      }
    }
    n
  }

  fn clear_buff(&mut self) {
    self.sorted = true;
    self.buffer.clear();
  }
}

/// Structure defining a simple BMOC.
/// Three different iterators are available:
/// - `bmoc.iter() -> Iterator<u64>` : iterates on the raw value stored in the BMOC (the ordering
///   follow the z-order-curve order).
/// - `bmoc.into_iter() -> Iterator<Cell>`: same a `iter()` except that it returns Cells,
///   i.e. decoded raw value containing the `depth`, `order` and `flag`.
/// - `bmoc.flat_iter() -> Iterator<u64>`: iterates on all the cell number at the maximum depth, in
///   ascending order (flag information is lost).
/// - `bmoc.flat_iter_cell() -> Iterator<Cell>` same as `flat_iter()` but conserving then `flag`
///   information (and the depth which must always equals the BMOC depth).
#[derive(Debug, PartialEq, Eq)]
pub struct BMOC {
  depth_max: u8,
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
    let is_full = (raw_value & 1_u64) == 1_u64;
    // Remove the flag bit, then divide by 2 (2 bits per level)
    let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
    // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
    let hash = raw_value >> (2 + (delta_depth << 1));
    let depth = depth_max - delta_depth;
    Cell {
      raw_value,
      depth,
      hash,
      is_full,
    }
  }
}

impl BMOC {
  pub fn new_empty(depth: u8) -> Self {
    let builder = BMOCBuilderUnsafe::new(depth, 0);
    builder.to_bmoc()
  }

  pub fn new_allsky(depth: u8) -> Self {
    let mut builder = BMOCBuilderUnsafe::new(depth, 12);
    builder.push_all(0, 0, 12, true);
    builder.to_bmoc()
  }

  pub fn size(&self) -> usize {
    self.entries.len()
  }

  /// We suppose here that the entries are already sorted (ASC natural ordering) with
  /// no duplicates and no small cells included into larger one's.
  pub(super) fn create_unsafe(depth_max: u8, entries: Box<[u64]>) -> BMOC {
    BMOC { depth_max, entries }
  }

  pub(super) fn create_unsafe_copying(depth_max: u8, entries: &[u64]) -> BMOC {
    let mut entries_copy = Vec::with_capacity(entries.len());
    for e in entries {
      entries_copy.push(*e);
    }
    BMOC {
      depth_max,
      entries: entries_copy.into_boxed_slice(),
    }
  }

  pub fn get_depth_max(&self) -> u8 {
    self.depth_max
  }

  pub fn equals(&self, other: &BMOC) -> bool {
    if self.depth_max == other.depth_max && self.entries.len() == other.entries.len() {
      for (r1, r2) in self.iter().zip(other.iter()) {
        if r1 != r2 {
          return false;
        }
      }
      return true;
    }
    false
  }

  pub fn assert_equals(&self, other: &BMOC) {
    if self.depth_max == other.depth_max {
      for (r1, r2) in self.iter().zip(other.iter()) {
        if *r1 != *r2 {
          panic!(
            "Left: {:?}; Right: {:?}",
            self.from_raw_value(*r1),
            other.from_raw_value(*r2)
          );
        }
      }
      if self.entries.len() != other.entries.len() {
        panic!("Lengths are different");
      }
    } else {
      panic!("Depths are different");
    }
  }

  /// Test the given point and return its "Status": in, out of the MOC or maybe.
  pub fn test_coo(&self, lon: f64, lat: f64) -> Status {
    let h_raw = build_raw_value(
      self.depth_max,
      hash(self.depth_max, lon, lat),
      true,
      self.depth_max,
    );
    match self.entries.binary_search(&h_raw) {
      Ok(i) => {
        if is_partial(self.entries[i]) {
          Status::UNKNOWN
        } else {
          Status::IN
        }
      }
      Err(i) => {
        let cell = Cell::new(h_raw, self.depth_max);
        // look in next or previous cels
        if i > 0 && is_in(&self.from_raw_value(self.entries[i - 1]), &cell) {
          if is_partial(self.entries[i - 1]) {
            Status::UNKNOWN
          } else {
            Status::IN
          }
        } else if i < self.entries.len() && is_in(&self.from_raw_value(self.entries[i]), &cell) {
          if is_partial(self.entries[i]) {
            Status::UNKNOWN
          } else {
            Status::IN
          }
        } else {
          Status::OUT
        }
      }
    }
  }

  /// Returns the BMOC complement:
  /// - cells with flag set to 1 (fully covered) are removed
  /// - cells with flag set to 0 (partially covered) are kept
  /// - empty cells are added with flag set to 1
  ///
  /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
  pub fn not(&self) -> BMOC {
    // Worst case: only 1 sub-cell by cell in the MOC (+11 for depth 0)
    let mut builder = BMOCBuilderUnsafe::new(self.depth_max, 3 * self.entries.len() + 12);
    // Empty MOC, easy
    if self.entries.is_empty() {
      for h in 0..12_u64 {
        builder.push(0_u8, h, true);
      }
      return builder.to_bmoc();
    }
    // Real case
    let mut d = 0_u8;
    let mut h = 0_u64;
    // Go down to first cell
    let mut cell = self.from_raw_value(self.entries[0]);
    go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
    if !cell.is_full {
      builder.push_raw_unsafe(cell.raw_value);
    }
    // Between first and last
    for i in 1..self.entries.len() {
      cell = self.from_raw_value(self.entries[i]);
      let dd = dd_4_go_up(d, h, cell.depth, cell.hash);
      go_up(&mut d, &mut h, dd, true, &mut builder);
      go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
      if !cell.is_full {
        builder.push_raw_unsafe(cell.raw_value);
      }
    }
    // After last
    let delta_depth = d;
    go_up(&mut d, &mut h, delta_depth, true, &mut builder); // go up to depth 0
    for h in h..12 {
      // Complete with base cells if needed
      builder.push(0_u8, h, true);
    }
    builder.to_bmoc()
  }

  /*
  /// Go to the next hash value:
  /// - if the input hash is not the last one of the super-cell
  ///   (the cell of depth deph - 1 the hash belongs to), the result is simply
  ///   - output_depth = input_depth
  ///   - output_hash = input_hash + 1
  /// - else, the depth is changed (we go up) until the hash is not the last of the super-cell
  ///   and the result is:
  ///   - output_depth < input_depth
  ///   - output_hash = input_hash_at_outpu_depth + 1
  fn go_next(&self, start_depth: &mut u8, start_hash: &mut u64) {
    while *start_depth > 0 && ((*start_hash & 3_u64) == 3_u64) {
      *start_depth -= 1;
      *start_hash >>= 2;
    }
    *start_hash += 1;
  }*/

  /// Returns the intersection of this BMOC with the given BMOC:
  /// - all non overlapping cells are removed
  /// - when two cells are overlapping, the overlapping part is kept
  ///   - the value of the flag is the result of a logical AND between the flags of the merged cells.
  ///
  /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
  pub fn and(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len()),
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
          match l.hash.cmp(&hr_at_dl) {
            Ordering::Less => left = it_left.next(),
            Ordering::Greater => right = it_right.next(),
            Ordering::Equal => {
              debug_assert_eq!(l.hash, hr_at_dl);
              builder.push(r.depth, r.hash, r.is_full && l.is_full);
              right = it_right.next()
            }
          }
        }
        Ordering::Greater => {
          let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
          match hl_at_dr.cmp(&r.hash) {
            Ordering::Less => left = it_left.next(),
            Ordering::Greater => right = it_right.next(),
            Ordering::Equal => {
              debug_assert_eq!(hl_at_dr, r.hash);
              builder.push(l.depth, l.hash, r.is_full && l.is_full);
              left = it_left.next()
            }
          }
        }
        Ordering::Equal => {
          debug_assert_eq!(l.depth, r.depth);
          match l.hash.cmp(&r.hash) {
            Ordering::Less => left = it_left.next(),
            Ordering::Greater => right = it_right.next(),
            Ordering::Equal => {
              debug_assert_eq!(l.hash, r.hash);
              builder.push(l.depth, l.hash, r.is_full && l.is_full);
              left = it_left.next();
              right = it_right.next()
            }
          }
        }
      }
    }
    builder.to_bmoc()
  }

  /* Try making operations with as few if as possible, playing on indices
  fn and_v2(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len())
    );
    let mut left = self.entries;
    let mut right = other.entries;
    let mut ileft = 0_usize;
    let mut iright = 0_usize;

  }
  */

  /// Returns the union of this BMOC with the given BMOC:
  /// - all non overlapping cells in both BMOCs are kept;
  /// - overlapping cells are merged, the value of the flag is the result of a logical OR between.
  ///
  /// the flags of the merged cells.
  /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
  pub fn or(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len()),
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
          if l.hash < hr_at_dl {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if l.hash > hr_at_dl {
            builder.push(r.depth, r.hash, r.is_full);
            right = it_right.next();
          } else if l.is_full {
            debug_assert_eq!(l.hash, hr_at_dl);
            builder.push(l.depth, l.hash, l.is_full);
            right = consume_while_overlapped(l, &mut it_right);
            left = it_left.next();
          } else {
            debug_assert_eq!(l.hash, hr_at_dl);
            debug_assert!(!l.is_full);
            let mut is_overlapped = false;
            right = consume_while_overlapped_and_partial(l, &mut it_right, &mut is_overlapped);
            if is_overlapped {
              right = self.not_in_cell_4_or(l, right.unwrap(), &mut it_right, &mut builder);
            } else {
              // all flags set to 0 => put large cell with flag  = 0
              builder.push(l.depth, l.hash, false);
            }
            left = it_left.next();
          }
        }
        Ordering::Greater => {
          let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
          if hl_at_dr < r.hash {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if hl_at_dr > r.hash {
            builder.push(r.depth, r.hash, r.is_full);
            right = it_right.next();
          } else if r.is_full {
            debug_assert_eq!(hl_at_dr, r.hash);
            builder.push(r.depth, r.hash, r.is_full);
            left = consume_while_overlapped(r, &mut it_left);
            right = it_right.next();
          } else {
            debug_assert_eq!(hl_at_dr, r.hash);
            debug_assert!(!r.is_full);
            let mut is_overlapped = false;
            left = consume_while_overlapped_and_partial(r, &mut it_left, &mut is_overlapped);
            if is_overlapped {
              left = self.not_in_cell_4_or(r, left.unwrap(), &mut it_left, &mut builder);
            } else {
              // all flags set to 0 => put large cell with flag  = 0
              builder.push(r.depth, r.hash, false);
            }
            right = it_right.next();
          }
        }
        Ordering::Equal => {
          debug_assert_eq!(l.depth, r.depth);
          match l.hash.cmp(&r.hash) {
            Ordering::Less => {
              builder.push(l.depth, l.hash, l.is_full);
              left = it_left.next();
            }
            Ordering::Greater => {
              builder.push(r.depth, r.hash, r.is_full);
              right = it_right.next();
            }
            Ordering::Equal => {
              debug_assert_eq!(l.hash, r.hash);
              builder.push(l.depth, l.hash, r.is_full || l.is_full);
              left = it_left.next();
              right = it_right.next();
            }
          }
        }
      }
    }
    while let Some(l) = &left {
      debug_assert!(right.is_none());
      builder.push(l.depth, l.hash, l.is_full);
      left = it_left.next();
    }
    while let Some(r) = &right {
      debug_assert!(left.is_none());
      builder.push(r.depth, r.hash, r.is_full);
      right = it_right.next();
    }
    builder.to_bmoc_packing()
  }

  fn not_in_cell_4_or(
    &self,
    low_resolution: &Cell,
    mut c: Cell,
    iter: &mut BMOCIter,
    builder: &mut BMOCBuilderUnsafe,
  ) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    debug_assert!(c.is_full);
    go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
    builder.push(c.depth, c.hash, true);
    let mut is_overlapped = false;
    let mut cell;
    while {
      cell = consume_while_overlapped_and_partial(low_resolution, iter, &mut is_overlapped);
      is_overlapped
    } {
      c = cell.unwrap(); // if flag => right is not None
      let dd = dd_4_go_up(d, h, c.depth, c.hash);
      go_up(&mut d, &mut h, dd, false, builder);
      go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
      builder.push(c.depth, c.hash, true);
    }
    let dd = d - low_resolution.depth;
    go_up(&mut d, &mut h, dd, false, builder);
    go_down(
      &mut d,
      &mut h,
      low_resolution.depth,
      low_resolution.hash + 1,
      false,
      builder,
    );
    cell
  }

  /// Returns the symmetric difference of this BMOC with the given BMOC:
  /// - all non overlapping cells in both BMOCs are kept
  /// - when two cells are overlapping, the overlapping part is:
  ///   - removed if both flags = 1
  ///   - kept if one of the flags = 0 (since 0 meas partially covered but O don't know which part)
  ///
  /// The method as been tested when all flags are `is_full` (i.e. regular MOC case).
  pub fn xor(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len()),
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
          if l.hash < hr_at_dl {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if l.hash > hr_at_dl {
            builder.push(r.depth, r.hash, r.is_full);
            right = it_right.next();
          } else if l.is_full {
            debug_assert_eq!(l.hash, hr_at_dl);
            right = self.not_in_cell_4_xor(l, r, &mut it_right, &mut builder);
            left = it_left.next();
          } else {
            debug_assert_eq!(l.hash, hr_at_dl);
            debug_assert!(!l.is_full);
            builder.push(l.depth, l.hash, l.is_full);
            right = consume_while_overlapped(l, &mut it_right);
            left = it_left.next();
          }
        }
        Ordering::Greater => {
          let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
          if hl_at_dr < r.hash {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if hl_at_dr > r.hash {
            builder.push(r.depth, r.hash, r.is_full);
            right = it_right.next();
          } else if r.is_full {
            debug_assert_eq!(hl_at_dr, r.hash);
            left = self.not_in_cell_4_xor(r, l, &mut it_left, &mut builder);
            right = it_right.next();
          } else {
            debug_assert_eq!(hl_at_dr, r.hash);
            debug_assert!(!r.is_full);
            builder.push(r.depth, r.hash, r.is_full);
            left = consume_while_overlapped(r, &mut it_left);
            right = it_right.next();
          }
        }
        Ordering::Equal => {
          debug_assert_eq!(l.depth, r.depth);
          match l.hash.cmp(&r.hash) {
            Ordering::Less => {
              builder.push(l.depth, l.hash, l.is_full);
              left = it_left.next();
            }
            Ordering::Greater => {
              builder.push(r.depth, r.hash, r.is_full);
              right = it_right.next();
            }
            Ordering::Equal => {
              debug_assert_eq!(l.hash, r.hash);
              let both_fully_covered = r.is_full && l.is_full;
              if !both_fully_covered {
                builder.push(l.depth, l.hash, both_fully_covered);
              }
              left = it_left.next();
              right = it_right.next();
            }
          }
        }
      }
    }
    while let Some(l) = &left {
      debug_assert!(right.is_none());
      builder.push(l.depth, l.hash, l.is_full);
      left = it_left.next();
    }
    while let Some(r) = &right {
      debug_assert!(left.is_none());
      builder.push(r.depth, r.hash, r.is_full);
      right = it_right.next();
    }
    builder.to_bmoc_packing()
  }

  // add elements of the low resolution cell which are not in the c cell
  fn not_in_cell_4_xor(
    &self,
    low_resolution: &Cell,
    c: &Cell,
    iter: &mut BMOCIter,
    builder: &mut BMOCBuilderUnsafe,
  ) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
    if !c.is_full {
      builder.push(c.depth, c.hash, false);
    }
    let mut cell = iter.next();
    while let Some(c) = &cell {
      if !is_in(low_resolution, c) {
        break;
      }
      let dd = dd_4_go_up(d, h, c.depth, c.hash);
      go_up(&mut d, &mut h, dd, true, builder);
      go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
      if !c.is_full {
        builder.push(c.depth, c.hash, false);
      }
      cell = iter.next()
    }
    let dd = d - low_resolution.depth;
    go_up(&mut d, &mut h, dd, true, builder);
    go_down(
      &mut d,
      &mut h,
      low_resolution.depth,
      low_resolution.hash + 1,
      true,
      builder,
    );
    cell
  }

  /// Returns the difference of this BMOC (left) with the given BMOC (right):
  /// - all non overlapping cells of this (left) BMOC are kept
  /// - non overlapping cells of the other (right) BMOC are removed
  ///   if full, and kept if partially covered (since A MINUS B = A AND (NOT(B))
  /// - when two cells are overlapping, the overlapping part is:
  ///   - removed if both flags = 1
  ///   - kept if one of the flags = 0 (since 0 meas partially covered but O don't know which part)
  ///
  /// Poor's man implementation: A MINUS B = A AND NOT(B).
  pub fn minus(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len()),
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      match l.depth.cmp(&r.depth) {
        Ordering::Less => {
          // The l cell is larger than the r cell
          // - degrade r cell at the l cell depth
          let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
          if l.hash < hr_at_dl {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if l.hash > hr_at_dl {
            right = it_right.next();
          } else if l.is_full {
            debug_assert_eq!(l.hash, hr_at_dl);
            // add elements of the l cell which are not in common with the r cell
            right = self.not_in_cell_4_xor(l, r, &mut it_right, &mut builder);
            left = it_left.next();
          } else {
            debug_assert_eq!(l.hash, hr_at_dl);
            debug_assert!(!l.is_full);
            builder.push(l.depth, l.hash, false);
            right = consume_while_overlapped(l, &mut it_right);
            left = it_left.next();
          }
        }
        Ordering::Greater => {
          // The r cell is larger than the l cell
          // - degrade l cell at the r cell depth
          let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
          if hl_at_dr < r.hash {
            builder.push(l.depth, l.hash, l.is_full);
            left = it_left.next();
          } else if hl_at_dr > r.hash {
            // remove the r cell
            right = it_right.next();
          } else if !r.is_full || !l.is_full {
            debug_assert_eq!(hl_at_dr, r.hash);
            builder.push(l.depth, l.hash, false);
            left = it_left.next();
          }
        }
        Ordering::Equal => {
          debug_assert_eq!(l.depth, r.depth);
          match l.hash.cmp(&r.hash) {
            Ordering::Less => {
              builder.push(l.depth, l.hash, l.is_full);
              left = it_left.next();
            }
            Ordering::Greater => {
              right = it_right.next();
            }
            Ordering::Equal => {
              debug_assert_eq!(l.hash, r.hash);
              let both_fully_covered = r.is_full && l.is_full;
              if !both_fully_covered {
                builder.push(l.depth, l.hash, both_fully_covered);
              }
              left = it_left.next();
              right = it_right.next();
            }
          }
        }
      }
    }
    while let Some(l) = &left {
      debug_assert!(right.is_none());
      builder.push(l.depth, l.hash, l.is_full);
      left = it_left.next();
    }
    builder.to_bmoc_packing()
  }

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

  /// Iterator on the BMOC **raw values**. Raw values contains:
  /// * the depth if the cell;
  /// * the cell number;
  /// * a flag telling if the cell is fully covered or not.
  /// WARNING: this is probably not the method you are interested in as a user,
  /// see [flat_iter](#method.flat_iter) instead.
  /// TODO: make a public method to extract information from the raw value
  pub fn iter(&self) -> Iter<u64> {
    self.entries.iter()
  }

  /// Returns an iterator iterating over all cells at the BMOC maximum depth
  /// (the iteration is made in the natural cell order).
  pub fn flat_iter(&self) -> BMOCFlatIter {
    BMOCFlatIter::new(self.depth_max, self.deep_size(), self.entries.iter())
  }

  /// Returns an iterator iterating over all cells at the BMOC maximum depth
  /// (the iteration is made in the natural cell order).
  pub fn into_flat_iter(self) -> BMOCIntoFlatIter {
    BMOCIntoFlatIter::new(
      self.depth_max,
      self.deep_size(),
      self.entries.into_vec().into_iter(),
    )
  }

  /// Returns an iterator iterating over all cells at the BMOC maximum depth
  /// (the iteration is made in the natural cell order).  
  /// Contrary to [flat_iter](fn.flat_iter.html), the full cell information (the raw BMOC value
  /// it belongs to, its flag) is kept.
  pub fn flat_iter_cell(&self) -> BMOCFlatIterCell {
    BMOCFlatIterCell::new(self.depth_max, self.deep_size(), self.entries.iter())
  }

  /// Returns an array containing all the BMOC cells flattened at the maximum depth.
  /// This is an utility methods basically calling `deep_size` to initialize an array
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

  fn get_depth_icell(&self, raw_value: u64) -> (u8, u64) {
    // Remove the flag bit, then divide by 2 (2 bits per level)
    let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
    // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
    let hash = raw_value >> (2 + (delta_depth << 1));
    let depth = self.depth_max - delta_depth;
    (depth, hash)
  }

  /// Transform this (B)MOC as a simple (sorted) array of ranges
  /// (WARNING: the ranges are at the MOC depth, not at the depth 29).
  /// During the operation, we loose the `flag` information attached to each BMOC cell.  
  pub fn to_ranges(&self) -> Box<[std::ops::Range<u64>]> {
    let mut ranges: Vec<std::ops::Range<u64>> = Vec::with_capacity(self.entries.len());
    let mut prev_min = 0_u64;
    let mut prev_max = 0_u64;
    for cell in self.into_iter() {
      if cell.depth < self.depth_max {
        let range = to_range(cell.hash, self.depth_max - cell.depth);
        if range.start != prev_max {
          if prev_min != prev_max {
            // false only at first call, then always true
            ranges.push(prev_min..prev_max);
          }
          prev_min = range.start;
        }
        prev_max = range.end;
      } else if cell.hash == prev_max {
        prev_max += 1;
      } else {
        if prev_min != prev_max {
          // false only at first call, then always true
          ranges.push(prev_min..prev_max);
        }
        prev_min = cell.hash;
        prev_max = cell.hash + 1;
      }
    }
    if prev_min != prev_max {
      // false only at first call, then always true
      ranges.push(prev_min..prev_max);
    }
    ranges.into_boxed_slice()
  }

  /// Transform this (B)MOC as a simple (sorted) array of ranges.
  /// Ranges containing different flag values are split in sub-ranges
  pub fn to_flagged_ranges(&self) -> Vec<(std::ops::Range<u64>, bool)> {
    let mut ranges: Vec<(std::ops::Range<u64>, bool)> = Vec::with_capacity(self.entries.len());
    let mut prev_min = 0_u64;
    let mut prev_max = 0_u64;
    let mut prev_flag = self
      .entries
      .first()
      .map(|v| (v & 1_u64) == 1_u64)
      .unwrap_or(false);
    for cell in self.into_iter() {
      if cell.depth < self.depth_max {
        let range = to_range(cell.hash, self.depth_max - cell.depth);
        if range.start == prev_max && (prev_max == 0 || cell.is_full == prev_flag) {
          prev_max = range.end;
        } else {
          if prev_min != prev_max {
            // false only at first call, then always true
            ranges.push((prev_min..prev_max, prev_flag));
          }
          prev_min = range.start;
          prev_max = range.end;
          prev_flag = cell.is_full;
        }
      } else if cell.hash == prev_max && cell.is_full == prev_flag {
        prev_max += 1;
      } else {
        if prev_min != prev_max {
          // false only at first call, then always true
          ranges.push((prev_min..prev_max, prev_flag));
        }
        prev_min = cell.hash;
        prev_max = cell.hash + 1;
        prev_flag = cell.is_full;
      }
    }
    if prev_min != prev_max {
      // false only at first call, then always true
      ranges.push((prev_min..prev_max, prev_flag));
    }
    ranges.shrink_to_fit();
    ranges
  }

  /// Transform this (B)MOC in a very compressed version. We call it `lossy` because
  /// during the operation, we loose the `flag` information attached to each BMOC cell.
  /// # Remark
  /// * If needed we could store the flag information!
  ///
  /// # Info
  /// * Original idea by F.-X. Pineau (see Java library), improved by M. Reinecke (through
  ///   private communication) leading to an even better compression factor.
  /// * Although its seems (c.f. M. Reinecke) that this is quite similar to `Interpolative coding`,
  ///   M. Reinecke tests show a slightly better compression factor. M. Reinecke raised the following
  ///   question: was it worth implementing this specific case instead of using an
  ///   `Interpolative coding` library?
  ///
  /// # Idea
  /// * The basic idea consists in...
  #[allow(clippy::many_single_char_names)]
  pub fn compress_lossy(&self) -> CompressedMOC {
    let n = self.entries.len();
    let dm = self.depth_max;
    let mut b = CompressedMOCBuilder::new(dm, 4 + 3 * n);
    if n == 0 {
      // Special case of empty MOC
      if dm == 0 {
        for _ in 0..12 {
          b.push_leaf_empty();
        }
      } else {
        for _ in 0..12 {
          b.push_node_empty();
        }
      }
      return b.to_compressed_moc();
    } else if dm == 0 {
      // Special case of other MOC at depth max = 0
      let (curr_d, _) = self.get_depth_icell(self.entries[0]);
      assert_eq!(curr_d, 0);
      let mut h = 0_u64;
      for (_, curr_h) in self.entries.iter().map(|e| self.get_depth_icell(*e)) {
        for _ in h..curr_h {
          b.push_leaf_empty();
        }
        b.push_leaf_full();
        h = curr_h + 1;
      }
      for _ in h..12 {
        b.push_leaf_empty();
      }
      return b.to_compressed_moc();
    }
    // Let's start serious things
    let mut d;
    let mut h = 0;
    let (curr_d, curr_h) = self.get_depth_icell(self.entries[0]);
    // go down to curr hash
    for dd in (0..=curr_d).rev() {
      let target_h = curr_h >> (dd << 1);
      if curr_d == dm && dd == 0 {
        for _ in h..target_h {
          b.push_leaf_empty();
        }
        b.push_leaf_full();
      } else {
        for _ in h..target_h {
          b.push_node_empty();
        }
        if dd == 0 {
          b.push_node_full()
        } else {
          b.push_node_partial()
        };
      }
      h = target_h << 2;
    }
    d = curr_d;
    h = curr_h;
    // middle, go up and down
    let mut i = 1_usize;
    while i < n {
      let (curr_d, curr_h) = self.get_depth_icell(self.entries[i]);
      // go up (if needed)!
      let target_h = if d > curr_d {
        // case previous hash deeper that current hash
        let dd = d - curr_d;
        curr_h << (dd << 1)
      } else {
        // case current hash deeper that previous hash, need to go up?
        let dd = curr_d - d;
        curr_h >> (dd << 1)
      };
      let mut dd = ((63 - (h ^ target_h).leading_zeros()) >> 1) as u8;
      if dd > d {
        dd = d;
      }
      // - go up to common depth
      if dd > 0 && d == dm {
        for _ in h & 3..3 {
          // <=> (h + 1) & 3 < 4
          b.push_leaf_empty();
        }
        h >>= 2;
        dd -= 1;
        d -= 1;
      }
      for _ in 0..dd {
        for _ in h & 3..3 {
          // <=> (h + 1) & 3 < 4
          b.push_node_empty();
        }
        h >>= 2;
      }
      d -= dd;
      h += 1;
      // - go down
      let dd = curr_d - d;
      for rdd in (0..=dd).rev() {
        let target_h = curr_h >> (rdd << 1);
        if curr_d == dm && rdd == 0 {
          for _ in h..target_h {
            b.push_leaf_empty();
          }
          b.push_leaf_full();
        } else {
          for _ in h..target_h {
            b.push_node_empty();
          }
          if rdd == 0 {
            b.push_node_full()
          } else {
            b.push_node_partial()
          };
        }
        h = target_h << 2;
      }
      d = curr_d;
      h = curr_h;
      i += 1;
    }
    // - go up to depth 0
    if d == dm {
      for _ in h & 3..3 {
        b.push_leaf_empty();
      }
      h >>= 2;
      d -= 1;
    }
    for _ in 0..d {
      for _ in h & 3..3 {
        b.push_node_empty();
      }
      h >>= 2;
    }
    // - complete till base cell 11
    if dm == 0 {
      for _ in h + 1..12 {
        b.push_leaf_empty();
      }
    } else {
      for _ in h + 1..12 {
        b.push_node_empty();
      }
    }
    b.to_compressed_moc()
  }

  /// FITS serialization aiming at been the most efficient (in space and read speed) and made first
  /// of all for internal usage.
  pub fn to_fits<W: Write>(&self, mut writer: W) -> Result<(), IoError> {
    self
      .write_fits_header(&mut writer)
      .and_then(|()| self.write_fits_data(&mut writer))
      .and_then(|n_bytes_written| write_final_padding_ioerr(writer, n_bytes_written))
  }

  fn write_fits_header<W: Write>(&self, mut writer: W) -> Result<(), IoError> {
    // Prepare the header
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    it.next().unwrap()[0..30].copy_from_slice(b"SIMPLE  =                    T"); // Conform to FITS standard
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8"); // We work on bytes, i.e. 8 bits
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2"); // Number of data axis
    if self.depth_max <= 13 {
      // Number of bytes in a u32 value
      it.next().unwrap()[0..30].copy_from_slice(b"NAXIS1  =                    4");
    } else {
      // Number of bytes in a u64 value
      it.next().unwrap()[0..30].copy_from_slice(b"NAXIS1  =                    8");
    }
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.entries.len() as u64); // Len of data axis 2 = n_hash(depth)+1
    it.next().unwrap()[0..30].copy_from_slice(b"EXTEND  =                    F"); // No extension allowed
    it.next().unwrap()[0..16].copy_from_slice(b"PRODTYPE= 'BMOC'"); // Product type
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"HPXORDER", self.depth_max as u64);
    write_str_keyword_record(it.next().unwrap(), b"VAL_NAME", "FLAGGED_ZUNIQ");
    it.next().unwrap()[0..20].copy_from_slice(b"DTENDIAN= 'LITTLE  '");
    write_keyword_record(
      it.next().unwrap(),
      b"DATE    ",
      format!(
        "'{}'",
        Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true)
      )
      .as_str(),
    );
    write_keyword_record(
      it.next().unwrap(),
      b"CREATOR ",
      format!(
        "'Rust crate {} {}'",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
      )
      .as_str(),
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..])
  }
  /// Returns the number of bytes written.
  fn write_fits_data<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    let elem_byte_size = if self.depth_max <= 13 {
      for fzuniq in &self.entries {
        writer.write_all(((*fzuniq) as u32).to_le_bytes().as_ref())?;
      }
      size_of::<u32>()
    } else {
      for fzuniq in &self.entries {
        writer.write_all(fzuniq.to_le_bytes().as_ref())?;
      }
      size_of::<u64>()
    };
    Ok(self.entries.len() * elem_byte_size)
  }

  pub fn from_fits_file<P: AsRef<Path>>(path: P) -> Result<Self, FitsError> {
    File::open(path)
      .map_err(FitsError::Io)
      .and_then(|file| Self::from_fits(BufReader::new(file)))
  }

  pub fn from_fits<R: BufRead>(mut reader: R) -> Result<Self, FitsError> {
    let mut raw_header = [b' '; 2880];
    // Parse header
    let mut raw_cards_it = next_36_chunks_of_80_bytes(&mut reader, &mut raw_header)?;
    // Parse mandatory, well ordered keywords
    let (zuniq_n_bytes, n_rows) =
      check_keyword_and_val(raw_cards_it.next().unwrap(), b"SIMPLE ", b"T")
        .and_then(|()| check_keyword_and_val(raw_cards_it.next().unwrap(), b"BITPIX ", b"8"))
        .and_then(|()| check_keyword_and_val(raw_cards_it.next().unwrap(), b"NAXIS  ", b"2"))
        .and_then(|()| {
          check_keyword_and_parse_uint_val::<u64>(raw_cards_it.next().unwrap(), b"NAXIS1  ")
        })
        .and_then(|n_bytes| {
          check_keyword_and_parse_uint_val::<u64>(raw_cards_it.next().unwrap(), b"NAXIS2  ")
            .map(move |n_rows| (n_bytes, n_rows))
        })?;
    // Parse other keywords
    let mut end_found = false;
    let mut prodtype_found = false;
    let mut dtendian_found = false;
    let mut hpxoder: Option<u8> = None;
    let mut date: Option<SystemTime> = None;
    for kw_record in &mut raw_cards_it {
      match &kw_record[0..8] {
        b"EXTEND  " => check_keyword_and_val(kw_record, b"EXTEND  ", b"F"),
        b"PRODTYPE" => {
          check_keyword_and_str_val(kw_record, b"PRODTYPE", b"BMOC").map(|()| prodtype_found = true)
        }
        b"HPXORDER" => parse_uint_val::<u8>(kw_record).map(|v| hpxoder = Some(v)),
        b"DTENDIAN" => check_keyword_and_str_val(kw_record, b"DTENDIAN", b"LITTLE")
          .map(|()| dtendian_found = true),
        b"DATE    " => get_str_val_no_quote(kw_record).map(|v| {
          date = unsafe { str::from_utf8_unchecked(v) }
            .parse::<DateTime<Utc>>()
            .ok()
            .map(|dt| dt.into())
        }),
        b"CREATOR " => continue,
        b"END     " => {
          end_found = true;
          break;
        }
        _ => {
          debug!("Ignored FITS card: {}", unsafe {
            str::from_utf8_unchecked(kw_record)
          });
          continue;
        }
      }?;
    }
    // Check keywords
    if !end_found {
      return Err(FitsError::new_custom(String::from(
        "'END' keyword not found in the first 36 primary header cards.",
      )));
    }
    if !(prodtype_found & dtendian_found) {
      return Err(FitsError::new_custom(String::from(
        "One of the BMOC mandatory cards is missing in the FITS header!",
      )));
    }
    let depth = hpxoder.ok_or_else(|| {
      FitsError::new_custom(String::from(
        "'HPXORDER' keyword not found in the primary header cards.",
      ))
    })?;
    match zuniq_n_bytes {
      4 => (0..n_rows)
        .into_iter()
        .map(|_| reader.read_u32::<LittleEndian>().map(|v| v as u64))
        .collect::<Result<Vec<u64>, IoError>>()
        .map_err(FitsError::Io),
      8 => {
        // We could have mmapped the data to handle large BMOCs... let's do it another day
        // (or make an iterator and perform operation from iterators, like in MOCs...).
        (0..n_rows)
          .into_iter()
          .map(|_| reader.read_u64::<LittleEndian>())
          .collect::<Result<Vec<u64>, IoError>>()
          .map_err(FitsError::Io)
      }
      _ => Err(FitsError::new_custom(format!(
        "Wrong 'NAXIS1'. Expected: 4 or 8. Actual: {}",
        zuniq_n_bytes
      ))),
    }
    .map(|entries| BMOC::create_unsafe(depth, entries.into_boxed_slice()))
  }

  /// FITS serialization aiming at been more inter-operable than its alternative.
  /// # Remark about generalization
  /// * could be the same as a ZUNIQ MOC with a extra boolean column
  /// * MOM could be the same with either:
  ///     + any primitive type column
  ///     + a column of pointer pointing to the FITS HEAP for variable size content (e.g CSV, gzipped CSV, VOTable, ...)
  /// * Very similar to EXPLICIT skymaps by at various orders?!!
  pub fn to_bintable_fits<W: Write>(&self, mut writer: W) -> Result<(), IoError> {
    self
      .write_bintable_fits_header(&mut writer)
      .and_then(|()| self.write_bintable_fits_data(&mut writer))
      .and_then(|n_bytes_written| write_final_padding_ioerr(writer, n_bytes_written))
  }
  fn write_bintable_fits_header<W: Write>(&self, mut writer: W) -> Result<(), IoError> {
    write_primary_hdu_ioerr(&mut writer)?;
    let mut header_block = [b' '; 2880];
    let mut it = header_block.chunks_mut(80);
    // Write BINTABLE specific keywords in the buffer
    it.next().unwrap()[0..20].copy_from_slice(b"XTENSION= 'BINTABLE'");
    it.next().unwrap()[0..30].copy_from_slice(b"BITPIX  =                    8");
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS   =                    2");
    it.next().unwrap()[0..30].copy_from_slice(b"NAXIS1  =                    9"); // 8 bytes for u64 + 1 byte for boolean
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"NAXIS2  ", self.entries.len() as u64);
    it.next().unwrap()[0..30].copy_from_slice(b"PCOUNT  =                    0");
    it.next().unwrap()[0..30].copy_from_slice(b"GCOUNT  =                    1");
    it.next().unwrap()[0..30].copy_from_slice(b"TFIELDS =                    2");
    it.next().unwrap()[0..20].copy_from_slice(b"TTYPE1  = 'ZUNIQ   '");
    it.next().unwrap()[0..20].copy_from_slice(b"TFORM1  = 'K       '");
    it.next().unwrap()[0..25].copy_from_slice(b"TTYPE2  = 'FULLY_COVERED'");
    it.next().unwrap()[0..20].copy_from_slice(b"TFORM2  = 'L       '");
    it.next().unwrap()[0..16].copy_from_slice(b"PRODTYPE= 'BMOC'");
    it.next().unwrap()[0..20].copy_from_slice(b"COORDSYS= 'C       '");
    it.next().unwrap()[0..20].copy_from_slice(b"EXTNAME = 'xtension'");
    write_uint_mandatory_keyword_record(it.next().unwrap(), b"MAXORDER", self.depth_max as u64);
    write_keyword_record(
      it.next().unwrap(),
      b"DATE    ",
      format!(
        "'{}'",
        Utc::now()
          .to_rfc3339_opts(SecondsFormat::Secs, true)
          .trim_end_matches('Z')
      )
      .as_str(),
    );
    write_keyword_record(
      it.next().unwrap(),
      b"CREATOR ",
      format!(
        "'Rust crate {} {}'",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION")
      )
      .as_str(),
    );
    it.next().unwrap()[0..3].copy_from_slice(b"END");
    // Do write the header
    writer.write_all(&header_block[..])
  }
  /// Returns the number of bytes written.
  fn write_bintable_fits_data<W: Write>(&self, mut writer: W) -> Result<usize, IoError> {
    let mut n = 0_usize;
    for fzuniq in self.entries.iter().cloned() {
      let is_full = if (fzuniq as u8 & 1_u8) == 1_u8 {
        b'T'
      } else {
        b'F'
      };
      let zuniq = fzuniq >> 1;
      writer
        .write_all(zuniq.to_be_bytes().as_ref())
        .and_then(|()| writer.write_u8(is_full))?;
      n += 1;
    }
    Ok(n * (size_of::<u64>() + size_of::<u8>()))
  }

  pub fn to_csv<W: Write>(&self, mut writer: W) -> Result<(), IoError> {
    write!(
      &mut writer,
      "# BMOC order max: {}; File sorted by ZUNIQ(order, ipix).\n",
      self.depth_max
    )?;
    write!(&mut writer, "order,ipix,is_fully_covered_flag\n")?;
    for fzuniq in self.entries.iter().cloned() {
      let is_full = fzuniq as u8 & 1_u8;
      let zuniq = fzuniq >> 1; // Remove flag bit
      let n_trailing_zero = zuniq.trailing_zeros(); // must be even
      let hash = zuniq >> (1 | n_trailing_zero); // remove trailing zero plus sentinel bit
      let depth = self.depth_max - (n_trailing_zero as u8 >> 1);
      write!(&mut writer, "{},{},{}\n", depth, hash, is_full)?;
    }
    Ok(())
  }

  // to_json
  // to_ascii => pack by depth and Full/Partial ?
}

#[inline]
fn consume_while_overlapped(low_resolution: &Cell, iter: &mut BMOCIter) -> Option<Cell> {
  let mut cell = iter.next();
  while {
    match &cell {
      Some(c) => is_in(low_resolution, c),
      None => false,
    }
  } {
    cell = iter.next();
  }
  cell
}

/// Returns boolean:
/// - false = returned cell do not overlap any more
/// - true =  returned cell overlap and its flag is 'full'
#[inline]
fn consume_while_overlapped_and_partial(
  low_resolution: &Cell,
  iter: &mut BMOCIter,
  res_is_overlapped: &mut bool,
) -> Option<Cell> {
  let mut cell = iter.next();
  while {
    match &cell {
      Some(c) => {
        if is_in(low_resolution, c) {
          if c.is_full {
            *res_is_overlapped = true;
            false
          } else {
            true
          }
        } else {
          false
        }
      }
      None => false,
    }
  } {
    cell = iter.next();
  }
  cell
  /*let mut cell = iter.next();
  while {
    match &cell {
      Some(c) => is_in(low_res_depth, low_res_hash,  c.depth, c.hash),
      None => false,
    }
  } {
    if cell.is_full {
      *res_is_overlapped = true;
      return cell;
    }
    cell = iter.next();
  }
  *res_is_overlapped = false;
  cell*/
}

#[inline]
fn dd_4_go_up(d: u8, h: u64, next_d: u8, next_h: u64) -> u8 {
  // debug_assert!(d != next_d || h != next_h);
  let target_h_at_d = if next_d < d {
    // previous hash deeper than current hash => need to go up
    next_h << ((d - next_d) << 1)
  } else {
    // current hash deeper then (or equal to) previous hash => need to go up only if current hash
    next_h >> ((next_d - d) << 1)
  };
  // - look at the difference to see if we have to go up to add lower level cells
  // We look at the depth of the deeper common cell (i.e. all most significant bits are the same)
  // With XOR (^), we only set to 1 the bits which are set to 1 in a value and 0 in the other.
  // If number of leading = 64 => the two cell are identical, WRONG :/
  // If number of leading zero = 63 or 62 => are in the same cell => dd = 0
  // If number of leading zero = 60 or 61 => dd = 1
  // We just have to add .min(d) since base cells are coded on 4 bits (not 2)
  let xor = h ^ target_h_at_d;
  if xor != 0 {
    ((63_u8 - (xor.leading_zeros() as u8)) >> 1).min(d)
  } else {
    0
  }
}

/// Returns `true` if the given high resolution cell is in the low resolution cell
#[inline]
fn is_in(low_resolution: &Cell, high_resolution: &Cell) -> bool {
  low_resolution.depth <= high_resolution.depth
    && low_resolution.hash
      == (high_resolution.hash >> ((high_resolution.depth - low_resolution.depth) << 1))
}
/*
fn is_in(low_res_depth: u8, low_res_hash: u64, high_res_depth: u8, high_res_hash: u64) -> bool {
  low_res_depth < high_res_depth
    && low_res_hash == (high_res_hash >> (high_res_depth - low_res_depth) << 1)
}*/

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

pub struct BMOCIntoFlatIter {
  depth_max: u8,
  deep_size: usize,
  raw_val_iter: IntoIter<u64>,
  curr_val: Option<u64>,
  curr_val_max: u64,
  n_returned: usize,
}

impl BMOCIntoFlatIter {
  fn new(depth_max: u8, deep_size: usize, raw_val_iter: IntoIter<u64>) -> BMOCIntoFlatIter {
    let mut flat_iter = BMOCIntoFlatIter {
      depth_max,
      deep_size,
      raw_val_iter,
      curr_val: None,
      curr_val_max: 0_u64,
      n_returned: 0_usize,
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
      Some(raw_value) => {
        // Remove the flag bit, then divide by 2 (2 bits per level)
        let delta_depth = ((raw_value >> 1).trailing_zeros() >> 1) as u8;
        let twice_delta_depth = delta_depth << 1;
        // Remove 2 bits per depth difference + 1 sentinel bit + 1 flag bit
        let hash = raw_value >> (2 + twice_delta_depth);
        let val = hash << twice_delta_depth;
        self.curr_val_max = val | ((1_u64 << twice_delta_depth) - 1_u64);
        self.curr_val.replace(val)
      }
    }
  }
}

impl Iterator for BMOCIntoFlatIter {
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
      depth_max,
      deep_size,
      raw_val_iter,
      curr_val: None,
      curr_val_max: 0_u64,
      n_returned: 0_usize,
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
      }
    }
  }
}

impl Iterator for BMOCFlatIter<'_> {
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

pub struct BMOCFlatIterCell<'a> {
  depth_max: u8,
  deep_size: usize,
  raw_val_iter: Iter<'a, u64>,

  //curr_raw_val: u64,
  //curr_flag: bool,
  curr_val: Option<Cell>,
  curr_val_max: u64,

  n_returned: usize,
}

impl<'a> BMOCFlatIterCell<'a> {
  fn new(depth_max: u8, deep_size: usize, raw_val_iter: Iter<'a, u64>) -> BMOCFlatIterCell<'a> {
    let mut flat_iter = BMOCFlatIterCell {
      depth_max,
      deep_size,
      raw_val_iter,
      curr_val: None,
      curr_val_max: 0_u64,
      n_returned: 0_usize,
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

  fn next_cell(&mut self) -> Option<Cell> {
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
        self.curr_val.replace(Cell {
          raw_value,
          depth: self.depth_max,
          hash: val,
          is_full: (raw_value & 1_u64) == 1_u64,
        })
      }
    }
  }
}

impl Iterator for BMOCFlatIterCell<'_> {
  type Item = Cell;

  fn next(&mut self) -> Option<Cell> {
    if let Some(cell) = &self.curr_val {
      self.n_returned += 1;
      if cell.hash < self.curr_val_max {
        let new_cell = Cell {
          raw_value: cell.raw_value,
          depth: self.depth_max,
          hash: cell.hash + 1,
          is_full: cell.is_full,
        };
        self.curr_val.replace(new_cell)
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

impl Iterator for BMOCIter<'_> {
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
    BMOCIter {
      depth_max: self.depth_max,
      iter: self.entries.iter(),
    }
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

/// Fill with all cells from `start_hash` at `start_depth` to `start_hash_at_target_depth + 1`.
/// with `target_depth` = `start_depth - delta_depth`.
/// - `flag`: value of the is_full flag to be set in cells while going up
///
/// The output depth is the input depth minus delta_depth
/// The output hash value is the input hash at the output depth, plus one
fn go_up(
  start_depth: &mut u8,
  start_hash: &mut u64,
  delta_depth: u8,
  flag: bool,
  builder: &mut BMOCBuilderUnsafe,
) {
  // let output_depth = *start_depth - delta_depth;       // For debug only
  // let output_hash = (*start_hash >> (delta_depth << 1)) + 1; // For debug only
  for _ in 0_u8..delta_depth {
    let target_hash = *start_hash | 3_u64;
    for h in (*start_hash + 1)..=target_hash {
      builder.push(*start_depth, h, flag);
    }
    *start_hash >>= 2;
    *start_depth -= 1;
  }
  *start_hash += 1;
  // debug_assert_eq!(*start_depth, output_depth);
  // debug_assert_eq!(*start_hash, output_hash);
}

fn go_down(
  start_depth: &mut u8,
  start_hash: &mut u64,
  target_depth: u8,
  target_hash: u64,
  flag: bool,
  builder: &mut BMOCBuilderUnsafe,
) {
  debug_assert!(target_depth >= *start_depth);
  let mut twice_dd = (target_depth - *start_depth) << 1;
  for d in *start_depth..=target_depth {
    //range(0, target_depth - start_depth).rev() {
    let target_h_at_d = target_hash >> twice_dd;
    for h in *start_hash..target_h_at_d {
      builder.push(d, h, flag);
    }
    if d != target_depth {
      *start_hash = target_h_at_d << 2;
      twice_dd -= 2;
    }
  }
  *start_depth = target_depth;
  *start_hash = target_hash;
}

pub struct CompressedMOCBuilder {
  moc: Vec<u8>,
  depth_max: u8,
  ibyte: usize,
  ibit: u8,
}

impl CompressedMOCBuilder {
  /// Capacity = number of bytes.
  fn new(depth_max: u8, capacity: usize) -> CompressedMOCBuilder {
    let mut moc = vec![0_u8; capacity + 1];
    moc[0] = depth_max;
    CompressedMOCBuilder {
      moc,
      depth_max,
      ibyte: 1,
      ibit: 0,
    }
  }

  #[allow(clippy::wrong_self_convention)]
  fn to_compressed_moc(mut self) -> CompressedMOC {
    self.moc.resize(
      if self.ibit == 0 {
        self.ibyte
      } else {
        self.ibyte + 1
      },
      0,
    );
    CompressedMOC {
      moc: self.moc.into_boxed_slice(),
      depth_max: self.depth_max,
    }
  }

  fn push_0(&mut self) {
    self.ibyte += (self.ibit == 7) as usize;
    self.ibit += 1;
    self.ibit &= 7;
  }
  fn push_1(&mut self) {
    self.moc[self.ibyte] |= 1_u8 << self.ibit;
    self.push_0();
  }
  fn push_node_empty(&mut self) {
    self.push_1();
    self.push_0();
  }
  fn push_node_full(&mut self) {
    self.push_1();
    self.push_1();
  }
  fn push_node_partial(&mut self) {
    self.push_0();
  }
  fn push_leaf_empty(&mut self) {
    self.push_0();
  }
  fn push_leaf_full(&mut self) {
    self.push_1();
  }
}

pub struct CompressedMOCDecompHelper<'a> {
  moc: &'a [u8],
  ibyte: usize,
  ibit: u8,
}

impl<'a> CompressedMOCDecompHelper<'a> {
  fn new(moc: &'a [u8]) -> CompressedMOCDecompHelper<'a> {
    CompressedMOCDecompHelper {
      moc,
      ibyte: 1,
      ibit: 0,
    }
  }

  fn get(&mut self) -> bool {
    let r = self.moc[self.ibyte] & (1_u8 << self.ibit) != 0;
    self.ibyte += (self.ibit == 7) as usize;
    self.ibit += 1;
    self.ibit &= 7;
    r
  }
}

/// First elements contains the maximum depth
pub struct CompressedMOC {
  moc: Box<[u8]>,
  depth_max: u8,
}

impl CompressedMOC {
  pub fn depth(&self) -> u8 {
    self.depth_max
  }

  pub fn byte_size(&self) -> usize {
    self.moc.len()
  }

  pub fn to_b64(&self) -> String {
    STANDARD.encode(&self.moc)
  }

  pub fn from_b64(b64_encoded: String) -> Result<CompressedMOC, DecodeError> {
    let decoded = STANDARD.decode(b64_encoded)?;
    let depth_max = decoded[0];
    Ok(CompressedMOC {
      moc: decoded.into_boxed_slice(),
      depth_max,
    })
  }

  pub fn self_decompress(&self) -> BMOC {
    CompressedMOC::decompress(&self.moc)
  }

  // TODO: create an iterator (to iterate on cells while decompressing)
  pub fn decompress(cmoc: &[u8]) -> BMOC {
    let depth_max = cmoc[0];
    let mut moc_builder = BMOCBuilderUnsafe::new(depth_max, 8 * (cmoc.len() - 1));
    let mut bits = CompressedMOCDecompHelper::new(cmoc);
    let mut depth = 0_u8;
    let mut hash = 0_u64;
    while depth != 0 || hash != 12 {
      if bits.get() {
        // bit = 1
        if depth == depth_max || bits.get() {
          moc_builder.push(depth, hash, true);
        }
        // go up if needed
        while hash & 3 == 3 && depth > 0 {
          hash >>= 2;
          depth -= 1;
        }
        // take next hash
        hash += 1;
      } else {
        // bit = 0
        if depth == depth_max {
          // go up if needed
          while hash & 3 == 3 && depth > 0 {
            hash >>= 2;
            depth -= 1;
          }
          // take next hash
          hash += 1;
        } else {
          debug_assert!(depth < depth_max);
          // go down of 1 level
          hash <<= 2;
          depth += 1;
        }
      }
    }
    moc_builder.to_bmoc()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::nested::{cone_coverage_approx_custom, polygon_coverage};

  fn build_compressed_moc_empty(depth: u8) -> CompressedMOC {
    let mut builder = BMOCBuilderFixedDepth::new(depth, true);
    builder.to_bmoc().unwrap().compress_lossy()
  }

  fn build_compressed_moc_full(depth: u8) -> CompressedMOC {
    let mut builder = BMOCBuilderFixedDepth::new(depth, true);
    for icell in 0..12 * (1 << (depth << 1)) {
      builder.push(icell)
    }
    let bmoc = builder.to_bmoc().unwrap();
    eprintln!("Entries: {}", bmoc.entries.len());
    bmoc.compress_lossy()
  }

  fn to_aladin_moc(bmoc: &BMOC) {
    print!("draw moc {}/", bmoc.get_depth_max());
    for cell in bmoc.flat_iter() {
      print!("{}, ", cell);
    }
  }

  #[test]
  fn testok_compressed_moc_empty_d0() {
    let compressed = build_compressed_moc_empty(0);
    assert_eq!(compressed.byte_size(), 1 + 2);
    assert_eq!(compressed.moc, vec![0_u8, 0_u8, 0_u8].into_boxed_slice());
    let b64 = compressed.to_b64();
    assert_eq!(b64, "AAAA");
    assert_eq!(
      CompressedMOC::decompress(&compressed.moc)
        .compress_lossy()
        .to_b64(),
      b64
    );
  }

  #[test]
  fn testok_compressed_moc_empty_d1() {
    let compressed = build_compressed_moc_empty(1);
    assert_eq!(compressed.byte_size(), 1 + 24 / 8);
    assert_eq!(
      compressed.moc,
      vec![1_u8, 85_u8, 85_u8, 85_u8].into_boxed_slice()
    );
    let b64 = compressed.to_b64();
    assert_eq!(b64, "AVVVVQ==");
    assert_eq!(
      CompressedMOC::decompress(&compressed.moc)
        .compress_lossy()
        .to_b64(),
      b64
    );
  }

  #[test]
  fn testok_compressed_moc_full_d0() {
    let compressed = build_compressed_moc_full(0);
    assert_eq!(compressed.byte_size(), 1 + 2);
    assert_eq!(compressed.moc, vec![0_u8, 255_u8, 15_u8].into_boxed_slice());
    // eprintln!("{}", compressed.to_b64());
    let b64 = compressed.to_b64();
    assert_eq!(b64, "AP8P");
    assert_eq!(
      CompressedMOC::decompress(&compressed.moc)
        .compress_lossy()
        .to_b64(),
      b64
    );
  }

  #[test]
  fn testok_compressed_moc_full_d1() {
    let compressed = build_compressed_moc_full(1);
    assert_eq!(compressed.byte_size(), 1 + 24 / 8);
    eprintln!("{:?}", compressed.moc);
    eprintln!("{}", compressed.to_b64());
    let b64 = compressed.to_b64();
    assert_eq!(b64, "Af///w==");
    assert_eq!(
      CompressedMOC::decompress(&compressed.moc)
        .compress_lossy()
        .to_b64(),
      b64
    );
  }

  #[test]
  fn test_ok_allsky_and_empty_bmoc() {
    let bmoc_allsky = BMOC::new_allsky(18);
    assert_eq!(bmoc_allsky.entries.len(), 12);
    let bmoc_empty = BMOC::new_empty(18);
    assert_eq!(bmoc_empty.entries.len(), 0);
    assert_eq!(bmoc_allsky.not(), bmoc_empty);
    assert_eq!(bmoc_allsky, bmoc_empty.not());
  }

  #[test]
  fn test_ok_bmoc_not_round_trip() {
    let mut poly_vertices_deg: [f64; 14] = [
      272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
      -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let vertices: Vec<(f64, f64)> = poly_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let moc_org = polygon_coverage(14, vertices.as_slice(), true);
    let moc_not = moc_org.not();
    let moc_out = moc_not.not();
    println!("len: {}", moc_org.size());
    println!("len: {}", moc_not.size());
    println!("len: {}", moc_out.size());
    assert_eq!(moc_org, moc_out);
  }

  #[test]
  fn test_ok_bmoc_union_and_not() {
    let mut poly1_vertices_deg: [f64; 14] = [
      272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
      -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly1_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let mut poly2_vertices_deg: [f64; 8] = [
      272.630446, -19.234210, 272.637274, -19.248542, 272.638942, -19.231476, 272.630868,
      -19.226364,
    ];
    for v in poly2_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let v1: Vec<(f64, f64)> = poly1_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly1_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let v2: Vec<(f64, f64)> = poly2_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly2_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let moc1 = polygon_coverage(14, v1.as_slice(), true);
    let moc2 = polygon_coverage(14, v2.as_slice(), true);
    let not_moc1 = moc1.not();
    let not_moc2 = moc2.not();
    let union = moc1.or(&moc2);
    let not_inter = not_moc1.and(&not_moc2);
    let union_bis = not_inter.not();
    //to_aladin_moc(&moc1);
    //println!("\n");
    //to_aladin_moc(&moc2);
    //println!("\n");
    // to_aladin_moc(&union);
    //println!("\n");
    // to_aladin_moc(&union_bis);
    //println!("\n");
    assert_eq!(union, union_bis);
  }

  #[test]
  fn test_ok_bmoc_xor_minus_coherency() {
    // No overlapping parts, so we do no test thoroughly XOR and MINUS!!
    let mut poly1_vertices_deg: [f64; 14] = [
      272.511081, -19.487278, 272.515300, -19.486595, 272.517029, -19.471442, 272.511714,
      -19.458837, 272.506430, -19.459001, 272.496401, -19.474322, 272.504821, -19.484924,
    ];
    for v in poly1_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let mut poly2_vertices_deg: [f64; 8] = [
      272.630446, -19.234210, 272.637274, -19.248542, 272.638942, -19.231476, 272.630868,
      -19.226364,
    ];
    for v in poly2_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let mut poly3_vertices_deg: [f64; 178] = [
      272.536719, -19.461249, 272.542612, -19.476380, 272.537389, -19.491509, 272.540192,
      -19.499823, 272.535455, -19.505218, 272.528024, -19.505216, 272.523437, -19.500298,
      272.514082, -19.503376, 272.502271, -19.500966, 272.488647, -19.490390, 272.481932,
      -19.490913, 272.476737, -19.486589, 272.487633, -19.455645, 272.500386, -19.444996,
      272.503003, -19.437557, 272.512303, -19.432436, 272.514132, -19.423973, 272.522103,
      -19.421523, 272.524511, -19.413250, 272.541021, -19.400024, 272.566264, -19.397500,
      272.564202, -19.389111, 272.569055, -19.383210, 272.588186, -19.386539, 272.593376,
      -19.381832, 272.596327, -19.370541, 272.624911, -19.358915, 272.629256, -19.347842,
      272.642277, -19.341020, 272.651322, -19.330424, 272.653174, -19.325079, 272.648903,
      -19.313708, 272.639616, -19.311098, 272.638128, -19.303083, 272.632705, -19.299839,
      272.627971, -19.289408, 272.628226, -19.276293, 272.633750, -19.270590, 272.615109,
      -19.241810, 272.614704, -19.221196, 272.618224, -19.215572, 272.630809, -19.209945,
      272.633540, -19.198681, 272.640711, -19.195292, 272.643028, -19.186751, 272.651477,
      -19.182729, 272.649821, -19.174859, 272.656782, -19.169272, 272.658933, -19.161883,
      272.678012, -19.159481, 272.689173, -19.176982, 272.689395, -19.183512, 272.678006,
      -19.204016, 272.671112, -19.206598, 272.664854, -19.203523, 272.662760, -19.211156,
      272.654435, -19.214434, 272.652969, -19.222085, 272.656724, -19.242136, 272.650071,
      -19.265092, 272.652868, -19.274296, 272.660871, -19.249462, 272.670041, -19.247807,
      272.675533, -19.254935, 272.673291, -19.273917, 272.668710, -19.279245, 272.671460,
      -19.287043, 272.667507, -19.293933, 272.669261, -19.300601, 272.663969, -19.307130,
      272.672626, -19.308954, 272.675225, -19.316490, 272.657188, -19.349105, 272.657638,
      -19.367455, 272.662447, -19.372035, 272.662232, -19.378566, 272.652479, -19.386871,
      272.645819, -19.387933, 272.642279, -19.398277, 272.629282, -19.402739, 272.621487,
      -19.398197, 272.611782, -19.405716, 272.603367, -19.404667, 272.586162, -19.422703,
      272.561792, -19.420008, 272.555815, -19.413012, 272.546500, -19.415611, 272.537427,
      -19.424213, 272.533081, -19.441402,
    ];
    for v in poly3_vertices_deg.iter_mut() {
      *v = (*v).to_radians();
    }
    let v1: Vec<(f64, f64)> = poly1_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly1_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let v2: Vec<(f64, f64)> = poly2_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly2_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let v3: Vec<(f64, f64)> = poly3_vertices_deg
      .iter()
      .copied()
      .step_by(2)
      .zip(poly3_vertices_deg.iter().copied().skip(1).step_by(2))
      .collect();
    let moc1 = polygon_coverage(18, v1.as_slice(), true);
    let moc2 = polygon_coverage(18, v2.as_slice(), true);
    let moc3 = polygon_coverage(18, v3.as_slice(), true);

    let union12 = moc1.or(&moc2);
    let res1 = moc3.xor(&union12);
    let res2 = moc3.minus(&union12);
    let res3 = moc3.and(&union12.not());

    assert_eq!(res1, res2);
    assert_eq!(res1, res3);
  }

  #[test]
  fn test_ok_bmoc_not() {
    let lon = 13.158329_f64.to_radians();
    let lat = -72.80028_f64.to_radians();
    let radius = 5.64323_f64.to_radians();
    let depth = 6;
    let delta_depth = 5;

    let bmoc = cone_coverage_approx_custom(depth, delta_depth, lon, lat, radius);

    let not_bmoc = bmoc.not();

    println!("BMOC");
    for (range, flag) in bmoc.to_flagged_ranges() {
      println!("flag: {}; range: {:?}", flag, range);
    }
    println!("NOT BMOC");
    for (range, flag) in not_bmoc.to_flagged_ranges() {
      println!("flag: {}; range: {:?}", flag, range);
    }
    // Asssert that the first range has the flag 'true'.
    assert!(not_bmoc.to_flagged_ranges().first().unwrap().1)
  }
}
