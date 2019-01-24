//! Definition of a BMOC, i.e. a MOC storing an additional flag telling if a cell is fully
//! or partially covered by the MOC.

use std::slice::Iter;
use std::cmp::max;
// use std::opt::Range;

use super::{to_range};
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
    let is_full = (raw_value & 1_u64) == 1_u64;
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
  
  pub fn equals(&self, other: &BMOC) -> bool {
    if self.depth_max == other.depth_max && self.entries.len() == other.entries.len() {
      for (r1, r2) in self.iter().zip(other.iter()) {
        if r1 != r2 {
          return false;
        }
      }
      return true;
    }
    return false;
  }
  
  /// Returns the BMOC complement:
  /// - cells with flag set to 1 (fully covered) are removed
  /// - cells with flag set to 0 (partially covered) are kept
  /// - empty cells are added with flag set to 1
  /// TODO: test this method
  pub fn not(&self) -> BMOC {
    // Worst case: only 1 sub-cell by cell in the MOC (+11 for depth 0)
    let mut builder = BMOCBuilderUnsafe::new(self.depth_max, 3 * self.entries.len() + 12);
    // Empty MOC, easy
    if self.entries.len() == 0 {
      for h in 0..12_u64 {
        builder.push(0_u8, h, true);
      }
      return  builder.to_bmoc();
    }
    // Real case
    let mut d = 0_u8;
    let mut h = 0_u64;
    // Go down to first cell
    let mut cell = self.from_raw_value(self.entries[0]);
// println!("Go down. From: d: {}, h: {}. To d: {}, h: {}", d, h, cell.depth, cell.hash);
    self.go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
// println!("done. d: {}, h: {}", &d, &h);

    if !cell.is_full {
      builder.push_raw_unsafe(cell.raw_value);
    }
    // self.go_next(&mut d, &mut h);
    // Between first and last
    for i in 1..self.entries.len() {
      cell = self.from_raw_value(self.entries[i]);
      // go up (if needed)!
      // - prepare variable
      /*let target_h_at_d = if cell.depth < d { // previous hash deeper that current hash
        cell.hash << ((d - cell.depth) << 1)
      } else {                                //  current hash deeper that previous hash, need to go up?
        cell.hash >> ((cell.depth - d) << 1)
      };
      // - look at the difference to see if we have to go up to add lower level cells
      let dd = h ^ target_h_at_d;
      let dd = if dd != 0_u64 { // min to handle case diff is in base cell
        ((63_u8 - (dd.leading_zeros() as u8)) >> 1).min(d)
      } else {
        0_u8
      };*/
// println!("DD. From: d: {}, h: {}. To d: {}, h: {}", d, h, cell.depth, cell.hash);

      let dd = dd_4_go_up(d, h, cell.depth, cell.hash);
      // - go up to common depth
// println!("Go up. From: d: {}, h: {}. To d: {}. To d: {}, h: {}", d, h, d - dd, cell.depth, cell.hash);
      self.go_up(&mut d, &mut h, dd, true, &mut builder);
// println!("done. d: {}, h: {}", &d, &h);
// println!("Go down. From: d: {}, h: {}. To d: {}, h: {}", d, h, cell.depth, cell.hash);
      // go down!
      self.go_down(&mut d, &mut h, cell.depth, cell.hash, true, &mut builder);
      if !cell.is_full {
        builder.push_raw_unsafe(cell.raw_value);
      }
      // self.go_next(&mut d, &mut h);
// println!("done. d: {}, h: {}", &d, &h);// Put right element with flag = 1
    }
    // After last
    let delta_depth = d;
// println!("Go up. From: d: {}, h: {}. To d: {}. To d: {}, h: {}", d, h, d - delta_depth, cell.depth, cell.hash);
    self.go_up(&mut d, &mut h, delta_depth, true, &mut builder); // go up to depth 0
// println!("done. d: {}, h: {}", &d, &h);
// println!("complete...");
    for h in h..12 { // Complete with base cells if needed
      builder.push(0_u8, h, true);
    }
// println!("done. d: {}, h: {}", &d, &h);
    builder.to_bmoc()
  }


  /// Starts at a given cell at a given depth, goes down to a deeper cell (target), adding all 
  /// non-overlapping cells in the builder.
  /// In output, the given cell in input points to the target cell.
  /// - `flag`: value of the is_full flag to be set in cells while going down
  /*fn go_down_first(&self, start_depth: &mut u8, start_hash: &mut u64,
             target_depth: u8, target_hash: u64, flag: bool, builder: &mut BMOCBuilderUnsafe) {
    debug_assert!(target_depth >= *start_depth);
    let mut twice_dd = (target_depth - *start_depth) << 1;
    for d in *start_depth..=target_depth { //range(0, target_depth - start_depth).rev() {
      let target_h_at_d = target_hash >> twice_dd;
      // println!(" - target_h_at_d: {}", target_h_at_d);
      for h in *start_hash..target_h_at_d {
        builder.push(d, h, flag);
      }
      if d != target_depth {
        *start_hash = target_h_at_d << 2;
        twice_dd -= 2;
      }
    }
    *start_depth = target_depth;
    *start_hash  = target_hash;
  }*/
  
  /// Starts at a given cell at a given depth, goes down to a deeper cell (target), adding all 
  /// non-overlapping cells (but the starting cell) in the builder.
  /// In output, the given cell in input points to the target cell.
  /// - `flag`: value of the is_full flag to be set in cells while going down
  fn go_down(&self, start_depth: &mut u8, start_hash: &mut u64, 
             target_depth: u8, target_hash: u64, flag: bool, builder: &mut BMOCBuilderUnsafe) {
    debug_assert!(target_depth >= *start_depth);
 // self.go_next(start_depth, start_hash);
    let mut twice_dd = (target_depth - *start_depth) << 1;
    for d in *start_depth..=target_depth { //range(0, target_depth - start_depth).rev() {
      let target_h_at_d = target_hash >> twice_dd;
// println!(" - target_h_at_d: {}", target_h_at_d);
      for h in *start_hash..target_h_at_d {
        builder.push(d, h, flag);
      }
      if d != target_depth {
        *start_hash = target_h_at_d << 2;
        twice_dd -= 2;
      }
    }
    *start_depth = target_depth;
    *start_hash  = target_hash;
    self.go_next(start_depth, start_hash);
  }
  
  fn go_next(&self, start_depth: &mut u8, start_hash: &mut u64) {
    while *start_depth > 0 && ((*start_hash & 3_u64) == 3_u64) {
      *start_depth -= 1;
      *start_hash >>= 2;
    }
    *start_hash += 1;
  }
  
  
  
  /// Fill with all cells from `start_hash` at  `start_depth` to `start_hash_at_target_depth + 1`.
  /// with `target_depth` = `start_depth - delta_depth`.
  /// - `flag`: value of the is_full flag to be set in cells while going up
  fn go_up(&self, start_depth: &mut u8, start_hash: &mut u64, delta_depth: u8, flag: bool, builder: &mut BMOCBuilderUnsafe) {
    for _ in 0_u8..delta_depth {
      
      let target_hash = (*start_hash | 3_u64) + 1_u64;
      // let depth = *start_depth - dd;
      for h in *start_hash..target_hash {
        builder.push(*start_depth, h, flag);
      }
      *start_hash = target_hash >> 2;
      *start_depth -= 1;
      
      /*let target_hash = *star348t_hash | 3_u64;
      // let depth = *start_depth - dd;
      for h in (*start_hash + 1)..=target_hash {
        builder.push(*start_depth, h, flag);
      }
      *start_hash = target_hash >> 2;
      *start_depth -= 1;*/
    }
  }
  
  /// Returns the intersection of this BMOC with the given BMOC:
  /// - all non overlapping cells are removed
  /// - when two cells are overlapping, the overlapping part is kept
  ///   - the value of the flag is the result of a logical AND between the flags of the merged cells.
  /// TODO: test the method
  pub fn and(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max), 
      max(self.entries.len(), other.entries.len())
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      if l.depth < r.depth {
        let hr_at_dl = r.hash >> ((r.depth - l.depth) << 1);
        if l.hash < hr_at_dl {
          left = it_left.next();
        } else if l.hash > hr_at_dl {
          right = it_right.next();
        } else {
          debug_assert_eq!(l.hash, hr_at_dl);
          builder.push(r.depth, r.hash, r.is_full && l.is_full);
          right = it_right.next();
        }
      } else if l.depth > r.depth {
        let hl_at_dr = l.hash >> ((l.depth - r.depth) << 1);
        if hl_at_dr < r.hash {
          left = it_left.next();
        } else if hl_at_dr > r.hash {
          right = it_right.next();
        } else {
          debug_assert_eq!(hl_at_dr, r.hash);
          builder.push(l.depth, l.hash, r.is_full && l.is_full);
          left = it_left.next();
        }
      } else {
        debug_assert_eq!(l.depth, r.depth);
        if l.hash < r.hash {
          left = it_left.next();
        } else if l.hash > r.hash  {
          right = it_right.next();
        } else {
          debug_assert_eq!(l.hash, r.hash);
          builder.push(l.depth, l.hash, r.is_full && l.is_full);
          left = it_left.next();
          right = it_right.next();
        }
      }
    }
    builder.to_bmoc()
  }

  /// Returns the union of this BMOC with the given BMOC:
  /// - all non overlapping cells in both BMOCs are kept
  /// - overlapping cells are merged, the value of the flag is the result of a logical OR between 
  /// the flags of the merged cells.
  /// TODO: test the method
  pub fn or(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len())
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      if l.depth < r.depth {
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
          } else { // all flags set to 0 => put large cell with flag  = 0
            builder.push(l.depth, l.hash, false);
          }
          left = it_left.next();
        }
      } else if l.depth > r.depth {
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
          } else { // all flags set to 0 => put large cell with flag  = 0
            builder.push(r.depth, r.hash, false);
          }
          right = it_right.next();
        }
      } else {
        debug_assert_eq!(l.depth, r.depth);
        if l.hash < r.hash {
          builder.push(l.depth, l.hash, l.is_full);
          left = it_left.next();
        } else if l.hash > r.hash  {
          builder.push(r.depth, r.hash, r.is_full);
          right = it_right.next();
        } else {
          debug_assert_eq!(l.hash, r.hash);
          builder.push(l.depth, l.hash, r.is_full || l.is_full);
          left = it_left.next();
          right = it_right.next();
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

  fn not_in_cell_4_or(&self, low_resolution: &Cell, mut c: Cell,  iter: &mut BMOCIter, builder: &mut BMOCBuilderUnsafe) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    debug_assert_eq!(true, c.is_full);
    self.go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
    builder.push(c.depth, c.hash, true);
    let mut is_overlapped = false;
    let mut cell = None;
    while {
      cell = consume_while_overlapped_and_partial(low_resolution, iter, &mut is_overlapped);
      is_overlapped
    } {
      c = cell.unwrap(); // if flag => right is not None
      let dd = dd_4_go_up(d, h, c.depth, c.hash);
      self.go_up(&mut d, &mut h, dd, false, builder);
      self.go_down(&mut d, &mut h, c.depth, c.hash, false, builder);
      builder.push(c.depth, c.hash, true);
    }
    let dd = d - low_resolution.depth;
    self.go_up(&mut d, &mut h, dd, false, builder);
    self.go_down(&mut d, &mut h, low_resolution.depth, low_resolution.hash + 1, false, builder);
    cell
  }
  
  
  /// Returns the symmetric difference of this BMOC with the given BMOC:
  /// - all non overlapping cells in both BMOCs are kept
  /// - when two cells are overlapping, the overlapping part is:
  ///   - removed if both flags = 1
  ///   - kept if one of the flags = 0 (since 0 meas partially covered but O don't know which part)
  /// TODO: test the method
  pub fn xor(&self, other: &BMOC) -> BMOC {
    let mut builder = BMOCBuilderUnsafe::new(
      max(self.depth_max, other.depth_max),
      max(self.entries.len(), other.entries.len())
    );
    let mut it_left = self.into_iter();
    let mut it_right = other.into_iter();
    let mut left = it_left.next();
    let mut right = it_right.next();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    while let (Some(l), Some(r)) = (&left, &right) {
      if l.depth < r.depth {
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
      } else if l.depth > r.depth {
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
      } else {
        debug_assert_eq!(l.depth, r.depth);
        if l.hash < r.hash {
          builder.push(l.depth, l.hash, l.is_full);
          left = it_left.next();
        } else if l.hash > r.hash  {
          builder.push(r.depth, r.hash, r.is_full);
          right = it_right.next();
        } else {
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

  fn not_in_cell_4_xor(&self, low_resolution: &Cell, mut c: &Cell,  iter: &mut BMOCIter, builder: &mut BMOCBuilderUnsafe) -> Option<Cell> {
    let mut d = low_resolution.depth;
    let mut h = low_resolution.hash;
    self.go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
    if !c.is_full {
      builder.push(c.depth, c.hash, false);
    }
    let mut cell = iter.next();
    while let Some(c) = &cell {
      /*match &cell {
        Some(cc) => {
          c = cc;
          is_in(low_res_depth, low_res_hash,  c.depth, c.hash)
        },
        None => false,
      }
    } {*/
      if !is_in(low_resolution,  c) {
        break;
      }
      let dd = dd_4_go_up(d, h, c.depth, c.hash);
      self.go_up(&mut d, &mut h, dd, true, builder);
      self.go_down(&mut d, &mut h, c.depth, c.hash, true, builder);
      if !c.is_full {
        builder.push(c.depth, c.hash, false);
      }
      cell = iter.next()
    }
    let dd = d - low_resolution.depth;
    self.go_up(&mut d, &mut h, dd, true, builder);
    self.go_down(&mut d, &mut h, low_resolution.depth, low_resolution.hash + 1, true, builder);
    cell
  }
  
  /* KEEP THIS METHOD FOR A MOC
  pub fn from_ranges_unsafe(depth_max: u8, ranges: &[range]) {
    final long dMask = (1L << hhCoder.nBits(d, d)) - 1; // 11 for HEALPix
    final long dNHash = hhCoder.nHash(d, d); // 4 for HEALLPix
    final SortedHHashSet res = new SortedHHashSet(hhCoder);
    // We recall that hM is exclusive!
    // We ignore successive ranges (range_i_max == range_i+1_min)
    for range in ranges {
      width = range.to - range.from;
      // Deal with level 0 range
      if ranges.len() == 1 && width == n_hash_unsafe(depth_max) {
        assert_eq!(range.from, 0_u64);
        for (hm = 0; hm < hhCoder.nHash(0); hm++) {
          res.add(new HHImpl(0, hm));
        }
        hm = hM;
      }
            // ++ pre part
            // While we don't start with the 2 last bits = 00. 
            while ((((hm & dMask) != 0 // Hash don't end by a series of 0
                    || (hM - hm) < dNHash))// Not enough distinct hash to fill the depth 
                    //|| (d == 0 && hM - hm == 12)) // Special case of d=0 and range=0-12
                    && hm < hM) { // Still elements in the range
                final HHash hh = new HHImpl(d, hm);
if (debug) System.out.println("Add1 " + hh + " --> " + hhCoder.hashRange(hh.depth(), hh.hashValue(), d));
                res.add(hh);
                ++hm;
            }
            while (hm < hM) { // Still elements in the range
                // ++ med part
                rangeWidth = hM - hm;
if (debug) System.out.println("hm: " + hm + "; bits: " + Long.toBinaryString(hm));
                // get number of depth fille by the last 0 bits 
                int ddFromMin = hhCoder.nFilledDepthOnLastBits(d, (int) Long.numberOfTrailingZeros(hm));
if (debug) System.out.println("ddFromMin: " + ddFromMin);
                int ddFromRangeWidth = hhCoder.nFilledDepth(d, rangeWidth);
if (debug) System.out.println("ddFromRangeWidth: " + ddFromRangeWidth);
                int dd = ddFromRangeWidth < ddFromMin ? ddFromRangeWidth : ddFromMin;
                final HHash hhh = new HHImpl(d - dd, hhCoder.hash(d, hm, d - dd));
                res.add(hhh);
if (debug) System.out.println("Add2 " + hhh + " --> " + hhCoder.hashRange(hhh.depth(), hhh.hashValue(), d));
if (debug) System.out.println("nHash (from=" + ( d - dd + 1) + ", to= --> " + d + ") =  " + hhCoder.nHash(( d - dd + 1), d));
                hm += hhCoder.nHash(d - dd + 1, d); // 1 << (d << 1) = 2^(2*d) = 4^d
                //System.out.println("dd: " + dd + "; new hm = " + hm);
                //++ post part
                if (hM - hm < dNHash) {
                    while (hm < hM) {
                        final HHash hh = new HHImpl(d, hm);
if (debug) System.out.println("Add3 " + hh);
                        res.add(hh);
                        ++hm;
                    }
                }
            }
        }

  }*/
  
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
  
  /// Transform this (B)MOC as a simple (sorted) array of ranges.
  /// During the operation, we loose the `flag` information attached to each BMOC cell.  
  pub fn to_ranges(&self) -> Box<[std::ops::Range<u64>]> {
    let mut ranges: Vec<std::ops::Range<u64>> = Vec::with_capacity(self.entries.len());
    let mut prev_min = 0_u64;
    let mut prev_max = 0_u64;
    for cell in self.into_iter() {
      if cell.depth < self.depth_max {
        let range = to_range(cell.hash, self.depth_max - cell.depth);
        if range.start == prev_max {
          prev_max = range.end;
        } else {
          if prev_min != prev_max { // false only at first call, then alway true
            ranges.push(prev_min..prev_max);
          }
          prev_min = range.start;
          prev_max = range.end;
        }
      }  else {
        if cell.hash == prev_max {
          prev_max += 1;
        } else {
          if prev_min != prev_max { // false only at first call, then alway true
            ranges.push(prev_min..prev_max);
          }
          prev_min = cell.hash;
          prev_max = cell.hash + 1;
        }
      }
    }
    if prev_min != prev_max { // false only at first call, then alway true
      ranges.push(prev_min..prev_max);
    }
    ranges.into_boxed_slice()
  }
  
}

#[inline]
fn consume_while_overlapped(low_resolution: &Cell, iter: &mut BMOCIter) -> Option<Cell> {
  let mut cell = iter.next();
  while {
    match &cell {
      Some(c) => is_in(low_resolution,  c),
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
fn consume_while_overlapped_and_partial(low_resolution: &Cell, iter: &mut BMOCIter, res_is_overlapped: &mut bool) -> Option<Cell> {
  let mut cell = iter.next();
  while {
    match &cell {
      Some(c) => {
        if is_in(low_resolution,  c) {
          if c.is_full {
            *res_is_overlapped = true;
            false
          } else {
            true
          }
        } else {
          false
        }
      },
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
    && low_resolution.hash == (high_resolution.hash >> ( high_resolution.depth - low_resolution.depth) << 1)
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
