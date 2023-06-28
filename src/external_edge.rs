//! Module simply defining a the structure storing the external edges of an HEALPix cell and
//! provided accesses according to the main wind directions.
//! This strcut is common to nested and ring scheme.

use crate::compass_point::{Cardinal, Ordinal};

/// Stores the external edges of an HEALPix cell, providing accessor using the main wind directions.
#[derive(Debug)]
pub struct ExternalEdge {
  // very simila to CardinalMap + OrdinalMap
  corners: [Option<u64>; 4],
  edges: [Box<[u64]>; 4],
}

impl ExternalEdge {
  pub(crate) fn new_empty() -> Self {
    ExternalEdge {
      corners: Default::default(),
      edges: Default::default(),
    }
  }

  pub(crate) fn set_corner(&mut self, cardinal_direction: &Cardinal, value: u64) {
    self.corners[cardinal_direction.index() as usize] = Some(value);
  }

  pub(crate) fn set_edge(&mut self, ordinal_direction: &Ordinal, values: Box<[u64]>) {
    self.edges[ordinal_direction.index() as usize] = values;
  }

  /// Returns the neighbour cell located at the given cardinal point (if it exists).
  pub fn get_corner(&self, cardinal_point: &Cardinal) -> Option<u64> {
    self.corners[cardinal_point.index() as usize]
  }

  /// Returns the neighbour cells along the given ordinal direction.
  pub fn get_edge(&self, ordinal_point: &Ordinal) -> &[u64] {
    &self.edges[ordinal_point.index() as usize]
  }
}
