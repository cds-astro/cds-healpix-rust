use crate::compass_point::{Cardinal, Ordinal};

pub struct ExternalEdge { // very simila to CardinalMap + OrdinalMap
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
  
  pub fn get_corner(&self, cardinal_point: &Cardinal) -> Option<u64> {
    match self.corners[cardinal_point.index() as usize] {
      Some(h) => Some(h),
      None => None,
    }
  }
  
  pub fn get_edge(&self, ordinal_point: &Ordinal) -> &[u64] {
    &self.edges[ordinal_point.index() as usize]
  }
  
}