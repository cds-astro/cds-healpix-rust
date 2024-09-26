use std::{iter::Map, slice::Iter};

use super::super::{super::skymap::SkyMapValue, Mom, ZUniqHashT};

/// Implementation of a MOM in an ordered vector of `(zuniq, values)` tuples.
pub struct MomVecImpl<Z, V>
where
  Z: ZUniqHashT,
  V: SkyMapValue,
{
  depth: u8,
  entries: Vec<(Z, V)>,
}
impl<'a, Z, V> Mom<'a> for MomVecImpl<Z, V>
where
  Z: ZUniqHashT,
  V: 'a + SkyMapValue,
{
  type ZUniqHType = Z;
  type ValueType = V;
  type ZuniqIt = Map<Iter<'a, (Z, V)>, fn(&'a (Z, V)) -> Z>;
  type EntriesIt = Map<Iter<'a, (Z, V)>, fn(&'a (Z, V)) -> (Z, &'a V)>;

  fn depth_max(&self) -> u8 {
    self.depth
  }

  fn get_cell_containing_unsafe(
    &'a self,
    hash_at_depth_max: Self::ZUniqHType,
  ) -> Option<(Self::ZUniqHType, &'a Self::ValueType)> {
    match self
      .entries
      .binary_search_by(|&(z, _)| z.cmp(&hash_at_depth_max))
    {
      Ok(i) => {
        let e = &self.entries[i];
        Some((e.0, &e.1))
      }
      Err(i) => {
        if i > 0 {
          // if array len is 0, i will be 0 so we do not enter here.
          let e = &self.entries[i - 1];
          if Z::are_overlapping(hash_at_depth_max, e.0) {
            return Some((e.0, &e.1));
          }
        }
        if i < self.entries.len() {
          let e = &self.entries[i];
          if Z::are_overlapping(hash_at_depth_max, e.0) {
            return Some((e.0, &e.1));
          }
        }
        None
      }
    }
  }

  fn get_overlapped_cells(
    &'a self,
    zuniq: Self::ZUniqHType,
  ) -> Vec<(Self::ZUniqHType, &'a Self::ValueType)> {
    let mut range = match self.entries.binary_search_by(|&(z, _)| z.cmp(&zuniq)) {
      Ok(i) => i..i + 1,
      Err(i) => i..i,
    };
    while range.start - 1 > 0 && Z::are_overlapping(zuniq, self.entries[range.start - 1].0) {
      range.start -= 1;
    }
    while range.end < self.entries.len() && Z::are_overlapping(zuniq, self.entries[range.end].0) {
      range.end += 1;
    }
    range
      .into_iter()
      .map(|i| {
        let (z, v) = &self.entries[i];
        (*z, v)
      })
      .collect()
  }

  fn zuniqs(&'a self) -> Self::ZuniqIt {
    self.entries.iter().map(|&(zuniq, _)| zuniq)
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self.entries.iter().map(|(z, v)| (*z, v))
  }
}
