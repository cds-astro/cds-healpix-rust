use crate::n_hash;
use crate::nested::map::skymap::implicit::ImplicitSkyMapArray;
use crate::nested::map::{
  skymap::{SkyMap, SkyMapValue},
  HHash,
};
use std::{
  collections::{
    btree_map::{IntoIter as BTreeMapIntoIter, Iter as BTreeMapIter, Values},
    BTreeMap,
  },
  iter::Map,
};

/// `BTreeMap` based structure: slow loading (tree creation), but fast queries.
/// To read from FITS file, an alternative is a sorted array with binary search.
#[derive(Debug)]
pub struct ExplicitSkyMapBTree<H: HHash, V: SkyMapValue> {
  depth: u8,
  entries: BTreeMap<H, V>,
  _zero: V,
}
impl<H: HHash, V: SkyMapValue> ExplicitSkyMapBTree<H, V> {
  pub fn new(depth: u8, entries: BTreeMap<H, V>) -> Self {
    Self {
      depth,
      entries,
      _zero: V::zero(),
    }
  }
  pub fn to_implicit_map(self) -> ImplicitSkyMapArray<H, V> {
    // We may implement an "Implicit iterator" to write a very large map on the disk
    // (with an implicit map possibly larger that the RAM).
    let depth = self.depth;
    let mut values = vec![V::zero(); n_hash(self.depth) as usize];
    for (k, v) in self.entries {
      values[k.to_usize().unwrap()] = v;
    }
    ImplicitSkyMapArray::new(depth, values.into_boxed_slice())
  }
}

impl<'a, H: HHash, V: SkyMapValue + 'a> SkyMap<'a> for ExplicitSkyMapBTree<H, V> {
  type HashType = H;
  type ValueType = V;
  type ValuesIt = Values<'a, H, V>;
  type EntriesIt = Map<BTreeMapIter<'a, H, V>, fn((&'a H, &'a V)) -> (H, &'a V)>;
  type OwnedEntriesIt = BTreeMapIntoIter<H, V>;

  fn depth(&self) -> u8 {
    self.depth
  }

  fn is_implicit(&self) -> bool {
    false
  }

  fn len(&self) -> usize {
    self.entries.len()
  }

  fn get(&self, hash: Self::HashType) -> &Self::ValueType {
    match self.entries.get(&hash) {
      Some(v) => v,
      None => &self._zero,
    }
  }

  fn values(&'a self) -> Self::ValuesIt {
    self.entries.values()
  }

  fn entries(&'a self) -> Self::EntriesIt {
    self.entries.iter().map(|(h, v)| (*h, v))
  }

  fn owned_entries(self) -> Self::OwnedEntriesIt {
    BTreeMap::into_iter(self.entries)
  }
}
