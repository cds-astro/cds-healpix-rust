use std::{convert::Infallible, error::Error};

use crate::nested::{get, iter::checksorted::NotSortedItError, Layer};

/// The internal Iterator must be an iterator on `(lon, lat)` tuples in radians.
/// This Iterator fails (i.e. returns an Error) if the computed value are not ordered
/// by increasing `hash` values.
pub struct SortedHashIt<I>
where
  I: Iterator<Item = (f64, f64)>,
{
  layer: &'static Layer,
  prev: u64,
  coo_it: I,
}

impl<I> SortedHashIt<I>
where
  I: Iterator<Item = (f64, f64)>,
{
  pub fn new(depth: u8, coo_it: I) -> Self {
    Self {
      layer: get(depth),
      prev: 0,
      coo_it,
    }
  }
  pub fn depth(&self) -> u8 {
    self.layer.depth
  }
}

impl<I> Iterator for SortedHashIt<I>
where
  I: Iterator<Item = (f64, f64)>,
{
  type Item = Result<u64, NotSortedItError<u64, Infallible>>;

  fn next(&mut self) -> Option<Self::Item> {
    match self.coo_it.next() {
      Some((lon, lat)) => {
        let curr = self.layer.hash(lon, lat);
        if self.prev <= curr {
          self.prev = curr;
          Some(Ok(curr))
        } else {
          Some(Err(NotSortedItError::NotSorted {
            prev: self.prev,
            curr,
          }))
        }
      }
      None => None,
    }
  }
}

/// The internal Iterator must be an iterator on `(lon, lat)` tuples in radians.
pub struct SortedHashItFromFallibleIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<(f64, f64), E>>,
{
  layer: &'static Layer,
  prev: u64,
  res_coo_it: I,
}

impl<E, I> SortedHashItFromFallibleIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<(f64, f64), E>>,
{
  pub fn new(depth: u8, res_coo_it: I) -> Self {
    Self {
      layer: get(depth),
      prev: 0,
      res_coo_it,
    }
  }
  pub fn depth(&self) -> u8 {
    self.layer.depth
  }
}

impl<E, I> Iterator for SortedHashItFromFallibleIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<(f64, f64), E>>,
{
  type Item = Result<u64, NotSortedItError<u64, E>>;

  fn next(&mut self) -> Option<Self::Item> {
    match self.res_coo_it.next() {
      Some(Ok((lon, lat))) => {
        let curr = self.layer.hash(lon, lat);
        if self.prev <= curr {
          self.prev = curr;
          Some(Ok(curr))
        } else {
          Some(Err(NotSortedItError::NotSorted {
            prev: self.prev,
            curr,
          }))
        }
      }
      Some(Err(e)) => Some(Err(NotSortedItError::Other(e))),
      None => None,
    }
  }
}
