use std::convert::Infallible;
use std::error::Error;

use crate::nested::iter::checksorted::{
  CheckIsSorted, CheckIsSortedFromFallible, NotSortedItError,
};
use crate::{
  n_hash,
  nested::iter::sortedcoo2sortedhash::{SortedHashIt, SortedHashItFromFallibleIt},
};

/// Transforms a **sorted** iterator of HEALPix hash values into an iterator of tuple `(hash, count)`.
/// It is equivalent to:
/// * `dedup_with_count` in `itertools`;
/// * the `uniq -c` linux command.
///
/// # WARNING:
/// * This iterator does not check whether the input is sorted or not, but the only 2 `From` implemented
///   are supposed to ensure so.
/// * If you implemented a `new` method (or `From`) from generic Iterators, you must first decorate
///  them to ensure the input iterator is sorted.
pub struct SortedHash2HashCountIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  prev_hash: u64,
  prev_count: u32,
  not_depleted: bool,
  it: I,
}
impl<E, I> SortedHash2HashCountIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  fn new(it: I) -> Self {
    Self {
      prev_hash: 0,
      prev_count: 0,
      not_depleted: true,
      it: it.into(),
    }
  }
}

impl<J: Iterator<Item = u64>> From<J>
  for SortedHash2HashCountIt<NotSortedItError<u64, Infallible>, CheckIsSorted<u64, J>>
{
  fn from(it: J) -> Self {
    Self::new(it.into())
  }
}

impl<F: Error, J: Iterator<Item = Result<u64, F>>> From<J>
  for SortedHash2HashCountIt<NotSortedItError<u64, F>, CheckIsSortedFromFallible<u64, F, J>>
{
  fn from(it: J) -> Self {
    Self::new(it.into())
  }
}

impl<J> From<SortedHashIt<J>>
  for SortedHash2HashCountIt<NotSortedItError<u64, Infallible>, SortedHashIt<J>>
where
  J: Iterator<Item = (f64, f64)>,
{
  fn from(it: SortedHashIt<J>) -> Self {
    Self::new(it)
  }
}

impl<F, J> From<SortedHashItFromFallibleIt<F, J>>
  for SortedHash2HashCountIt<NotSortedItError<u64, F>, SortedHashItFromFallibleIt<F, J>>
where
  F: Error,
  J: Iterator<Item = Result<(f64, f64), F>>,
{
  fn from(it: SortedHashItFromFallibleIt<F, J>) -> Self {
    Self::new(it)
  }
}

impl<E, I> Iterator for SortedHash2HashCountIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  type Item = Result<(u64, u32), E>;

  fn next(&mut self) -> Option<Self::Item> {
    if self.not_depleted {
      while let Some(res_hash) = self.it.next() {
        match res_hash {
          Ok(hash) => {
            if hash == self.prev_hash {
              self.prev_count += 1;
            } else if self.prev_count != 0 {
              let res = Some(Ok((self.prev_hash, self.prev_count)));
              self.prev_hash = hash;
              self.prev_count = 1;
              return res;
            } else {
              // Called only at first iteration if first hash != 0
              self.prev_hash = hash;
              self.prev_count = 1;
            }
          }
          Err(err) => return Some(Err(err)),
        }
      }
      // We use 'not_depleted' to avoid having to call self.it once it already returned None.
      self.not_depleted = false;
      // No more elements in the iterator
      if self.prev_count != 0 {
        let res = Some(Ok((self.prev_hash, self.prev_count)));
        self.prev_count = 0; // useless now that we are using 'not_depleted'
        res
      } else {
        // Called only if the input iterator is empty
        None
      }
    } else {
      None
    }
  }
}

/// Transforms a **sorted** iterator of HEALPix hash values into an iterator of tuple `(hash, count)`.
/// It is equivalent to:
/// * `dedup_with_count` in `itertools`;
/// * the `uniq -c` linux command.
/// Except that all hash values from 0 to `n_cells(depth)` are returned, in order, with count = 0 for
/// hash value not present in the input iterator.
///
/// # WARNING:
/// * This iterator does not check whether the input is sorted or not, but the only 2 `From` implemented
///   are supposed to ensure so.
/// * If you implemented a `new` method (or `From`) from generic Iterators, you must first decorate
///  them to ensure the input iterator is sorted.
pub struct SortedHash2HashCountIncludingZeroIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  max_hash: u64,
  prev_hash: u64,
  curr_hash: u64,
  curr_count: u32,
  it: I,
}
impl<E, I> SortedHash2HashCountIncludingZeroIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  fn new(depth: u8, it: I) -> Self {
    Self {
      max_hash: n_hash(depth),
      prev_hash: 0,
      curr_hash: 0,
      curr_count: 0,
      it: it.into(),
    }
  }
}

impl<J>
  SortedHash2HashCountIncludingZeroIt<NotSortedItError<u64, Infallible>, CheckIsSorted<u64, J>>
where
  J: Iterator<Item = u64>,
{
  pub fn new_from_infallible(depth: u8, it: J) -> Self {
    Self::new(depth, it.into())
  }
}

impl<F, J>
  SortedHash2HashCountIncludingZeroIt<
    NotSortedItError<u64, F>,
    CheckIsSortedFromFallible<u64, F, J>,
  >
where
  F: Error,
  J: Iterator<Item = Result<u64, F>>,
{
  pub fn new_from_fallible(depth: u8, it: J) -> Self {
    Self::new(depth, it.into())
  }
}

impl<J> From<SortedHashIt<J>>
  for SortedHash2HashCountIncludingZeroIt<NotSortedItError<u64, Infallible>, SortedHashIt<J>>
where
  J: Iterator<Item = (f64, f64)>,
{
  fn from(it: SortedHashIt<J>) -> Self {
    Self::new(it.depth(), it.into())
  }
}

impl<F, J> From<SortedHashItFromFallibleIt<F, J>>
  for SortedHash2HashCountIncludingZeroIt<
    NotSortedItError<u64, F>,
    SortedHashItFromFallibleIt<F, J>,
  >
where
  F: Error,
  J: Iterator<Item = Result<(f64, f64), F>>,
{
  fn from(it: SortedHashItFromFallibleIt<F, J>) -> Self {
    Self::new(it.depth(), it.into())
  }
}

impl<E, I> Iterator for SortedHash2HashCountIncludingZeroIt<E, I>
where
  E: Error,
  I: Iterator<Item = Result<u64, E>>,
{
  type Item = Result<(u64, u32), E>;

  fn next(&mut self) -> Option<Self::Item> {
    if self.prev_hash < self.curr_hash {
      let res = Some(Ok((self.prev_hash, 0)));
      self.prev_hash += 1;
      res
    } else if self.curr_hash != self.max_hash {
      while let Some(res_hash) = self.it.next() {
        match res_hash {
          Ok(hash) => {
            if hash == self.curr_hash {
              self.curr_count += 1;
            } else {
              let res = Some(Ok((self.curr_hash, self.curr_count)));
              self.prev_hash = self.curr_hash + 1;
              self.curr_count = 1;
              self.curr_hash = hash;
              return res;
            }
          }
          Err(err) => return Some(Err(err)),
        }
      }
      let res = Some(Ok((self.curr_hash, self.curr_count)));
      self.prev_hash = self.curr_hash + 1;
      self.curr_count = 0; // useless, but does not cost much
      self.curr_hash = self.max_hash;
      res
    } else {
      None
    }
  }
}

#[cfg(test)]
mod tests {
  use std::iter;

  use super::*;

  #[test]
  fn test_count_iterators_ok() {
    // Test coutn from empty iterators
    for depth in 0..10 {
      println!("test depth: {}", depth);
      let count_it = SortedHash2HashCountIncludingZeroIt::new_from_infallible(depth, iter::empty());
      assert_eq!(
        count_it.filter(|it| it.is_ok()).count(),
        n_hash(depth) as usize
      )
    }
    // Test full iterators in with count = hash value
    for depth in 0..5 {
      for (h, c) in
        SortedHash2HashCountIt::from((0..n_hash(depth)).map(|i| vec![i; i as usize]).flatten())
          .map(|r| r.unwrap())
      {
        assert_eq!(h as u32, c);
      }
    }
    // Test full iterators in with count = hash value
    for depth in 0..5 {
      for (h, c) in SortedHash2HashCountIncludingZeroIt::new_from_infallible(
        depth,
        (0..n_hash(depth)).map(|i| vec![i; i as usize]).flatten(),
      )
      .map(|r| r.unwrap())
      {
        assert_eq!(h as u32, c);
      }
    }
    // Basic test
    let data = [2, 2, 2, 6, 6, 10, 99, 99];
    let hash_count_vec = SortedHash2HashCountIt::from(data.into_iter())
      .collect::<Result<Vec<(u64, u32)>, _>>()
      .unwrap();
    assert_eq!(hash_count_vec, [(2, 3), (6, 2), (10, 1), (99, 2)].to_vec());
    let hash_count_vec =
      SortedHash2HashCountIncludingZeroIt::new_from_infallible(4, data.into_iter())
        .collect::<Result<Vec<(u64, u32)>, _>>()
        .unwrap();
    assert_eq!(hash_count_vec.len(), n_hash(4) as usize);
    assert_eq!(
      hash_count_vec
        .iter()
        .filter(|(_, count)| *count == 0)
        .count(),
      n_hash(4) as usize - 4
    );
    assert_eq!(
      hash_count_vec
        .into_iter()
        .filter(|(_, count)| *count > 0)
        .collect::<Vec<(u64, u32)>>(),
      [(2, 3), (6, 2), (10, 1), (99, 2)].to_vec()
    );
  }

  #[test]
  fn test_count_iterators_error() {
    let data = [2, 2, 2, 6, 6, 99, 99, 10];
    let res_hash_count =
      SortedHash2HashCountIt::from(data.into_iter()).collect::<Result<Vec<(u64, u32)>, _>>();
    assert!(matches!(res_hash_count, Err(_)));
    let res_hash_count =
      SortedHash2HashCountIncludingZeroIt::new_from_infallible(4, data.into_iter())
        .collect::<Result<Vec<(u64, u32)>, _>>();
    assert!(matches!(res_hash_count, Err(_)));
  }
}
