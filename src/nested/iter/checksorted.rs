use std::{
  convert::Infallible,
  error::Error,
  fmt::{Debug, Display, Formatter},
};

#[derive(Clone, Copy)]
pub enum NotSortedItError<V, E>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  E: Error,
{
  NotSorted { prev: V, curr: V },
  Other(E),
}

impl<V, E> Debug for NotSortedItError<V, E>
where
  E: Error,
  V: Clone + Copy + Display + PartialEq + PartialOrd,
{
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match self {
      Self::NotSorted { prev, curr } => write!(
        f,
        "Elements in the iterator are not sorted: {} > {}.",
        &prev, &curr
      ),

      Self::Other(e) => Debug::fmt(e, f),
    }
  }
}

impl<V, E> Display for NotSortedItError<V, E>
where
  E: Error,
  V: Clone + Copy + Display + PartialEq + PartialOrd,
{
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match self {
      Self::NotSorted { prev, curr } => {
        write!(
          f,
          "Elements in the iterator are not sorted: {} > {}.",
          &prev, &curr
        )
      }
      Self::Other(e) => Display::fmt(e, f),
    }
  }
}

impl<V, E> Error for NotSortedItError<V, E>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  E: Error,
{
}

/// Decorates an iterator, checking whether it is sorted or not while iterating.
pub struct CheckIsSorted<V, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  I: Iterator<Item = V>,
{
  prev: Option<V>,
  it: I,
}
impl<V, I> From<I> for CheckIsSorted<V, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  I: Iterator<Item = V>,
{
  fn from(it: I) -> Self {
    Self { prev: None, it }
  }
}

impl<V, I> Iterator for CheckIsSorted<V, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  I: Iterator<Item = V>,
{
  type Item = Result<V, NotSortedItError<V, Infallible>>;

  fn next(&mut self) -> Option<Self::Item> {
    match self.it.next() {
      Some(curr) => match self.prev.replace(curr) {
        Some(prev) => {
          if prev.le(&curr) {
            Some(Ok(curr))
          } else {
            Some(Err(NotSortedItError::NotSorted { prev, curr }))
          }
        }
        None => Some(Ok(curr)),
      },
      None => None,
    }
  }
}

/// Decorates an iterator, checking whether it is sorted or not while iterating.
pub struct CheckIsSortedFromFallible<V, E, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  E: Error,
  I: Iterator<Item = Result<V, E>>,
{
  prev: Option<V>,
  it: I,
}

impl<V, E, I> From<I> for CheckIsSortedFromFallible<V, E, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  E: Error,
  I: Iterator<Item = Result<V, E>>,
{
  fn from(it: I) -> Self {
    Self { prev: None, it }
  }
}

impl<V, E, I> Iterator for CheckIsSortedFromFallible<V, E, I>
where
  V: Display + Clone + Copy + PartialEq + PartialOrd,
  E: Error,
  I: Iterator<Item = Result<V, E>>,
{
  type Item = Result<V, NotSortedItError<V, E>>;

  fn next(&mut self) -> Option<Self::Item> {
    match self.it.next() {
      Some(Ok(curr)) => match self.prev.replace(curr) {
        Some(prev) => {
          if prev.le(&curr) {
            Some(Ok(curr))
          } else {
            Some(Err(NotSortedItError::NotSorted { prev, curr }))
          }
        }
        None => Some(Ok(curr)),
      },
      Some(Err(e)) => Some(Err(NotSortedItError::Other(e))),
      None => None,
    }
  }
}
