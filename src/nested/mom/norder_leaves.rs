//! A data structure for storing the leaves of a tree of a given depth.

use crate::nested::mom::error::Error;
use std::ops::Range;

/// A data structure for storing the leaves of a tree of a given depth. The leaves are stored as a
/// pair of vectors: indexes and values. The indexes are stored in the ascending order.
#[derive(PartialEq, Debug)]
pub struct NorderLeaves<T> {
    indexes: Vec<usize>,
    values: Vec<T>,
}

impl<T> Default for NorderLeaves<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> NorderLeaves<T> {
    /// Create a new empty [NorderLeaves].
    pub fn new() -> Self {
        Self {
            indexes: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Inserts a value at the given index.
    ///
    /// The method panics if the index is not greater than the last inserted index.
    pub fn insert(&mut self, index: usize, value: T) -> Result<(), Error> {
        if !self.indexes.is_empty() && index <= self.indexes[0] {
            return Err(Error::IndexError);
        }
        self.indexes.push(index);
        self.values.push(value);
        Ok(())
    }

    /// Removes the last inserted value.
    ///
    /// The method panics if the [NorderLeaves] is empty.
    pub fn pop(&mut self) {
        self.indexes.pop();
        self.values.pop();
    }

    /// Returns the number of values stored in the [NorderLeaves].
    pub fn len(&self) -> usize {
        self.indexes.len()
    }

    /// Returns true if the [NorderLeaves] is empty.
    pub fn is_empty(&self) -> bool {
        self.indexes.is_empty()
    }

    /// Get the last inserted values according to the given indexes.
    ///
    /// The method returns [Err] if the given indexes are not equal to the last inserted indexes,
    /// and [Ok] with a slice of values otherwise.
    pub fn get_last_checked(&self, indexes: Range<usize>) -> Result<&[T], Error> {
        let last = self.len();
        if indexes.len() > last {
            return Err(Error::IndexError);
        }
        let first = last - indexes.len();
        if self.indexes[first..last]
            .iter()
            .zip(indexes.clone())
            .all(|(&a, b)| a == b)
        {
            let last = self.len();
            let first = last - indexes.len();
            Ok(&self.values[first..last])
        } else {
            Err(Error::IndexError)
        }
    }

    /// Converts the [NorderLeaves] into a pair of vectors: indexes and values.
    pub fn into_tuple(self) -> (Vec<usize>, Vec<T>) {
        (self.indexes, self.values)
    }
}
