//! Tree and related types.

use crate::nested::mom::norder_leaves::NorderLeaves;

/// Tree is a vector of [NorderLeaves] representing leaves at a given depth. Index denotes the norder
/// (depth) of the tree.
pub type Tree<S> = Vec<NorderLeaves<S>>;

/// Mutable reference to a tree.
pub type TreeMutRef<'a, S> = &'a mut [NorderLeaves<S>];
