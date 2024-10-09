//! Tree configuration.

use serde::{Deserialize, Serialize};

/// Tree configuration.
///
/// Actually, it is a configuration of the forest of trees, but for simplicity we call it
/// a multi-root tree.
#[derive(Clone, Serialize, Deserialize)]
pub struct TreeConfig {
    /// Number of roots.
    #[allow(dead_code)]
    n_root: usize,
    /// Number of children per node.
    n_children: usize,
    /// Maximum depth of the tree.
    max_norder: usize,
    /// Number of leaves on the maximum depth. This number represents the "full" tree, so it is the
    /// maximum number of leaves. The actual number of leaves could be less if some leaves are
    /// been merged or missed.
    max_norder_nleaves: usize,
}

impl TreeConfig {
    /// Creates a new [TreeConfig] with  the given number of roots, number of children per node and
    /// maximum depth.
    pub fn new(
        n_root: impl Into<usize>,
        n_children: impl Into<usize>,
        max_norder: impl Into<usize>,
    ) -> Self {
        let n_root = n_root.into();
        let n_children = n_children.into();
        let max_norder = max_norder.into();
        let max_norder_nleaves = n_root * n_children.pow(max_norder as u32);

        Self {
            n_root,
            n_children,
            max_norder,
            max_norder_nleaves,
        }
    }

    /// Returns the number of roots.
    pub fn n_children(&self) -> usize {
        self.n_children
    }

    /// Returns the number of children per node.
    pub fn max_norder(&self) -> usize {
        self.max_norder
    }

    /// Returns the maximum number of leaves the tree can have.
    pub fn max_norder_nleaves(&self) -> usize {
        self.max_norder_nleaves
    }
}
