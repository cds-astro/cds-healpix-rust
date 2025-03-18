//! Build tree from the leaf states of the same maximum depth.
//!
//! ### Index convention
//!
//! ```text
//!        __0__           __1__
//!       /     \         /     \
//!      0       1       2      3
//!     / \     / \     / \    / \
//!    0   1   2  3    4  5    6  7
//! ```  
//! (An example of a "tree" with `n_children=2`, `max_norder=2`, `n_root=2`)
//!
//! The index of the root nodes is arbitrary, but the index of the children is determined by the
//! index of the parent node, so the for the node of index `parent_index` the children have indexes
//! `[n_children * parent_index, n_children * (parent_index + 1))`.
//!
//! ### Data structure
//!
//! The tree is represented by a vector of [NorderLeaves], where each [NorderLeaves] represents
//! indexes and states of the leaves at a given depth (norder). The index of the vector is the
//! norder of the tree.
//!
//! ### Algorithm overview
//! 1. Tree builder consumes an iterator over the leaf states. Some of the leaves could be missing,
//!    but the iterator must be sorted by strictly increasing index. If the children group is
//!    complete, the builder tries to merge the states into the parent node. If the group is
//!    incomplete or the merge state is invalid the states are inserted into the tree as is.
//! 2. If the merge is successful, the builder adds the merged state to the parent node.
//! 3. Check if the parent node is a last child of the grandparent node. For example, for the sketch
//!    above, the middle level nodes 1 and 3 are last children of the root nodes 0 and 1. If we have
//!    just merged the last child of the grandparent node, we run attempt a recursive merge
//!    of the parent node with its siblings.
//! 4. The recursion interrupts when we cannot merge (children group is incomplete or merge is
//!    invalid) or we are at the root nodes.

use super::exclusive_option::ExclusiveOption;
use super::norder_leaves::NorderLeaves;
use super::state::merge_states::MergeStates;
use super::tree::{Tree, TreeMutRef};
use super::tree_config::TreeConfig;
use itertools::Itertools;

/// A structure to build a tree from a fallible iterator over the leaf states.
///
/// The tree is built bottom-up, i.e. the leaves are merged to their parent nodes recursively.
/// The merge rules are defined by the user, see [crate::state] module for the details. For now it
/// is exclusively used by [build_tree] function.
#[derive(Clone)]
pub struct TreeBuilder<Merger> {
    pub merger: Merger,
    pub config: TreeConfig,
}

impl<S, Merger> TreeBuilder<Merger>
where
    Merger: MergeStates<State = S>,
    S: Copy + std::fmt::Debug,
{
    /// Build a tree from a fallible iterator over the leaf states.
    /// See [build_tree] for the convention about `max_norder_states` argument.
    fn build<E>(
        &self,
        max_norder_states: impl IntoIterator<Item = Result<(usize, S), E>>,
    ) -> Result<Tree<S>, E> {
        let mut tree: Tree<S> = (0..=self.config.max_norder())
            .map(|_| NorderLeaves::new())
            .collect();
        let (max_norder_leaves, subtree) = tree.split_last_mut().expect("tree should not be empty");

        // Group states by Some(parent_index) if Ok, or a single-element group of None if Err
        // We use here custom ExclusiveOption instead of Option to make sure that errors are not
        // being grouped and the very first error is returned.
        let groups = max_norder_states
            .into_iter()
            .chunk_by(|result| match result {
                Ok((index, _state)) => ExclusiveOption::Some(*index / self.config.n_children()),
                Err(_) => ExclusiveOption::None,
            });

        for (parent_index, group) in groups.into_iter() {
            let (indexes, states) = {
                let it = group.into_iter();
                itertools::process_results(it, |iter| iter.unzip::<_, _, Vec<_>, Vec<_>>())?
            };
            let parent_index = parent_index.expect("Error has been processed before");

            if states.len() > self.config.n_children() {
                panic!("Too many states in one group, check if index is correct");
            }

            let merged = if states.len() == self.config.n_children() {
                // Run recursive merge
                self.merge_states(subtree, &states, parent_index)
            } else {
                false
            };

            // If not merged then insert the states into the current norder, max_norder
            if !merged {
                for (index, state) in indexes.into_iter().zip(states.into_iter()) {
                    max_norder_leaves
                        .insert(index, state)
                        .expect("Tree insertion failed");
                }
            }
        }

        Ok(tree)
    }

    /// Merges a group of states into a parent node with [MergeStates].
    /// If the states are mergeable [MergeStates.merge] returns [Some],
    /// the merged state is inserted into the tree
    /// and the recursive merge is run with [Self::process_group_by_top_title_index] for the
    /// parent node if it is the last child of the grandparent node.
    ///
    /// Arguments:
    /// - `subtree`: a mutable reference to the parent subtree, which is a vector of [NorderLeaves]
    ///   representing leaves at depth one less than the `states` one.
    /// - `states`: a slice of states to be merged. The length of the slice must be equal to the
    ///   number of children per node.
    /// - `parent_index`: the index of the parent node.
    ///
    /// Returns:
    /// - `true` if the states were merged, `false` otherwise.
    fn merge_states(&self, subtree: TreeMutRef<S>, states: &[S], parent_index: usize) -> bool {
        assert_eq!(states.len(), self.config.n_children());
        assert_ne!(subtree.len(), 0);

        if let Some(merged_state) = self.merger.merge(states) {
            subtree
                .last_mut()
                .expect("tree should not be empty")
                .insert(parent_index, merged_state)
                .expect("Tree insertion failed");

            if parent_index % self.config.n_children() == self.config.n_children() - 1 {
                self.process_group_by_top_title_index(subtree, parent_index);
            }
            true
        } else {
            false
        }
    }

    /// Process a group of child states determined by the index of the last child.
    ///
    /// It returns immediately if the tree has only one node (root node) and stops the recursion.
    /// Otherwise, it collects the states of the last child and its siblings and tries to merge them
    /// into the parent node. If the merge is successful, the pre-merged states are deleted.
    /// If the group is not complete, the function returns and no merge is performed.
    ///
    /// Arguments:
    /// - `tree`: a mutable reference to the subtree where the states are located.
    /// - `top_index`: the index of the last child.
    fn process_group_by_top_title_index(&self, tree: TreeMutRef<S>, top_index: usize) {
        // If we at norder=0 (root nodes), we are done
        if tree.len() == 1 {
            return;
        }

        let offset_from_bottom = self.config.n_children() - 1;
        assert_eq!(top_index % self.config.n_children(), offset_from_bottom);

        let (norder_leaves, subtree) = tree.split_last_mut().expect("tree should not be empty");

        // Collect pre-merged states
        // If any of the states is missing we cannot merge further, so return
        let states_result = norder_leaves
            .get_last_checked(top_index - offset_from_bottom..top_index + 1)
            .map(|sl| sl.to_vec());
        let states = match states_result {
            Ok(states) => states,
            Err(_) => return,
        };

        let parent_index = top_index / self.config.n_children();

        // Merge leaves if possible
        let merged = self.merge_states(subtree, &states, parent_index);
        // Delete pre-merged states if they were merged
        if merged {
            for _ in 0..self.config.n_children() {
                norder_leaves.pop();
            }
        }
    }
}

/// Build a tree from a fallible iterator over the leaf states.
///
/// Arguments:
/// - `state_builder`: a [StateBuilder] object, which specifies the merge rules.
/// - `tree_config`: a [TreeConfig] object, which specifies the tree structure: number of roots and
///   number of children per node
/// - `max_norder_states`: an object implementing [IntoIterator]. The iterator must yield
///   [Result]s, where [Ok] contains a tuple of `(index, state)`, and [Err] contains
///   an error. Leaves are assumed to be indexed over the range of
///   `[offset + 0, offset + max_norder_leaves)`, where `offset = const * max_norder_leaves`.
///   The iterator must be sorted by strictly increasing index, but it is allowed to skip some
///   indexes. If a leaf is missing, the leaf group will never be merged into a parent node and will
///   be kept as is.
///
/// Returns:
/// - [Ok] with a [Tree] object, if all the states are valid, and [Err] if one of the states has
///   failed. The function never returns [Err] for its own errors, but could panic.
pub fn build_tree<S, Merger, E>(
    merger: Merger,
    tree_config: TreeConfig,
    max_norder_states: impl IntoIterator<Item = Result<(usize, S), E>>,
) -> Result<Tree<S>, E>
where
    Merger: MergeStates<State = S>,
    S: Copy + std::fmt::Debug,
{
    TreeBuilder {
        merger,
        config: tree_config,
    }
    .build(max_norder_states)
}
