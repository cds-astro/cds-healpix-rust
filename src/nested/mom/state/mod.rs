//! Leaf state and merge rules.
//!
//! Merge logic is represented by [MergeStates] trait, which merges leaf states to their parent
//! nodes, or returns [None] if the states cannot be merged.
//!
//! Currently, the only type of the leaf state is implemented, [MinMaxMeanState],
//! with [MinMaxMeanStateMerger] merge rule,
//! which uses [MinMaxMeanStateValidator] to check if the result state is valid.
//! The state is represented by three values: minimum, maximum and mean.
//! [MinMaxMeanStateMerger] merges states by taking minimum and maximum of the states
//! and calculating the mean value.
//! [MinMaxMeanStateValidator] checks if the relative difference between minimum and maximum is less
//! than a given threshold.

pub mod merge_is_valid;
pub mod merge_states;
pub mod min_max_mean;
pub mod value;
