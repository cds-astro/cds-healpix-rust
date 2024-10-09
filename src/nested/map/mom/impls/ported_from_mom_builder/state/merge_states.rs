/// Trait for merging leaf states to their parent nodes.
pub trait MergeStates {
    /// Type of the leaf state.
    type State: Sized;

    /// Merges the given states to a single state or returns [None] if the states cannot be merged.
    fn merge(&self, states: &[Self::State]) -> Option<Self::State>;
}
