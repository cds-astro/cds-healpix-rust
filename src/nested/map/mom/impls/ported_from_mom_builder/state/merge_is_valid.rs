/// Trait for checking if the merged state is valid.
pub trait MergeIsValid {
    /// Type of the leaf state.
    type State;

    /// Checks if merge from the given state to the merged state is valid.
    fn merge_is_valid(&self, original_states: &[Self::State], merged_state: &Self::State) -> bool;
}
