use crate::nested::mom::state::merge_states::MergeStates;
use crate::nested::mom::state::value::value_state::ValueState;
use serde::{Deserialize, Serialize};
use std::marker::PhantomData;

/// Merges leaf states if values are exactly equal.
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct ExactlyEqualMerger<T> {
    phantom: PhantomData<T>,
}

impl<T> ExactlyEqualMerger<T> {
    fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}

impl<T> Default for ExactlyEqualMerger<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> MergeStates for ExactlyEqualMerger<T>
where
    T: PartialEq + Copy,
{
    type State = ValueState<T>;

    fn merge(&self, states: &[Self::State]) -> Option<Self::State> {
        assert!(!states.is_empty());

        // Output any state if all states are equal, otherwise return None
        let first_state = states[0];
        if states.iter().skip(1).all(|&state| state == first_state) {
            Some(first_state)
        } else {
            None
        }
    }
}
