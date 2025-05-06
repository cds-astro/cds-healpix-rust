use crate::nested::map::mom::impls::ported_from_mom_builder::state::merge_is_valid::MergeIsValid;
use crate::nested::map::mom::impls::ported_from_mom_builder::state::merge_states::MergeStates;
use super::ValueState;
use num_traits::Zero;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Sum values until merge state is valid.
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SumUntilMerger<V> {
    pub validator: V,
}

impl<V> SumUntilMerger<V> {
    pub fn new(validator: V) -> Self {
        Self { validator }
    }
}

impl<V> Default for SumUntilMerger<V>
where
    V: Default,
{
    fn default() -> Self {
        Self::new(V::default())
    }
}

impl<T, V> MergeStates for SumUntilMerger<V>
where
    T: Copy + std::iter::Sum,
    V: MergeIsValid<State = ValueState<T>>,
{
    type State = ValueState<T>;

    fn merge(&self, states: &[Self::State]) -> Option<Self::State> {
        assert!(!states.is_empty());

        let sum_value = states.iter().map(|state| state.0).sum();
        let merged_state = ValueState(sum_value);

        if self.validator.merge_is_valid(states, &merged_state) {
            Some(merged_state)
        } else {
            None
        }
    }
}

/// Validator for maximum value, which checks if the merge state is less than or equal to the maximum value.
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MaximumValueValidator<T>(pub T);

impl<T> MergeIsValid for MaximumValueValidator<T>
where
    T: PartialOrd,
{
    type State = ValueState<T>;

    /// Checks if the merged state is valid, i.e. less than or equal to the maximum value.
    fn merge_is_valid(&self, _original_states: &[Self::State], merged_state: &Self::State) -> bool {
        merged_state.0 <= self.0
    }
}

/// Validator for maximum value, which also checks if there are zero-value children states.
///
/// This validator allows merge if 1) the merged state is less than or equal to the maximum value
/// and 2) number of empty siblings is less than or equal to the maximum number of empty siblings.
/// Also, you can optionally allow merge if all children states are zero. We refer to the children
/// states with zero value as "empty siblings".
#[derive(Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MaximumValueEmptySiblingValidator<T> {
    /// Maximum value of the merged state.
    pub threshold: T,
    /// Maximum number of empty siblings to allow merge.
    /// Setting it larger than the number of children would lead to the same behavior as
    /// [MaximumValueValidator].
    pub max_empty_siblings: usize,
    /// Allow merge if all children states are zero, even if max_empty_siblings is exceeded.
    pub allow_empty_merge: bool,
}

impl<T> MergeIsValid for MaximumValueEmptySiblingValidator<T>
where
    T: PartialOrd + Zero,
{
    type State = ValueState<T>;

    /// Checks if the merged state is valid, i.e. less than or equal to the maximum value and number
    /// of empty siblings is less than or equal to the maximum number of empty siblings. If
    /// `allow_empty_merge` is true, it also allows merge if all children states are zero.
    fn merge_is_valid(&self, original_states: &[Self::State], merged_state: &Self::State) -> bool {
        let empty_siblings = original_states
            .iter()
            .filter(|state| state.0 == T::zero())
            .count();

        let bellow_threshold = merged_state.0 <= self.threshold;

        if self.allow_empty_merge && empty_siblings == original_states.len() {
            return bellow_threshold;
        } else {
            // This else here is to show that this if should always go after the previous if
            if empty_siblings > self.max_empty_siblings {
                return false;
            }
        }

        bellow_threshold
    }
}
