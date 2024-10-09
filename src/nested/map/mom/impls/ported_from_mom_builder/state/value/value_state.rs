use serde::{Deserialize, Serialize};

/// Simple leaf state with a single value.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValueState<T>(pub T);

macro_rules! impl_from_value_state {
    ($($t:ty),*) => {
        $(
            impl From<ValueState<$t>> for $t {
                fn from(state: ValueState<$t>) -> Self {
                    state.0
                }
            }
        )*
    };
    () => {};
}

impl_from_value_state!(bool, i8, i16, i32, i64, i128, u8, u16, u32, u64, u128, f32, f64, String);

impl<T> From<T> for ValueState<T> {
    fn from(val: T) -> Self {
        Self(val)
    }
}
