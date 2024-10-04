//! [Option] reimplementation with [PartialEq] implementation to make [None] variants distinct

/// Naive [Option] reimplementation with [PartialEq] implementation to make [None] variants
/// to be non-equal (exclusive).
///
/// We use it to make sure that errors are not being grouped and the very first error is returned.
pub enum ExclusiveOption<T> {
    Some(T),
    None,
}

impl<T> ExclusiveOption<T> {
    pub fn expect(self, msg: &str) -> T {
        match self {
            Self::Some(val) => val,
            Self::None => panic!("{}", msg),
        }
    }
}

impl<T> PartialEq for ExclusiveOption<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Some(a), Self::Some(b)) => a == b,
            _ => false,
        }
    }
}
