use super::Base;
use std::fmt;

#[derive(PartialEq, Eq, Copy, Clone)]
pub struct Base4(u8);

pub const B0: Base4 = Base4(0);
pub const B1: Base4 = Base4(1);
pub const B2: Base4 = Base4(2);
pub const B3: Base4 = Base4(3);

impl fmt::Debug for Base4 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Base for Base4 {
    fn succ(self) -> Self {
        Base4((self.0 + 1) & 3)
    }

    fn from_char(c: char) -> Option<Self> {
        match c {
            '0' => Some(Base4(0)),
            '1' => Some(Base4(1)),
            '2' => Some(Base4(2)),
            '3' => Some(Base4(3)),
            _ => None,
        }
    }
}
