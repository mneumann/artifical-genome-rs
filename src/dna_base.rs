use super::Base;
use rand::{Rand, Rng};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
#[repr(u8)]
pub enum DNABase {
    A,
    T,
    G,
    C,
}

impl Rand for DNABase {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        *rng.choose(&[DNABase::A, DNABase::T, DNABase::G, DNABase::C]).unwrap()
    }
}

impl Base for DNABase {
    fn succ(self) -> Self {
        match self {
            DNABase::A => DNABase::T,
            DNABase::T => DNABase::G,
            DNABase::G => DNABase::C,
            DNABase::C => DNABase::A,
        }
    }

    fn from_char(c: char) -> Option<Self> {
        match c {
            'A' => Some(DNABase::A),
            'T' => Some(DNABase::T),
            'G' => Some(DNABase::G),
            'C' => Some(DNABase::C),
            _ => None,
        }
    }
}

#[test]
fn test_dnabase_succ() {
    assert_eq!(DNABase::A, DNABase::C.succ());
}
