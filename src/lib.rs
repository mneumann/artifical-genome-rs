/// Represents the bases used in the genome string.
/// For example the bases of the DNA are adenine (A),
/// thymine (T), guanine (G) and cytosine (C).
trait Base: Sized {
    /// Returns the "next" base, wrapping around. Used
    /// to produce the gene product.
    fn next(self) -> Self;
}

#[derive(Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum DNABase {
    A,
    T,
    G,
    C
}

impl Base for DNABase {
    fn next(self) -> Self {
        match self {
            DNABase::A => DNABase::T,
            DNABase::T => DNABase::G,
            DNABase::G => DNABase::C,
            DNABase::C => DNABase::A
        }
    }
}

#[test]
fn test_dnabase_next() {
    assert_eq!(DNABase::A, DNABase::C.next());
}
