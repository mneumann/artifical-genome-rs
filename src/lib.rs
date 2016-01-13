/// Represents the bases used in the genome string.
/// For example the bases of the DNA are adenine (A),
/// thymine (T), guanine (G) and cytosine (C).
pub trait Base: Sized + PartialEq + Eq {
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
    C,
}

impl Base for DNABase {
    fn next(self) -> Self {
        match self {
            DNABase::A => DNABase::T,
            DNABase::T => DNABase::G,
            DNABase::G => DNABase::C,
            DNABase::C => DNABase::A,
        }
    }
}

pub struct GenomeSlicer<'a, 'b, B: Base + 'a + 'b> {
    start_pos: usize,
    current: &'a [B],
    promoter: &'b [B],
}

fn split_slice_at_promoter<'a, T: Eq>(s: &'a [T],
                                      start_pos: usize,
                                      promoter: &[T])
                                      -> (&'a [T], Option<usize>) {
    let plen = promoter.len();
    assert!(plen > 0);

    if start_pos >= s.len() {
        return (&s[0..0], None);
    }

    let s = &s[start_pos..];
    for (i, window) in s.windows(plen).enumerate() {
        if window == promoter {
            return (&s[..i], Some(start_pos + i + plen));
        }
    }

    return (s, None);
}


// cuts the genome into slices
impl<'a, 'b, B: Base + 'a + 'b> Iterator for GenomeSlicer<'a, 'b, B> {
    type Item = &'a [B];

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match split_slice_at_promoter(self.current, self.start_pos, self.promoter) {
                (s, None) => {
                    if s.is_empty() {
                        return None;
                    }
                    self.start_pos = self.current.len();
                    return Some(s);
                }
                (s, Some(next_pos)) => {
                    self.start_pos = next_pos;
                    if s.len() > 0 {
                        return Some(s);
                    }
                }
            }
        }
    }
}

/// A Genome is a string of Base
pub struct Genome<B: Base> {
    pub genome: Vec<B>,
}

// Convert genome into sections, i.e. Split at the promoter.
impl<B: Base> Genome<B> {
    pub fn slices<'a, 'b>(&'a self, promoter: &'b [B]) -> GenomeSlicer<'a, 'b, B> {
        GenomeSlicer {
            start_pos: 0,
            current: &self.genome,
            promoter: promoter,
        }
    }
}

#[test]
fn test_dnabase_next() {
    assert_eq!(DNABase::A, DNABase::C.next());
}
