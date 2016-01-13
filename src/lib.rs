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

fn locate_promoter<'a, T: Eq>(s: &'a [T], promoter: &[T]) -> Option<usize> {
    let plen = promoter.len();
    assert!(plen > 0);

    if s.len() < plen {
        return None;
    }

    for (i, window) in s.windows(plen).enumerate() {
        if window == promoter {
            return Some(i);
        }
    }

    return None;
}

#[derive(Debug)]
pub struct Gene<'a, B: Base + 'a> {
    pub regulatory_region: &'a [B],
    pub gene: &'a [B],
}

pub struct GeneIterator<'a, 'b, B: Base + 'a + 'b> {
    start_pos: usize,
    // Genes have fixed length
    length_of_gene: usize,
    genome: &'a [B],
    promoter: &'b [B],
}

impl<'a, 'b, B: Base + 'a + 'b> Iterator for GeneIterator<'a, 'b, B> {
    type Item = Gene<'a, B>;

    fn next(&mut self) -> Option<Self::Item> {
        match locate_promoter(&self.genome[self.start_pos..], self.promoter) {
            Some(pos) => {
                let reg_start = self.start_pos;
                let reg_end = reg_start + pos;
                let gene_start = reg_end + self.promoter.len();
                let gene_end = gene_start + self.length_of_gene;

                // gene is not complete
                if gene_end > self.genome.len() {
                    return None;
                }

                let gene = Gene {
                    regulatory_region: &self.genome[reg_start..reg_end],
                    gene: &self.genome[gene_start..gene_end],
                };
                self.start_pos = gene_end;
                return Some(gene);
            }
            None => {
                return None;
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
    pub fn iter_genes<'a, 'b>(&'a self,
                              promoter: &'b [B],
                              length_of_gene: usize)
                              -> GeneIterator<'a, 'b, B> {
        GeneIterator {
            start_pos: 0,
            length_of_gene: length_of_gene,
            genome: &self.genome,
            promoter: promoter,
        }
    }
}

#[test]
fn test_dnabase_next() {
    assert_eq!(DNABase::A, DNABase::C.next());
}
