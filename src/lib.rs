pub mod dna_base;
pub mod base4;

/// Represents the bases used in the genome string.
/// For example the bases of the DNA are adenine (A),
/// thymine (T), guanine (G) and cytosine (C).
pub trait Base: Sized + PartialEq + Eq + Copy + Clone {
    /// Returns the "successor" base, wrapping around. Used
    /// to produce the gene product.
    fn succ(self) -> Self;

    fn from_char(c: char) -> Option<Self>;
}

// Locate ```substr``` in ```s```.
fn locate_substr<T: Eq>(s: &[T], substr: &[T]) -> Option<usize> {
    assert!(substr.len() > 0);

    if s.len() < substr.len() {
        return None;
    }

    for (i, window) in s.windows(substr.len()).enumerate() {
        if window == substr {
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

impl<'a, B: Base + 'a> Gene<'a, B> {
    pub fn gene_product(&self) -> Vec<B> {
        self.gene.iter().map(|b| b.succ()).collect()
    }

    pub fn find_product_in_regulatory_region(&self, product: &[B]) -> bool {
        locate_substr(self.regulatory_region, product).is_some()
    }
}

pub struct GeneIterator<'a, 'b, B: Base + 'a + 'b> {
    start_pos: usize,
    // Genes have fixed length
    length_of_gene: usize,
    sequence: &'a [B],
    promoter: &'b [B],
}

impl<'a, 'b, B: Base + 'a + 'b> Iterator for GeneIterator<'a, 'b, B> {
    type Item = Gene<'a, B>;

    fn next(&mut self) -> Option<Self::Item> {
        match locate_substr(&self.sequence[self.start_pos..], self.promoter) {
            Some(pos) => {
                let reg_start = self.start_pos;
                let reg_end = reg_start + pos;
                let gene_start = reg_end + self.promoter.len();
                let gene_end = gene_start + self.length_of_gene;

                // gene is not complete
                if gene_end > self.sequence.len() {
                    return None;
                }

                let gene = Gene {
                    regulatory_region: &self.sequence[reg_start..reg_end],
                    gene: &self.sequence[gene_start..gene_end],
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
            sequence: &self.genome,
            promoter: promoter,
        }
    }
}
