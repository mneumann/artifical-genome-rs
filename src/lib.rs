extern crate fixedbitset;
extern crate rand;

pub mod dna_base;
pub mod base4;
pub mod graph;

use std::str::FromStr;
use std::ops::Deref;
use std::fmt::{self, Debug};
use fixedbitset::FixedBitSet;
use rand::{Rng, Rand};

/// Represents the bases used in the genome string.
/// For example the bases of the DNA are adenine (A),
/// thymine (T), guanine (G) and cytosine (C).
pub trait Base: Sized + PartialEq + Eq + Copy + Clone + Debug + Rand {
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

// Count occurrence of ```substr``` in ```s```.
fn count_substr<T: Eq>(s: &[T], substr: &[T]) -> usize {
    assert!(substr.len() > 0);

    let mut cnt = 0;
    for window in s.windows(substr.len()) {
        if window == substr {
            cnt += 1;
        }
    }

    return cnt;
}


#[derive(Debug)]
pub struct Gene<'a, B: Base + 'a> {
    pub regulatory_region: &'a [B],
    pub gene: &'a [B],
}

impl<'a, B: Base + 'a> Gene<'a, B> {
    /// The gene product
    pub fn product(&self) -> BaseString<B> {
        BaseString { v: self.gene.iter().map(|b| b.succ()).collect() }
    }

    pub fn find_product_in_regulatory_region(&self, product: &[B]) -> bool {
        locate_substr(self.regulatory_region, product).is_some()
    }
    pub fn count_product_in_regulatory_region(&self, product: &[B]) -> usize {
        count_substr(self.regulatory_region, product)
    }
}

pub struct GeneIterator<'a, 'b, B: Base + 'a + 'b> {
    // Genes have fixed length
    length_of_gene: usize,
    sequence: &'a [B],
    promoter: &'b [B],
}

impl<'a, 'b, B: Base + 'a + 'b> Iterator for GeneIterator<'a, 'b, B> {
    type Item = Gene<'a, B>;

    fn next(&mut self) -> Option<Self::Item> {
        match locate_substr(self.sequence, self.promoter) {
            Some(pos) => {
                let gene_start = pos + self.promoter.len();
                let gene_end = gene_start + self.length_of_gene;

                // gene is not complete
                if gene_end > self.sequence.len() {
                    return None;
                }

                let gene = Gene {
                    regulatory_region: &self.sequence[..pos],
                    gene: &self.sequence[gene_start..gene_end],
                };
                self.sequence = &self.sequence[gene_end..];
                return Some(gene);
            }
            None => {
                return None;
            }
        }
    }
}

// A positive value enhances, a negative inhibits the expression of a gene.
#[derive(Debug)]
pub struct ProteinRegulator(i32);

impl ProteinRegulator {
    pub fn enhance() -> ProteinRegulator {
        ProteinRegulator(1)
    }

    pub fn inhibit() -> ProteinRegulator {
        ProteinRegulator(-1)
    }
}

#[derive(Clone)]
pub struct BaseString<B: Base> {
    v: Vec<B>,
}

impl<B: Base> FromStr for BaseString<B> {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BaseString { v: s.chars().filter_map(|c| B::from_char(c)).collect() })
    }
}

impl<B: Base> Deref for BaseString<B> {
    type Target = [B];

    fn deref(&self) -> &Self::Target {
        &self.v
    }
}

impl<B: Base> Debug for BaseString<B> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "["));
        for base in &self.v {
            try!(write!(f, "{:?}", base));
        }
        write!(f, "]")
    }
}

impl<B: Base> BaseString<B> {
    pub fn random<R: Rng>(rng: &mut R, n: usize) -> BaseString<B> {
        assert!(n > 0);
        BaseString { v: (0..n).map(|_| rng.gen()).collect() }
    }
}


/// A Genome is a string of Base
#[derive(Debug, Clone)]
pub struct Genome<B: Base> {
    genome: BaseString<B>,
}

impl<B: Base> Genome<B> {
    pub fn from_vec(v: Vec<B>) -> Genome<B> {
        Genome { genome: BaseString { v: v } }
    }

    pub fn random<R: Rng>(rng: &mut R, n: usize) -> Genome<B> {
        Genome { genome: BaseString::random(rng, n) }
    }
}

impl<B: Base> Deref for Genome<B> {
    type Target = [B];

    fn deref(&self) -> &Self::Target {
        &self.genome
    }
}

#[derive(Debug)]
struct Edge {
    src: usize,
    weight: ProteinRegulator,
}

#[derive(Debug)]
pub struct Node {
    incoming_edges: Vec<Edge>,
}

impl Node {
    fn new() -> Node {
        Node { incoming_edges: Vec::new() }
    }

    pub fn sum_edges(&self, network_state: &GeneNetworkState) -> i32 {
        let mut sum = 0;
        for edge in self.incoming_edges.iter() {
            let factor = if network_state.state.contains(edge.src) {
                1
            } else {
                0
            };
            sum += factor * edge.weight.0;
        }

        return sum;
    }
}

#[derive(Debug)]
pub struct GeneNetwork {
    nodes: Vec<Node>,
}

#[derive(Debug, Clone)]
pub struct GeneNetworkState {
    pub state: FixedBitSet,
}

impl GeneNetwork {
    fn new(num_nodes: usize) -> GeneNetwork {
        assert!(num_nodes > 0);
        GeneNetwork { nodes: (0..num_nodes).map(|_| Node::new()).collect() }
    }

    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    fn add_edge(&mut self, src: usize, dst: usize, weight: ProteinRegulator) {
        assert!(src < self.nodes.len());
        assert!(dst < self.nodes.len());
        self.nodes[dst].incoming_edges.push(Edge {
            src: src,
            weight: weight,
        });
    }

    pub fn new_state(&self) -> GeneNetworkState {
        GeneNetworkState { state: FixedBitSet::with_capacity(self.nodes.len()) }
    }
}

// Convert genome into sections, i.e. Split at the promoter.
impl<B: Base> Genome<B> {
    pub fn iter_genes<'a, 'b>(&'a self,
                              promoter: &'b [B],
                              length_of_gene: usize)
                              -> GeneIterator<'a, 'b, B> {
        GeneIterator {
            length_of_gene: length_of_gene,
            sequence: &self.genome,
            promoter: promoter,
        }
    }

    // Construct a dependency network between the genes
    pub fn construct_network<F>(&self,
                                promoter: &[B],
                                length_of_gene: usize,
                                protein_regulation: &F)
                                -> Option<GeneNetwork>
        where F: Fn(&[B]) -> ProteinRegulator
    {
        let genes: Vec<_> = self.iter_genes(promoter, length_of_gene).collect();
        let num_genes = genes.len();

        if num_genes == 0 {
            return None;
        }

        // each gene is a node in the boolean network
        let mut network = GeneNetwork::new(num_genes);

        for (src, gene) in genes.iter().enumerate() {
            let product = gene.product();
            // A gene product either enhances (> 0) or inyhibits (< 0) the expression of
            // another gene.
            let regulator = protein_regulation(&product);

            // determine which other genes ```gene``` regulates
            for (dst, gene2) in genes.iter().enumerate() {
                // XXX: Can a gene regulate itself?
                // if src != dst {
                let factor = gene2.count_product_in_regulatory_region(&product);
                if factor > 0 {
                    network.add_edge(src, dst, ProteinRegulator(regulator.0 * factor as i32));
                }
                // }
            }
        }

        return Some(network);
    }
}

impl<B: Base> FromStr for Genome<B> {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        FromStr::from_str(s).map(|bs| Genome { genome: bs })
    }
}
