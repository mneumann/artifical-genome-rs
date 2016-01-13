extern crate artificial_genome;

use artificial_genome::{Genome, DNABase};

fn main() {
    let genome = Genome {
        genome: vec![
            DNABase::A, DNABase::T, DNABase::A, DNABase::T,
            DNABase::C, DNABase::G,
        ],
    };

    // TATA box
    let promoter = &[DNABase::A, DNABase::T, DNABase::A, DNABase::T];

    for genslice in genome.slices(promoter) {
        println!("{:?}", genslice);
    }
}
