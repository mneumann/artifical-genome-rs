extern crate artificial_genome;

use artificial_genome::{Base, Genome};
// use artificial_genome::dna_base::DNABase;
use artificial_genome::base4::Base4;

fn main() {
    let genome = Genome {
        // genome: "A ATAT CG".chars().filter_map(|c| DNABase::from_char(c)).collect(),
        genome: "--gene A -- 0101 1111 --gene B-- 2222 0101 3101"
                    .chars()
                    .filter_map(|c| Base4::from_char(c))
                    .collect(),
    };

    // TATA box
    // let promoter = &[DNABase::A, DNABase::T, DNABase::A, DNABase::T];
    let promoter: Vec<_> = "0101".chars().filter_map(|c| Base4::from_char(c)).collect();

    let genes: Vec<_> = genome.iter_genes(&promoter, 4).collect();

    println!("{:?}", genes);
}
