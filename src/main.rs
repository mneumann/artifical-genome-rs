extern crate artificial_genome;

use artificial_genome::{BaseString, Genome};
use artificial_genome::base4::Base4;
use std::str::FromStr;

fn main() {
    let genome = Genome::<Base4>::from_str("--gene A -- 0101 1111 --gene B-- 2222 0101 3101")
                     .unwrap();
    let promoter = BaseString::<Base4>::from_str("0101").unwrap();

    let genes: Vec<_> = genome.iter_genes(&promoter, 4).collect();

    println!("{:?}", genes);
}
