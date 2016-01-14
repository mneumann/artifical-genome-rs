#![feature(slice_patterns)]

extern crate artificial_genome;

use artificial_genome::{Genome, ProteinRegulator};
use artificial_genome::base4::{Base4, B0, B1, B2, B3};
use std::str::FromStr;

fn main() {
    // let genome = Genome::<Base4>::from_str("<0101> 1111 | 2222 <0101> 3101")
    let genome = Genome::<Base4>::from_str("...11 _0320_23 <0101> T:0311 2...3 _1022_ 133 <0101> \
                                            W:3213 121...")
                     .unwrap();
    // let promoter = BaseString::<Base4>::from_str("0101").unwrap();
    let promoter = [B0, B1, B0, B1];

    let genes: Vec<_> = genome.iter_genes(&promoter, 4).collect();
    println!("{:?}", genes);

    let network = genome.construct_network(&promoter,
                                           4,
                                           &|product| {
                                               if product.last() == Some(&B0) {
                                                   // Inhibitor
                                                   ProteinRegulator::inhibit()
                                               } else {
                                                   // Inhibitor
                                                   ProteinRegulator::enhance()
                                               }
                                           },
                                           &|gene| {
                                               // XXX: Convert into u8
                                               match gene {
                                                   [B0, B3, B1, B1] => 'T',
                                                   [B3, B2, B1, B3] => 'W',
                                                   _ => 'X',
                                               }
                                           });

    println!("{:?}", network);

}
