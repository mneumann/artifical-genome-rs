extern crate artificial_genome;
extern crate rand;

use artificial_genome::Genome;
use artificial_genome::base4::Base4;
use artificial_genome::graph::graph_from_base4_genome;

fn main() {
    let mut rng = rand::thread_rng();

    let genome = Genome::<Base4>::random(&mut rng, 8 * 256);
    let graph = graph_from_base4_genome(&genome, 5);
    println!("{:#?}", graph);
}
