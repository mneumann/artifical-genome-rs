#![feature(slice_patterns)]
#![allow(dead_code)]

extern crate artificial_genome;
extern crate rand;

use artificial_genome::{Genome, ProteinRegulator, GeneNetwork, GeneNetworkState};
use artificial_genome::base4::{Base4, B0, B1};
// use std::str::FromStr;
use std::mem;
use std::borrow::Cow;
use std::io::{self, Write};

#[derive(Debug, Clone)]
struct Edge {
    src_node: usize,
    dst_node: usize,
    length: f32,

    // there is a structural node, and each time it is activated,
    // this count is increased.
    type_count: usize,

    // Every edge has it's own regulatory network state embedded.
    network_state: GeneNetworkState,
}

const EDGE_DIFFERENTIATION: usize = 0;
const EDGE_SPLIT: usize = 1;
const EDGE_DUPLICATE: usize = 2;
const EDGE_SWAP: usize = 3;
const EDGE_GROW: usize = 4;
const EDGE_SHRINK: usize = 5;
const EDGE_TYPE: usize = 6;

impl Edge {
    fn transition_state(&mut self, network: &GeneNetwork) {
        // create an empty new state.
        let mut new_state = network.new_state();

        for (i, node) in network.nodes().iter().enumerate() {
            // determine the new state of ```node``` (position i in new_state)
            if node.sum_edges(&self.network_state) > 0 {
                // node is enabled
                new_state.state.insert(i);
            }
        }
        self.network_state = new_state;
    }

    fn develop(&mut self,
               network: &GeneNetwork,
               next_node_id: &mut usize,
               new_edges: &mut Vec<Edge>) {

        // first perform the state transition.
        self.transition_state(network);

        // perform the actions of all active nodes in the gene network.

        // first node is used for differentiation.
        // in a split/duplicate node, this will be set 0 or 1 in the children.
        // starting from second node, the nodes are mapped to graph grammar rules.

        if self.network_state.state.contains(EDGE_SPLIT) {
            let new_node = *next_node_id;
            *next_node_id += 1;

            // differentiate
            let mut child_state = self.network_state.clone();
            child_state.state.set(EDGE_DIFFERENTIATION,
                                  !self.network_state.state.contains(EDGE_DIFFERENTIATION));

            let new_edge = Edge {
                src_node: new_node,
                dst_node: self.dst_node,
                length: 0.5 * self.length,
                network_state: child_state,
                type_count: self.type_count, /* XXX: start with 0 or with the count of the current edge? */
            };
            self.dst_node = new_node;
            self.length /= 2.0;
            new_edges.push(new_edge);
        }

        if self.network_state.state.contains(EDGE_DUPLICATE) {
            // differentiate
            let mut child_state = self.network_state.clone();
            child_state.state.set(EDGE_DIFFERENTIATION,
                                  !self.network_state.state.contains(EDGE_DIFFERENTIATION));

            let new_edge = Edge {
                src_node: self.dst_node,
                dst_node: self.src_node,
                length: self.length,
                network_state: child_state,
                type_count: self.type_count, /* XXX: start with 0 or with the count of the current edge? */
            };
            new_edges.push(new_edge);
        }

        if self.network_state.state.contains(EDGE_SWAP) {
            mem::swap(&mut self.dst_node, &mut self.src_node);
        }

        if self.network_state.state.contains(EDGE_GROW) {
            self.length *= 1.25;
        }

        if self.network_state.state.contains(EDGE_SHRINK) {
            self.length *= 0.75;
        }

        if self.network_state.state.contains(EDGE_TYPE) {
            self.type_count += 1;
        }
    }
}

#[derive(Debug)]
struct GraphBuilder {
    edges: Vec<Edge>,
    next_node_id: usize,
    network: GeneNetwork,
}

impl GraphBuilder {
    fn new(network: GeneNetwork, zygote: GeneNetworkState) -> GraphBuilder {
        let initial_edge = Edge {
            src_node: 0,
            dst_node: 1,
            length: 1.0,
            type_count: 0,
            network_state: zygote,
        };

        GraphBuilder {
            edges: vec![initial_edge],
            next_node_id: 2,
            network: network,
        }
    }

    // During the process, some edges will be added (split), some others will be modified.
    fn next(&mut self) {
        let mut new_edges = Vec::new();
        for edge in self.edges.iter_mut() {
            edge.develop(&self.network, &mut self.next_node_id, &mut new_edges);
        }
        println!("next_node_id: {}", self.next_node_id);
        println!("new edges: {:?}", new_edges);
        self.edges.extend(new_edges);
    }

    fn write_dot<W: Write>(&self, wr: &mut W) -> io::Result<()> {
        try!(writeln!(wr, "digraph artificial {{"));

        for edge in self.edges.iter() {
            try!(writeln!(wr,
                          "{} -> {} [weight={} label=\"{}\"]",
                          edge.src_node,
                          edge.dst_node,
                          edge.length,
                          edge.type_count));
        }

        try!(writeln!(wr, "}}"));


        Ok(())
    }
}


fn main() {
    use std::fs::File;
    // let genome = Genome::<Base4>::from_str("...11 _0320_23 <0101> T:0311 2...3 _1022_ 133 <0101> \
    // W:3213 121...")
    let mut rng = rand::thread_rng();

    let genome = Genome::<Base4>::random(&mut rng, 10 * 256);

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
                                           });

    println!("{:#?}", network);

    let mut zygote = network.new_state();
    zygote.state.set(0, true);
    zygote.state.set(1, true);

    let mut gb = GraphBuilder::new(network, zygote);
    println!("{:#?}", gb);

    for i in 0..3 {
        gb.next();
    }
    println!("{:#?}", gb);

    gb.write_dot(&mut File::create("example1.dot").unwrap());

    // XXX: Convert each edge into a node and then connect these new nodes.
}
