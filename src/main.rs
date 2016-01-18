#![feature(slice_patterns)]
#![allow(dead_code)]

extern crate artificial_genome;
extern crate rand;
extern crate fixedbitset;

use artificial_genome::{Genome, ProteinRegulator, GeneNetwork, GeneNetworkState};
use artificial_genome::base4::{Base4, B0, B1};
use std::mem;
use std::io::{self, Write};
use std::collections::{BTreeMap, BTreeSet};
use fixedbitset::FixedBitSet;

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

const RESIZE_FACTOR: f32 = 0.25;

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
            self.length += RESIZE_FACTOR * self.length;
        }

        if self.network_state.state.contains(EDGE_SHRINK) {
            self.length -= RESIZE_FACTOR * self.length;
        }

        if self.network_state.state.contains(EDGE_TYPE) {
            self.type_count += 1;
        }
    }
}

struct NodeGraph {
    nodes: Vec<NodeInfo>,
    edges: BTreeSet<(usize, usize)>,
}

struct NodeInfo {
    length: f32,
    type_count: usize,
}

#[derive(Debug)]
struct FoundPath {
    target_node: usize,
    length: f32,
}

// Find all paths from current_node to connected nodes.
fn path_finder(current_node: usize,
               current_length: f32,
               neighborhood: &Vec<Vec<usize>>,
               lengths: &Vec<f32>,
               targets: &BTreeSet<usize>,
               visited: &mut FixedBitSet,
               paths: &mut Vec<FoundPath>) {
    for &ni in neighborhood[current_node].iter() {
        if targets.contains(&ni) {
            // we found a target
            paths.push(FoundPath {
                target_node: ni,
                length: current_length,
            });
        } else if !visited.contains(ni) {
            // the node ``ni`` is not a target node and we haven't
            // visited yet.

            // we haven't visited ```ni``` yet. recurse down
            visited.set(ni, true);
            path_finder(ni,
                        current_length + lengths[ni],
                        neighborhood,
                        lengths,
                        targets,
                        visited,
                        paths);
            visited.set(ni, false);
        }
    }
}

#[derive(Debug)]
struct StructuredNode {
    length: f32,
    type_count: usize,
    connections: Vec<FoundPath>,
}

// There are two types of nodes. Processing nodes (e.g. Neuron) or
// connecting nodes (e.g. Synapses). Keep the processing nodes as
// nodes in the graph, while the connecting nodes become edges.
#[derive(Debug)]
struct StructuredGraph {
    nodes: BTreeMap<usize, StructuredNode>,
}

impl NodeGraph {
    fn write_dot<W: Write>(&self, wr: &mut W) -> io::Result<()> {
        try!(writeln!(wr, "digraph artificial {{"));

        // the edges are nodes in this graph.
        for (i, node) in self.nodes.iter().enumerate() {
            try!(writeln!(wr,
                          "{} [label=\"{}:{:.2}:{}\"]",
                          i,
                          i,
                          node.length,
                          node.type_count));
        }

        // now connect them
        for &(src_edge, dst_edge) in self.edges.iter() {
            try!(writeln!(wr, "{} -> {}", src_edge, dst_edge));
        }

        try!(writeln!(wr, "}}"));

        Ok(())
    }

    fn into_structured_graph(&self, type_processing: usize) -> StructuredGraph {
        // determine max type_count.
        // lets say every type_count >= 3 is a neuron for now.
        // everything else is a synapse
        // If two synapses are directly connected to each other,
        // treat them as one synapse with the summation of it's lengths.
        // XXX: remove synapses which are too long.
        // If two neurons are next to each other, what to do?
        // a) Treat them as one big neuron
        // b) Connect them with a minimum synapse?
        // We do b) for now, because it's easier.

        let mut neighbors: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for &(src, dst) in self.edges.iter() {
            neighbors.entry(src).or_insert_with(|| BTreeSet::new()).insert(dst);
        }

        let mut processing_nodes: BTreeSet<usize> = BTreeSet::new();

        for (i, node) in self.nodes.iter().enumerate() {
            if node.type_count >= type_processing {
                processing_nodes.insert(i);
            }
        }

        let neighborhood: Vec<Vec<usize>> = self.nodes
                                                .iter()
                                                .enumerate()
                                                .map(|(i, _)| {
                                                    neighbors.get(&i)
                                                             .map(|v| v.iter().cloned().collect())
                                                             .unwrap_or(Vec::new())
                                                })
                                                .collect();

        let lengths: Vec<f32> = self.nodes
                                    .iter()
                                    .map(|n| n.length)
                                    .collect();

        // now follow all possible directed paths from each processing_node, until
        // we arrive at another processing_node.
        // (Or repetetively merge consecutive elements which have same type).

        let mut visited = FixedBitSet::with_capacity(self.nodes.len());

        let mut g = StructuredGraph { nodes: BTreeMap::new() };

        for &pnode in processing_nodes.iter() {
            // find all processing elements connected either directly or through intermediate
            // synapses.
            let mut snode = StructuredNode {
                length: self.nodes[pnode].length,
                type_count: self.nodes[pnode].type_count,
                connections: Vec::new(),
            };

            visited.clear();
            path_finder(pnode,
                        0.0,
                        &neighborhood,
                        &lengths,
                        &processing_nodes,
                        &mut visited,
                        &mut snode.connections);

            g.nodes.insert(pnode, snode);
        }

        g
    }
}

impl StructuredGraph {
    fn write_dot<W: Write>(&self, wr: &mut W) -> io::Result<()> {
        try!(writeln!(wr, "digraph artificial {{"));

        for (&i, node) in self.nodes.iter() {
            if node.connections.is_empty() {
                // XXX
                continue;
            }
            try!(writeln!(wr,
                          "{} [label=\"{} {:.2} {}\"]",
                          i,
                          i,
                          node.length,
                          node.type_count));
        }

        // now connect them
        for (src, node) in self.nodes.iter() {
            for dst in node.connections.iter() {
                if dst.length == 0.0 {
                    // XXX
                    continue;
                }
                try!(writeln!(wr,
                              "{} -> {} [weight={} label=\"{}\"]",
                              src,
                              dst.target_node,
                              dst.length,
                              dst.length));
            }
        }

        try!(writeln!(wr, "}}"));

        Ok(())
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

    // The result of the GraphBuilder is a graph where every edge represents an element (either a
    // Neuron or a Synapse), while the nodes represent connection points (they carry purely
    // structural information). Transform this "edged" graph into a graph where the edges become
    // nodes.
    fn into_node_graph(&self) -> NodeGraph {
        // only keep one directed edge between each pair of nodes.
        // keep the edge with highest type_count.

        let mut e: Vec<Edge> = Vec::new();
        {
            // stores (src_node, dst_node) => edge idx
            let mut selected_edges: BTreeMap<(usize, usize), usize> = BTreeMap::new();
            for (i, edge) in self.edges.iter().enumerate() {
                let idx_ref = selected_edges.entry((edge.src_node, edge.dst_node)).or_insert(i);
                let old = &self.edges[*idx_ref];
                if edge.type_count > old.type_count {
                    *idx_ref = i;
                }
            }
            for (_, &edge_idx) in selected_edges.iter() {
                e.push(self.edges[edge_idx].clone());
            }
        }

        let mut node_graph = NodeGraph {
            nodes: Vec::with_capacity(e.len()),
            edges: BTreeSet::new(),
        };

        let mut node_out_edges: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

        // every edge becomes a node.
        // note that ```i``` matches the indices  we use for ```new_nodes```.
        for (i, edge) in e.iter().enumerate() {
            let j = node_graph.nodes.len();
            node_graph.nodes.push(NodeInfo {
                length: edge.length,
                type_count: edge.type_count,
            });
            assert!(i == j);
            node_out_edges.entry(edge.src_node).or_insert(Vec::new()).push(i);
        }
        assert!(node_graph.nodes.len() == e.len());

        for (i, edge) in e.iter().enumerate() {
            // connect i with all outgoing edges of dst_node

            if let Some(list) = node_out_edges.get(&edge.dst_node) {
                for &out_i in list.iter() {
                    node_graph.edges.insert((i, out_i));
                }
            }
        }

        node_graph
    }

    fn write_dot<W: Write>(&self, wr: &mut W) -> io::Result<()> {
        try!(writeln!(wr, "digraph artificial {{"));

        for edge in self.edges.iter() {
            try!(writeln!(wr,
                          "{} -> {} [weight={} label=\"{} / {}\"]",
                          edge.src_node,
                          edge.dst_node,
                          edge.length,
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
    // zygote.state.set(1, true);

    let mut gb = GraphBuilder::new(network, zygote);
    println!("{:#?}", gb);

    for _ in 0..5 {
        gb.next();
    }
    println!("{:#?}", gb);

    gb.write_dot(&mut File::create("example1.dot").unwrap()).unwrap();

    let node_graph = gb.into_node_graph();
    node_graph.write_dot(&mut File::create("example1_node.dot").unwrap()).unwrap();
    let g = node_graph.into_structured_graph(1);
    println!("g: {:?}", g);
    g.write_dot(&mut File::create("example1_struct.dot").unwrap()).unwrap();
}
