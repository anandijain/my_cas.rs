#![feature(is_sorted)]

extern crate ndarray;
extern crate petgraph;
extern crate rand;

use ndarray::Array2;
use petgraph::dot::{Config, Dot};
use petgraph::graph::Node;
use petgraph::visit::GetAdjacencyMatrix;
use petgraph::Graph;
use petgraph::Undirected;
use std::collections::{HashMap, HashSet};
use std::path;

use my_cas::Ex::*;
use my_cas::*;
use petgraph::stable_graph::{NodeIndex, StableDiGraph, StableGraph, StableUnGraph};

type G_TYPE = StableGraph<Ex, (), Undirected>;

pub fn pend_sys() -> System {
    let x = var("x");
    let y = var("y");
    let T = var("T");
    let r = par("r");
    let g = par("g");

    let eq1 = T.clone() * x.clone() - der(x.clone(), 2);
    let eq2 = T.clone() * y.clone() - g - der(y.clone(), 2);
    let eq3 = pow(x.clone(), 2.0) + pow(y.clone(), 2.0) - pow(r, 2.0);

    let sys = System {
        equations: vec![eq1, eq2, eq3],
        defaults: vec![],
        tspan: (0.0, 1.0),
    };
    sys
}

pub fn order_1_pend_sys() -> System {
    let x = var("x");
    let y = var("y");
    let T = var("T");
    let x_t = var("x_t");
    let y_t = var("y_t");
    let r = par("r");
    let g = par("g");

    let mut eqs = vec![
        der(x.clone(), 1) - x_t.clone(),
        der(y.clone(), 1) - y_t.clone(),
        der(x_t.clone(), 1) - T.clone() * x.clone(),
        der(y_t.clone(), 1) - (T.clone() * y.clone() - g.clone()),
        pow(x.clone(), 2.0) + pow(y.clone(), 2.0) - pow(r, 2.0),
    ];
    eqs.sort();

    System {
        equations: eqs,
        defaults: vec![],
        tspan: (0.0, 1.0),
    }
}

fn create_association_list(dvars_vec: &[Ex]) -> HashMap<NodeIndex, NodeIndex> {
    let mut association_list = HashMap::new();

    for (i, ex) in dvars_vec.iter().enumerate() {
        match ex {
            Ex::Var(name) => {
                if let Some(index) = dvars_vec
                    .iter()
                    .position(|other_ex| extract_var_from_derivative(other_ex) == Some(name))
                {
                    association_list.insert(NodeIndex::new(i), NodeIndex::new(index));
                }
            }
            _ => {}
        }
    }

    association_list
}

// todo fix the u32 / nodeindex weirddness
// assign is var -> eq

fn augment_path(
    graph: &G_TYPE,
    i: NodeIndex,
    path_found: &mut bool,
    v_nodes: &HashSet<NodeIndex>, // dont need to be mutable here
    e_nodes: &HashSet<NodeIndex>,
    assign: &mut HashMap<NodeIndex, NodeIndex>,
    colored: &mut HashSet<NodeIndex>,
) {
    // (1) Colour i
    colored.insert(i);

    // (2) Check V-nodes
    for &j in v_nodes.iter() {
        println!("j: {:?}", j);
        if graph.contains_edge(i.into(), j.into()) && assign.get(&j) == None {
            // (2a) Set PATHFOUND TRUE
            *path_found = true;
            // (2b) Set ASSIGN (j) to i (assuming this is the correct assignment)
            assign.insert(j, i);
            // (2c) Return
            return;
        }
    }

    // (3) Check all uncolored j nodes
    for &j in e_nodes.iter() {
        // graph.edges(j)
        if graph.contains_edge(i.into(), j.into()) && !colored.contains(&j) {
            // (3a) Colour j
            colored.insert(j);
            // (3b) Set k to ASSIGN(j)
            let k = *assign.get(&j).unwrap(); //.unwrap_or(&NodeIndex::new(0)); // sus of default.
                                              // (3c) Recursive call
            augment_path(graph, k, path_found, v_nodes, e_nodes, assign, colored);
            // (3d) If PATHFOUND then
            if *path_found {
                // (3d-1) Set ASSIGN (j) to i (again assuming this is the correct assignment)
                assign.insert(j, i);
                // (3d-2) Return
                return;
            }
        }
    }

    // (4) Return
    return;
}

fn detect_subsets_to_be_differentiated(
    mut num_equations: usize, // N
    mut num_variables: usize, // M
    g: &mut G_TYPE,
    // variable_association_list: &[usize], // A. maybe should be HashMap<NodeIndex, NodeIndex> st all keys and values are Vnodes
    var_diff_map: &mut HashMap<NodeIndex, NodeIndex>, // A. maybe should be HashMap<NodeIndex, NodeIndex> st all keys and values are Vnodes
    v_nodes: &mut HashSet<NodeIndex>,
    e_nodes: &mut HashSet<NodeIndex>,
) {
    // Step 1: Initialization
    // let mut assign = vec![0; num_variables]; // ASSIGN(j) = 0 for j = 1(1)M
    let mut assign: HashMap<NodeIndex, NodeIndex> = HashMap::new(); // Vnode to Enode
    let mut b: HashMap<NodeIndex, NodeIndex> = HashMap::new(); // Enode to Enode
    let mut colored: HashSet<NodeIndex> = HashSet::new(); // COLOURED = {}
    let mut path_found = false;

    // Step 2: Set N' = N
    let mut n_prime = num_equations;

    // Step 3: Loop over N'
    for k in 1..=n_prime {
        // Step 3a: Set i = k
        let mut i = k;

        // Step 3b: Repeat (TODO: Implement the repeat logic)
        // This is typically a while loop, but you need more specifics on the exit condition.
        loop {
            // this might need to also modify num variables
            three_b_1(g, &*var_diff_map, v_nodes); // MODIFIES GRAPH.
                                                   // num_variables = v_nodes.len();
                                                   // println!("{:?}", Dot::new(g));
            path_found = false;
            colored.clear();

            augment_path(
                &*g,
                NodeIndex::new(i),
                &mut path_found,
                v_nodes,
                e_nodes,
                &mut assign,
                &mut colored,
            );
            if path_found {
                break;
            } else {
                //  i) For every coloured V-node j do
                for v_node in v_nodes.intersection(&colored) {
                    // SetM M+1
                    num_variables += 1;
                    // Create new V-node M
                    let newvar = Var(format!("var_{}", num_variables)); // Suseyyy
                    let new_vn = g.add_node(newvar);
                    // SetA(j) M
                    var_diff_map.insert(*v_node, new_vn);
                }

                // i)ForeverycolouredE-node do
                for enode in e_nodes.intersection(&colored) {
                    // SetN N+1
                    num_equations += 1;
                    // Create new E-node N

                    let new_enode = g.add_node(Var("TODO".to_owned())); // node weight needs to be the diff'd eq

                    // Create edges from E-node N to all V-nodes j and A(j) such that edge (l-j) exists.
                    // let es = g.edges(*enode).;
                    let ns = g.neighbors(*enode).into_iter().collect::<Vec<_>>();
                    for n_ in ns { // these are the `j`s 
                        assert!(g.contains_edge(*enode, n_));
                        let e_ = g.find_edge(*enode, n_).unwrap();
                        
                        g.add_edge(new_enode, n_, ());
                        g.add_edge(new_enode, *var_diff_map.get(&n_).unwrap(), ());
                    }

                    // Set b(l)= N
                    b.insert(*enode, NodeIndex::new(num_equations));
                }
                // (iii) For every coloured V-node j
                for v_node in v_nodes.intersection(&colored) {
                    // set ASSIGN (A(j))= B(ASSIGN (j))
                    let a_j = var_diff_map.get(&v_node).unwrap();
                    let assign_j = assign.get(&v_node).unwrap();

                    assign.insert(*a_j, *assign_j);
                }

                i = b
                    .get(&NodeIndex::new(i.try_into().unwrap()))
                    .unwrap()
                    .index();
            }
        }
    }
}

// vnode to vnode
pub fn build_a(
    g: &G_TYPE,
    v_nodes: &HashSet<NodeIndex>,
    // e_nodes: &HashSet<NodeIndex>,
) -> HashMap<NodeIndex, NodeIndex> {
    let mut a = HashMap::new();
    // let es = g.edge_indices();
    for v in v_nodes {
        let w = g.node_weight(*v).unwrap();
        match w {
            Var(x) => {
                for v2 in v_nodes {
                    let w2 = g.node_weight(*v2).unwrap();
                    match w2 {
                        Der(inner, _) => match *inner.clone() {
                            Var(var) => {
                                if x == &var {
                                    a.insert(*v, *v2);
                                }
                            }
                            _ => {}
                        },
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
    a
}

/// Delete all V-nodes with A(. )+0 and all their incident edges from the graph
pub fn three_b_1(
    g: &mut G_TYPE,
    a: &HashMap<NodeIndex, NodeIndex>,
    v_nodes: &mut HashSet<NodeIndex>,
) {
    for (k, _) in a {
        v_nodes.remove(k);
        g.remove_node(*k);
    }
}
fn main() {
    let sys = pend_sys();
    // println!("{:#?}", sys);
    let (mut graph, mut _v_nodes, mut _e_nodes) = build_bipartite_graph(&sys);

    println!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
    let isb = petgraph::algo::is_bipartite_undirected(&graph, 0.into());
    println!("Is bipartite: {}", isb);
    let adj_mat = graph.adjacency_matrix();
    let array = convert_to_ndarray(&adj_mat, graph.node_count());
    let int_array: Array2<i32> = array.mapv(|x| if x { 1 } else { 0 });

    println!("{:#?}", int_array);

    let diff_idxs = sys.differential_indices();
    println!("differential indices: {:?}", diff_idxs);

    let sys = order_1_pend_sys(); // already calls sort on eqs
    assert!(sys.equations.is_sorted());
    let dvs = extract_diff_vars_system(&sys);
    assert_eq!(dvs.len(), 9);
    println!("X = {:?}", dvs);

    sys.equations
        .iter()
        .for_each(|x| println!("{}", x.pretty_print()));

    let (mut g, mut v_nodes, mut e_nodes) = build_bipartite_graph2_electric_boogaloo(&sys);
    println!("{:?}", Dot::with_config(&g, &[Config::EdgeNoLabel]));
    println!("nv:{:?}, ne:{:?}", g.node_count(), g.edge_count());
    println!("v_nodes:{:?} \ne_nodes:{:?}", v_nodes, e_nodes);

    let a = build_a(&g, &v_nodes);
    println!("A = {:?}", a);
    for (k, v) in &a {
        // println!("k: {:?}, v: {:?}", k, v);
        let kw = g.node_weight(*k).unwrap();
        let vw = g.node_weight(*v).unwrap();
        println!("k: {:?}, v: {:?}", kw, vw);
    }

    let mut g2 = g.clone();
    let mut v_nodes2 = v_nodes.clone();

    three_b_1(&mut g2, &a, &mut v_nodes2);
    println!("{:?}", Dot::with_config(&g2, &[Config::EdgeNoLabel]));
}

fn stable_graph_example() {
    let mut graph: StableDiGraph<&str, &str> = StableDiGraph::new();

    // Add some nodes and edges to the graph
    let node_a = graph.add_node("A");
    let node_b = graph.add_node("B");
    let node_c = graph.add_node("C");
    let edge_1 = graph.add_edge(node_a, node_b, "edge1");
    let edge_2 = graph.add_edge(node_b, node_c, "edge2");

    // Print the graph
    println!("Initial graph:");
    println!("{:?}", Dot::new(&graph));

    // Remove node B and edge (A, B)
    graph.remove_node(node_b);
    // graph.remove_edge(edge_1);

    // Print the graph after removals
    println!("\nGraph after removing node B and edge (A, B):");
    println!("{:?}", Dot::new(&graph));

    // Add node B and edge (A, B) back
    let node_b_new = graph.add_node("B");
    let node_d_new = graph.add_node("D");
    let edge_1_new = graph.add_edge(node_a, node_b_new, "edge1");

    // Print the graph after additions
    println!("\nGraph after adding node B and edge (A, B) back:");
    println!("{:?}", Dot::new(&graph));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pend_nv() {
        let sys = order_1_pend_sys();
        let (g, v_nodes, e_nodes) = build_bipartite_graph2_electric_boogaloo(&sys);
        assert_eq!(g.node_count(), 14);
        assert_eq!(g.edge_count(), 12);
        assert_eq!(v_nodes.len(), 9);
        assert_eq!(e_nodes.len(), 5);
    }
}
