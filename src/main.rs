extern crate ndarray;
extern crate petgraph;
extern crate rand;

use std::collections::{HashMap, HashSet};
use std::ffi::c_void;
use std::fs;
use std::future::pending;
use std::process::Command;
use std::thread::sleep;
use std::time::{Duration, Instant};

use fixedbitset::FixedBitSet;
use ndarray::Array2;
use petgraph::Undirected;
use petgraph::dot::{Config, Dot};
use petgraph::visit::GetAdjacencyMatrix;
use petgraph::{graph::UnGraph, Graph};
use rand::Rng;

use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use serde::{Deserialize, Serialize};
use serde_json;
use sundials_sys::*;

use crate::BinOpType::*;
// use crate::*;
use my_cas::Ex::*;
use my_cas::*;

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

    let eqs = vec![
        der(x.clone(), 1) - x_t.clone(),
        der(y.clone(), 1) - y_t.clone(),
        der(x_t.clone(), 1) - T.clone() * x.clone(),
        der(y_t.clone(), 1) - (T.clone() * y.clone() - g.clone()),
        pow(x.clone(), 2.0) + pow(y.clone(), 2.0) - pow(r, 2.0),
    ];
    System {
        equations: eqs,
        defaults: vec![],
        tspan: (0.0, 1.0),
    }
}

pub fn extract_diff_variables(expr: &Ex, variables: &mut HashSet<Ex>) {
    match expr {
        Ex::Var(name) => {
            variables.insert(Ex::Var(name.clone()));
        }
        Ex::Der(inner, _) => {
            variables.insert(expr.clone());
            extract_diff_variables(inner, variables);
        }
        Ex::BinaryOp(_, left, right) => {
            extract_diff_variables(left, variables);
            extract_diff_variables(right, variables);
        }
        _ => {}
    }
}

pub fn extract_diff_vars_system(system: &System) -> HashSet<Ex> {
    let mut variables = HashSet::new();
    for equation in &system.equations {
        extract_diff_variables(equation, &mut variables);
    }
    variables
}

fn extract_var_from_derivative(ex: &Ex) -> Option<&str> {
    if let Ex::Der(boxed_ex, _) = ex {
        if let Ex::Var(name) = &**boxed_ex {
            return Some(name);
        }
    }
    None
}

fn create_association_list(dvars_vec: &[Ex]) -> Vec<usize> {
    let mut association_list = vec![];

    for ex in dvars_vec.iter() {
        match ex {
            Ex::Var(name) => {
                if let Some(index) = dvars_vec
                    .iter()
                    .position(|other_ex| extract_var_from_derivative(other_ex) == Some(name))
                {
                    association_list.push(index);
                } else {
                    association_list.push(0);
                }
            }
            _ => {
                association_list.push(0);
            }
        }
    }

    association_list
}

// todo!();
// fn detect_subsets_to_be_differentiated(
//     initial_num_equations: usize,        // N
//     num_variables: usize,                // M
//     bipartite_graph: &Graph<String, (), Undirected>,    // Assumed type, you can replace with the actual type
//     variable_association_list: &[usize], // A
// ) {
//     // Step 1: Initialization
//     let mut assign = vec![0; num_variables]; // ASSIGN(j) = 0 for j = 1(1)M
//     let mut b = vec![0; initial_num_equations]; // B(i) = 0 for i = 1(1)N

//     // Step 2: Set N' = N
//     let mut n_prime = initial_num_equations;

//     // Step 3: Loop over N'
//     for k in 1..=n_prime {
//         // Step 3a: Set i = k
//         let i = k;

//         // Step 3b: Repeat (TODO: Implement the repeat logic)
//         // This is typically a while loop, but you need more specifics on the exit condition.
//         while true {
//             // TODO: Implement the inside of the repeat loop.
//         }
//     }
//     // TODO: Complete the rest of the algorithm
// }

fn main() {
    let sys = pend_sys();
    println!("{:#?}", sys);
    let graph = build_bipartite_graph(&sys);

    println!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
    let isb = petgraph::algo::is_bipartite_undirected(&graph, 0.into());
    println!("Is bipartite: {}", isb);
    let adj_mat = graph.adjacency_matrix();
    let array = convert_to_ndarray(&adj_mat, graph.node_count());
    let int_array: Array2<i32> = array.mapv(|x| if x { 1 } else { 0 });

    println!("{:#?}", int_array);

    let diff_idxs = sys.differential_indices();
    println!("{:?}", diff_idxs);

    let sys = order_1_pend_sys();
    sys.equations
        .iter()
        .for_each(|x| println!("{}", x.pretty_print()));
    let g = build_bipartite_graph(&sys);
    println!("sys \n{:?}", Dot::with_config(&g, &[Config::EdgeNoLabel]));

    let pend_vars = extract_variables_system(&sys);
    println!("{:?}", pend_vars);

    let pend_dvars = extract_diff_vars_system(&sys);
    let mut dvars_vec: Vec<_> = pend_dvars.into_iter().collect();
    dvars_vec.sort();
    dvars_vec
        .iter()
        .for_each(|x| println!("{}", x.pretty_print()));
    let vd = create_association_list(&dvars_vec);
    println!("vd:{:?}", vd);

    // vec![1,1,3].dedup();
    let e_nodes = vec![0, 3, 6, 8, 9];
    let v_nodes = g
        .node_indices()
        .filter(|x| !e_nodes.contains(&x.index()))
        .collect::<Vec<_>>();

    println!("{:?}", v_nodes);
}
