extern crate ndarray;
extern crate petgraph;
extern crate rand;

use std::collections::{HashMap, HashSet};
use std::ffi::c_void;
use std::fs;
use std::process::Command;
use std::thread::sleep;
use std::time::{Duration, Instant};

use fixedbitset::FixedBitSet;
use ndarray::Array2;
use petgraph::dot::{Config, Dot};
use petgraph::visit::GetAdjacencyMatrix;
use petgraph::{Graph, graph::UnGraph};
use rand::Rng;

use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use serde::{Deserialize, Serialize};
use serde_json;
use sundials_sys::*;

use crate::BinOpType::*;
use crate::Ex::*;
use my_cas::*;


fn main() {
    let eq1 = Ex::binop(
        Sub,
        Ex::binop(Mul, Ex::var("T"), Ex::var("x")),
        Ex::der(Ex::der(Ex::var("x"))),
    );

    let eq2 = Ex::binop(
        Sub,
        Ex::binop(
            Sub,
            Ex::binop(Mul, Ex::var("T"), Ex::var("y")),
            Ex::par("g"),
        ),
        Ex::der(Ex::der(Ex::var("y"))),
    );

    let eq3 = Ex::binop(
        Sub,
        Ex::binop(
            Add,
            Ex::binop(Pow, Ex::var("x"), Ex::constant(2.0)),
            Ex::binop(Pow, Ex::var("y"), Ex::constant(2.0)),
        ),
        Ex::binop(Pow, Ex::par("r"), Ex::constant(2.0)),
    );

    let sys = System {
        equations: vec![eq1, eq2, eq3],
        defaults: vec![],
        tspan: (0.0, 1.0),
    };

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
}
