extern crate petgraph;

use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use ndarray;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::ffi::c_void;
use std::fs;
use std::process::Command;
use std::thread::sleep;
use std::time::{Duration, Instant};
use sundials_sys::*;

use petgraph::dot::{Config, Dot};
use petgraph::Graph;

use my_cas::*;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum NodeType {
    Equation(String),
    Variable(String),
}

fn construct_bipartite_graph(system: &System) -> Graph<NodeType, ()> {
    let mut graph = Graph::<NodeType, ()>::new();

    // Iterate over each equation
    for eq in &system.equations {
        let eq_node = graph.add_node(NodeType::Equation(eq.variable.clone()));

        // Get variables from the equation's expression
        let variables = get_variables_from_expression(&eq.expression);

        for var in variables {
            let var_node = graph.add_node(NodeType::Variable(var));

            // If this variable node isn't already connected to this equation, create an edge
            if !graph.contains_edge(eq_node, var_node) {
                graph.add_edge(eq_node, var_node, ());
            }
        }
    }

    graph
}

fn get_variables_from_expression(expr: &Ex) -> Vec<String> {
    match expr {
        Ex::Const(_) => vec![],
        Ex::Var(v) => vec![v.clone()],
        Ex::Par(_) => vec![],
        Ex::Mul(a, b) | Ex::Add(a, b) | Ex::Pow(a, b) => {
            let mut vars = get_variables_from_expression(a);
            vars.extend(get_variables_from_expression(b));
            vars
        }
        Ex::Der(e) => get_variables_from_expression(e),
    }
}

fn main() {
    // Example system (populate this)
    // let (system, u0, p) = lorenz_sys();
    let sys = pendulum_sys();
    let graph = construct_bipartite_graph(&sys);

    // Print out the graph in DOT format
    println!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
    let isb = petgraph::algo::is_bipartite_undirected(&graph, 0.into());
    println!("Is bipartite: {}", isb);
}
