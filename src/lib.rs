use std::collections::HashSet;

use crate::BinOpType::*;
use crate::Ex::*;
use fixedbitset::FixedBitSet;
use ndarray::Array2;
use petgraph::dot::{Config, Dot};
use petgraph::visit::GetAdjacencyMatrix;
use petgraph::{graph::UnGraph, Graph};
use serde::{Deserialize, Serialize};

type Node = String;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub enum Ex {
    Const(f64),
    Var(String), // E.g., "x", "y" - usually for variables
    Par(String), // E.g., "a", "b" - for parameters
    BinaryOp(BinOpType, Box<Ex>, Box<Ex>),
    UnaryOp(UnaryOpType, Box<Ex>),
    Der(Box<Ex>), // Represents differentiation wrt time
}

impl Ex {
    pub fn constant(val: f64) -> Self {
        Ex::Const(val)
    }

    pub fn var(name: &str) -> Self {
        Ex::Var(name.to_string())
    }

    pub fn par(name: &str) -> Self {
        Ex::Par(name.to_string())
    }

    pub fn binop(op: BinOpType, lhs: Ex, rhs: Ex) -> Self {
        Ex::BinaryOp(op, Box::new(lhs), Box::new(rhs))
    }

    pub fn unop(op: UnaryOpType, operand: Ex) -> Self {
        Ex::UnaryOp(op, Box::new(operand))
    }

    pub fn der(expr: Ex) -> Self {
        Ex::Der(Box::new(expr))
    }

    pub fn differential_index(&self) -> usize {
        match self {
            Ex::Const(_) | Ex::Var(_) | Ex::Par(_) => 0,
            Ex::BinaryOp(_, left, right) => {
                std::cmp::max(left.differential_index(), right.differential_index())
            }
            Ex::UnaryOp(_, operand) => operand.differential_index(),
            Ex::Der(operand) => 1 + operand.differential_index(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub enum BinOpType {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Mod, // Modulus
    Log, // Logarithm base left of right: e.g., log_2(8)
    Max, // Maximum of two numbers
    Min, // Minimum of two numbers
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub enum UnaryOpType {
    Sin,
    Cos,
    Tan,
    ASin,  // Arcsine
    ACos,  // Arccosine
    ATan,  // Arctangent
    Sinh,  // Hyperbolic Sine
    Cosh,  // Hyperbolic Cosine
    Tanh,  // Hyperbolic Tangent
    Exp,   // Exponential function
    Ln,    // Natural Logarithm
    Abs,   // Absolute Value
    Sqrt,  // Square Root
    Ceil,  // Ceiling function
    Floor, // Floor function
    Neg,   // Negation
    Fact,  // Factorial
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct System {
    pub equations: Vec<Ex>, // rhs - lhs
    pub defaults: Vec<(String, f64)>,
    pub tspan: (f64, f64),
}

impl System {
    /// Computes the differential index for each equation in the system
    pub fn differential_indices(&self) -> Vec<usize> {
        self.equations
            .iter()
            .map(|eq| eq.differential_index())
            .collect()
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Equation {
    pub lhs: Ex,
    pub rhs: Ex,
}

pub fn build_bipartite_graph(system: &System) -> UnGraph<Node, ()> {
    let mut graph = UnGraph::<Node, ()>::new_undirected();

    // Add equations to the graph as nodes.
    for (index, eq) in system.equations.iter().enumerate() {
        let eq_node = format!("Eq{}", index + 1);
        let eq_node_idx = graph.add_node(eq_node);

        // Extract variables from the equation.
        let mut variables = HashSet::new();
        extract_variables(&eq, &mut variables);

        for var in &variables {
            let var_node_idx = if let Some(idx) = graph.node_indices().find(|&n| graph[n] == *var) {
                idx
            } else {
                graph.add_node(var.clone())
            };
            graph.add_edge(eq_node_idx, var_node_idx, ());
        }
    }

    graph
}

// A recursive function to extract variables from an expression.
pub fn extract_variables(expr: &Ex, variables: &mut HashSet<String>) {
    match expr {
        Var(name) => {
            variables.insert(name.clone());
        }
        BinaryOp(_, left, right) => {
            extract_variables(left, variables);
            extract_variables(right, variables);
        }
        Der(inner) => {
            extract_variables(inner, variables);
        }
        _ => {}
    }
}

pub fn convert_to_ndarray(fbs: &FixedBitSet, n: usize) -> Array2<bool> {
    let mut array = Array2::from_elem((n, n), false);

    for i in 0..n {
        for j in 0..n {
            let index = i * n + j;
            if fbs[index] {
                array[[i, j]] = true;
            }
        }
    }

    array
}
