#![feature(is_sorted)]

use std::collections::HashSet;
use std::fmt;

use crate::BinOpType::*;
use crate::Ex::*;
use colored::*;
use fixedbitset::FixedBitSet;
use ndarray::Array2;
use ordered_float::NotNan;
use petgraph::Undirected;
use petgraph::dot::{Config, Dot};
use petgraph::stable_graph::NodeIndex;
use petgraph::stable_graph::StableGraph;
use petgraph::stable_graph::StableUnGraph;
use petgraph::visit::GetAdjacencyMatrix;
use petgraph::{graph::UnGraph, Graph};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Ex {
    Const(NotNan<f64>),
    Var(String), // E.g., "x", "y" - usually for variables
    Par(String), // E.g., "a", "b" - for parameters
    BinOp(BinOpType, Box<Ex>, Box<Ex>),
    UnOp(UnOpType, Box<Ex>),
    Der(Box<Ex>, usize), // Represents differentiation wrt time with order
}

impl Ex {
    pub fn differential_index(&self) -> usize {
        match self {
            Ex::Const(_) | Ex::Var(_) | Ex::Par(_) => 0,
            Ex::BinOp(_, left, right) => {
                std::cmp::max(left.differential_index(), right.differential_index())
            }
            Ex::UnOp(_, operand) => operand.differential_index(),
            Ex::Der(operand, n) => n + operand.differential_index(),
        }
    }

    pub fn pretty_print(&self) -> String {
        match self {
            Ex::Const(val) => format!("{}", val).cyan().to_string(), // Colored cyan
            Ex::Var(name) => format!("{}", name).green().to_string(), // Colored green
            Ex::Par(name) => format!("{}", name).blue().to_string(), // Colored blue
            Ex::BinOp(op, left, right) => {
                let l = left.pretty_print();
                let r = right.pretty_print();
                format!("({} {} {})", l, op, r)
            }
            Ex::UnOp(op, operand) => {
                let op_str = operand.pretty_print();
                match op {
                    UnOpType::Sin => format!("sin({})", op_str),
                    // ... Add other unary operations as needed
                    _ => format!("{}({})", op, op_str), // Placeholder
                }
            }
            Ex::Der(operand, n) => {
                let op_str = operand.pretty_print();
                if *n == 1 {
                    format!("der({})", op_str)
                } else {
                    format!("der({}, {})", op_str, n)
                }
            }
        }
    }
}

pub fn c(val: f64) -> Ex {
    Ex::Const(NotNan::new(val).unwrap())
}

pub fn var(name: &str) -> Ex {
    Ex::Var(name.to_string())
}

pub fn par(name: &str) -> Ex {
    Ex::Par(name.to_string())
}

pub fn binop(op: BinOpType, lhs: Ex, rhs: Ex) -> Ex {
    Ex::BinOp(op, Box::new(lhs), Box::new(rhs))
}

pub fn unop(op: UnOpType, operand: Ex) -> Ex {
    Ex::UnOp(op, Box::new(operand))
}

pub fn der(expr: Ex, n: usize) -> Ex {
    Ex::Der(Box::new(expr), n)
}

// i dont like
pub fn pow(ex: Ex, exponent: f64) -> Ex {
    Ex::BinOp(BinOpType::Pow, Box::new(ex), Box::new(c(exponent)))
}

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
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

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum UnOpType {
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

macro_rules! define_op_overloads {
    ($($trait:ident, $method:ident, $op:path);* $(;)?) => {
        $(
            impl std::ops::$trait for Ex {
                type Output = Ex;

                fn $method(self, rhs: Ex) -> Ex {
                    binop($op, self, rhs)
                }
            }
        )*
    }
}

impl fmt::Display for BinOpType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            BinOpType::Add => write!(f, "+"),
            BinOpType::Sub => write!(f, "-"),
            BinOpType::Mul => write!(f, "*"),
            BinOpType::Div => write!(f, "/"),
            BinOpType::Pow => write!(f, "^"),
            // ... Add other binary operations as needed
            _ => write!(f, "?"), // Placeholder for unhandled operations
        }
    }
}

impl fmt::Display for UnOpType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            UnOpType::Sin => write!(f, "sin"),
            UnOpType::Cos => write!(f, "cos"),
            // ... Add other unary operations as needed
            _ => write!(f, "?"), // Placeholder for unhandled operations
        }
    }
}

// std::ops::
define_op_overloads! {
    Mul, mul, BinOpType::Mul;
    Add, add, BinOpType::Add;
    Sub, sub, BinOpType::Sub;
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

pub fn build_bipartite_graph(system: &System) -> (UnGraph<Ex, ()>, HashSet<NodeIndex>, HashSet<NodeIndex>) {
    let mut graph = UnGraph::<Ex, ()>::new_undirected();
    // add all enodes
    // then get_dvars_list, sort it, then add all vnodes
    

    let mut v_nodes = HashSet::new();
    let mut e_nodes = HashSet::new();

    // Add equations to the graph as nodes.
    for (index, eq) in system.equations.iter().enumerate() {
        // let eq_node = format!("Eq{}", index + 1);
        let eq_node_idx = graph.add_node(eq.clone());
        e_nodes.insert(eq_node_idx);

        // Extract variables from the equation.
        let mut variables = HashSet::new();
        extract_variables(&eq, &mut variables);

        for var in &variables {
            let var_node_idx = if let Some(idx) = graph.node_indices().find(|&n| graph[n] == Var(var.clone())) {
                idx
            } else {
                let vnode_idx = graph.add_node(Var(var.clone()));
                v_nodes.insert(vnode_idx);
                vnode_idx
            };
            graph.add_edge(eq_node_idx, var_node_idx, ());
        }
    }

    (graph, v_nodes, e_nodes)
}

pub fn build_bipartite_graph2_electric_boogaloo(system: &System) -> (StableUnGraph<Ex, ()>, HashSet<NodeIndex>, HashSet<NodeIndex>) {
    
    // StableUnGraph::with_capacity(u32::MAX, u32::MAX)
    let mut g: StableUnGraph<Ex, ()>= StableUnGraph::default();

    let mut v_nodes = HashSet::new();
    let mut e_nodes = HashSet::new();

    let mut dvars_list = extract_diff_vars_system(system).into_iter().collect::<Vec<_>>();
    dvars_list.sort(); // vnodes
    let eqs = &system.equations; // e nodes
    assert!(eqs.is_sorted());

    for (i, eq) in eqs.iter().enumerate() {
        let eq_node_idx = g.add_node(eq.clone());
        e_nodes.insert(eq_node_idx);
    }

    for (i, var) in dvars_list.iter().enumerate() {
        let var_node_idx = g.add_node(var.clone());
        v_nodes.insert(var_node_idx);
    }

    for (i, eq) in eqs.iter().enumerate() {
        let mut variables = HashSet::new();
        extract_diff_variables(&eq, &mut variables);

        for var in &variables {
            let var_node_idx = g.node_indices().find(|&n| g[n] == var.clone()).unwrap();
            g.add_edge((i as u32).into(), var_node_idx, ());
        }
    }
    
    (g, v_nodes, e_nodes)
}


// A recursive function to extract variables from an expression.
pub fn extract_variables(expr: &Ex, variables: &mut HashSet<String>) {
    match expr {
        Var(name) => {
            variables.insert(name.clone());
        }
        UnOp(_, inner) => {
            extract_variables(inner, variables);
        }
        BinOp(_, left, right) => {
            extract_variables(left, variables);
            extract_variables(right, variables);
        }
        Der(inner, _) => {
            extract_variables(inner, variables);
        }
        _ => {}
    }
}

pub fn extract_variables_system(system: &System) -> HashSet<String> {
    let mut variables = HashSet::new();

    for equation in &system.equations {
        extract_variables(equation, &mut variables);
    }

    variables
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

pub fn extract_diff_variables(expr: &Ex, variables: &mut HashSet<Ex>) {
    match expr {
        Ex::Var(name) => {
            variables.insert(Ex::Var(name.clone()));
        }
        Ex::Der(inner, _) => {
            variables.insert(expr.clone());
            // extract_diff_variables(inner, variables);
        }
        Ex::UnOp(_, inner) => {
            variables.insert(expr.clone());
            extract_diff_variables(inner, variables);
        }
        Ex::BinOp(_, left, right) => {
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

pub fn extract_var_from_derivative(ex: &Ex) -> Option<&str> {
    if let Ex::Der(boxed_ex, _) = ex {
        if let Ex::Var(name) = &**boxed_ex {
            return Some(name);
        }
    }
    None
}
