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
use petgraph::{graph::UnGraph, Graph};
use rand::Rng;

use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use serde::{Deserialize, Serialize};
use serde_json;
use sundials_sys::*;

use crate::BinOpType::*;
// use crate::*;
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

pub fn ode_order_lowering(expr: Ex, equations: &mut Vec<Ex>) -> Ex {
    match expr {
        Ex::Der(inner, order) => {
            if order > 1 {
                let base_var = match &*inner {
                    Ex::Var(v) => v.clone(),
                    _ => panic!("Expected a variable!"),
                };

                let mut prev_var = base_var.clone();
                for _ in 1..order {
                    let aux_var = format!("{}_t", prev_var);
                    let diff_equation = binop(
                        BinOpType::Sub,
                        Ex::Der(Box::new(Ex::Var(prev_var.clone())), 1),
                        Ex::Var(aux_var.clone()),
                    );
                    equations.push(diff_equation);
                    prev_var = aux_var;
                }

                return Ex::Var(prev_var);
            } else {
                Ex::Der(Box::new(ode_order_lowering(*inner, equations)), 1)
            }
        }
        Ex::BinaryOp(op, left, right) => {
            let lowered_left = ode_order_lowering(*left, equations);
            let lowered_right = ode_order_lowering(*right, equations);
            Ex::BinaryOp(op, Box::new(lowered_left), Box::new(lowered_right))
        }
        Ex::UnaryOp(op, operand) => {
            let lowered_operand = ode_order_lowering(*operand, equations);
            Ex::UnaryOp(op, Box::new(lowered_operand))
        }
        _ => expr,
    }
}

pub fn lower_equations(equations: Vec<Ex>) -> Vec<Ex> {
    let mut lowered_equations = equations.clone();

    for equation in equations.iter() {
        let lowered_equation = ode_order_lowering(equation.clone(), &mut lowered_equations);
        if let Some(index) = lowered_equations.iter().position(|x| *x == *equation) {
            lowered_equations[index] = lowered_equation;
        } else {
            lowered_equations.push(lowered_equation);
        }
    }

    lowered_equations
}

pub fn lower_system(system: System) -> System {
    System {
        equations: lower_equations(system.equations),
        defaults: system.defaults,
        tspan: system.tspan,
    }
}

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

    let sys = pend_sys();
    let lowered_sys = lower_system(sys);

    for equation in lowered_sys.equations {
        println!("{:?}", equation);
    }

    let ex = der(var("x"), 3) - c(1.0);
    println!("{:?}", ex);
    let exs = lower_equations(vec![ex]);
    println!("{:?}", exs);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lower_system() {
        let sys = pend_sys();
        let lowered_sys = lower_system(sys);

        let expected_eqs = vec![
            binop(Sub, binop(Mul, var("T"), var("x")), var("x_t")),
            binop(
                Sub,
                binop(Sub, binop(Mul, var("T"), var("y")), par("g")),
                var("y_t"),
            ),
            binop(
                Sub,
                binop(
                    Add,
                    binop(Pow, var("x"), c(2.0)),
                    binop(Pow, var("y"), c(2.0)),
                ),
                binop(Pow, par("r"), c(2.0)),
            ),
            binop(Sub, Ex::Der(Box::new(var("x")), 1), var("x_t")),
            binop(Sub, Ex::Der(Box::new(var("y")), 1), var("y_t")),
        ];

        for (given, expected) in lowered_sys.equations.iter().zip(expected_eqs.iter()) {
            assert_eq!(given, expected);
        }
    }
}
