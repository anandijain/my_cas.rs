use cvode_wrap::{self, AbsTolerance, LinearMultistepMethod, Realtype, RhsResult, SolverNoSensi};
use diffeq;
use num_traits::{Float, Num};
use std::collections::{BTreeSet, HashSet};
use sundials_sys;
// std::ops
// A symbolic variable with potential dependencies on independent variables
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Var {
    name: String,
    independent_vars: BTreeSet<String>,
}

// A generic enum to represent mathematical Exprs, parameterized over integer and float types
#[derive(Debug, Clone)]
enum Expr<I: Num + Clone, R: Float + Clone> {
    Variable(Var),
    Int(I),
    Real(R),
    Pi, // Representation for π
    E,  // Representation for e
    Add(Box<Expr<I, R>>, Box<Expr<I, R>>),
    Mul(Box<Expr<I, R>>, Box<Expr<I, R>>),
    Div(Box<Expr<I, R>>, Box<Expr<I, R>>),
    Pow(Box<Expr<I, R>>, Box<Expr<I, R>>),
    Der(Box<Expr<I, R>>, Box<Var>, u32),
}

#[derive(Debug, Clone)]
struct Equation<I: Num + Clone, R: Float + Clone> {
    lhs: Expr<I, R>,
    rhs: Expr<I, R>,
}

/// do we want to make a new, so that we can ensure that the variables in the equation exactly the same as the dvs
#[derive(Debug, Clone)]
struct ODESystem<I: Num + Clone, R: Float + Clone> {
    equations: Vec<Equation<I, R>>,
    independent_vars: BTreeSet<Var>,
    dependent_vars: BTreeSet<Var>,
    defaults : HashMap<String, R>,
    tspan: (R, R),

}

impl<I: Num + Clone + Into<f64>, R: Float + Clone> ODESystem<I, R> {
    pub fn evaluate_expr(
        &self,
        expr: &Expr<I, R>,
        t: Realtype,
        y: &[Realtype],
        p: &[Realtype],
    ) -> R {
        match expr {
            Expr::Variable(var) => {
                if var.name == "t" {
                    return R::from(t).unwrap();
                }
                let pos = self
                    .dependent_vars
                    .iter()
                    .position(|v| *v == *var)
                    .expect("Variable not found in dependent vars");
                return R::from(y[pos]).unwrap();
            }
            Expr::Int(i) => R::from(i.clone().into()).unwrap(),
            Expr::Real(r) => r.clone(),
            Expr::Pi => R::from(std::f64::consts::PI).unwrap(),
            Expr::E => R::from(std::f64::consts::E).unwrap(),
            Expr::Add(left, right) => {
                self.evaluate_expr(left, t, y, p) + self.evaluate_expr(right, t, y, p)
            }
            Expr::Mul(left, right) => {
                self.evaluate_expr(left, t, y, p) * self.evaluate_expr(right, t, y, p)
            }
            Expr::Div(left, right) => {
                self.evaluate_expr(left, t, y, p) / self.evaluate_expr(right, t, y, p)
            }
            Expr::Pow(base, exp) => self
                .evaluate_expr(base, t, y, p)
                .powf(self.evaluate_expr(exp, t, y, p)),
            // For the Derivative, we just evaluate the inner expression for now.
            // Actual differentiation would be more involved.
            Expr::Der(expr, _, _) => self.evaluate_expr(expr, t, y, p),
        }
    }
}

impl<I: Num + Clone + Into<f64>, R: Float + Clone> ODESystem<I, R> {
    pub fn to_sundials_f(
        &self,
    ) -> OdeFn<impl Fn(Realtype, &[Realtype], &mut [Realtype], &[Realtype]) -> RhsResult> {
        let equations = self.equations.clone();
        let func = move |t: Realtype, y: &[Realtype], ydot: &mut [Realtype], p: &[Realtype]| {
            for (idx, eq) in equations.iter().enumerate() {
                // Use a separate evaluate function that doesn't rely on `self`
                ydot[idx] = self.evaluate_expr(&eq.rhs, t, y, p).to_f64().unwrap();
            }
            RhsResult::Ok
        };
        OdeFn { func }
    }
}

// impl<I: Num + Clone, R: Float + Clone> ODESystem<I, R> {
//     // to_diffeq_f returns a function that can be used by the diffeq crate
//     // to_iip_f returns in place function handle (in my own handle format) Fn(dy, t, y, p) -> ()
//     // to_oop_f returns out of place function handle (in my own handle format) Fn(t, y, p) -> dy
//     // to_sundials_f returns a function that can be used by the sundials crate
struct OdeFn<F>
where
    F: Fn(Realtype, &[Realtype], &mut [Realtype], &[Realtype]) -> RhsResult + '_,
{
    func: F,
}

trait ExtractVars {
    fn extract_vars(&self) -> BTreeSet<Var>;
}

impl<I: Num + Clone, R: Float + Clone> ExtractVars for Expr<I, R> {
    fn extract_vars(&self) -> BTreeSet<Var> {
        match self {
            Expr::Variable(var) => {
                let mut set = BTreeSet::new();
                set.insert(var.clone());
                set
            }
            Expr::Int(_) | Expr::Real(_) | Expr::Pi | Expr::E => BTreeSet::new(),
            Expr::Add(left, right)
            | Expr::Mul(left, right)
            | Expr::Div(left, right)
            | Expr::Pow(left, right) => {
                let mut set = left.extract_vars();
                set.extend(right.extract_vars());
                set
            }
            Expr::Der(expr, var, _) => {
                let mut set = expr.extract_vars();
                set.insert(*var.clone());
                set
            }
        }
    }
}

impl<I: Num + Clone, R: Float + Clone> ExtractVars for Equation<I, R> {
    fn extract_vars(&self) -> BTreeSet<Var> {
        let mut vars = self.lhs.extract_vars();
        vars.extend(self.rhs.extract_vars());
        vars
    }
}
// prob = ODEProblem(sys, u0, tspan, p, jac = true)
// struct ODEProblem {
//     sys: ODESystem,
//     u0: Vec<f64>,
//     tspan: (Real, Real)
//     p: Vec<f64>,
// }

fn main() {
    let t = Var {
        name: "t".to_string(),
        independent_vars: BTreeSet::new(),
    };

    let mut t_set = BTreeSet::new();
    t_set.insert(t.name.clone());

    let x = Var {
        name: "x".to_string(),
        independent_vars: t_set.clone(),
    };
    let y = Var {
        name: "y".to_string(),
        independent_vars: t_set.clone(),
    };
    let z = Var {
        name: "z".to_string(),
        independent_vars: t_set.clone(),
    };

    let sig = Var {
        name: "sig".to_string(),
        independent_vars: BTreeSet::new(),
    };
    let rho = Var {
        name: "rho".to_string(),
        independent_vars: BTreeSet::new(),
    };
    let beta = Var {
        name: "beta".to_string(),
        independent_vars: BTreeSet::new(),
    };

    // dx/dt = sig(y - x)
    let lorenz1 = Equation {
        lhs: Expr::Der(
            Box::new(Expr::Variable::<i32, f64>(x.clone())),
            Box::new(t.clone()),
            1,
        ),
        rhs: Expr::Mul(
            Box::new(Expr::Variable(sig.clone())),
            Box::new(Expr::Add(
                Box::new(Expr::Variable(y.clone())),
                Box::new(Expr::Mul(
                    Box::new(Expr::Int(-1)),
                    Box::new(Expr::Variable(x.clone())),
                )),
            )),
        ),
    };

    // dy/dt = x(rho - z) - y
    let lorenz2 = Equation {
        lhs: Expr::Der(
            Box::new(Expr::Variable::<i32, f64>(y.clone())),
            Box::new(t.clone()),
            1,
        ),
        rhs: Expr::Add(
            Box::new(Expr::Mul(
                Box::new(Expr::Variable(x.clone())),
                Box::new(Expr::Add(
                    Box::new(Expr::Variable(rho.clone())),
                    Box::new(Expr::Mul(
                        Box::new(Expr::Int(-1)),
                        Box::new(Expr::Variable(z.clone())),
                    )),
                )),
            )),
            Box::new(Expr::Mul(
                Box::new(Expr::Int(-1)),
                Box::new(Expr::Variable(y.clone())),
            )),
        ),
    };

    // dz/dt = xy - beta*z
    let lorenz3 = Equation {
        lhs: Expr::Der(
            Box::new(Expr::Variable::<i32, f64>(z.clone())),
            Box::new(t.clone()),
            1,
        ),
        rhs: Expr::Add(
            Box::new(Expr::Mul(
                Box::new(Expr::Variable(x.clone())),
                Box::new(Expr::Variable(y.clone())),
            )),
            Box::new(Expr::Mul(
                Box::new(Expr::Mul(
                    Box::new(Expr::Int(-1)),
                    Box::new(Expr::Variable(beta.clone())),
                )),
                Box::new(Expr::Variable(z.clone())),
            )),
        ),
    };
    let ivs = vec![t];
    let dvs = vec![x, y, z, sig, rho, beta];
    let independent_vars: BTreeSet<Var> = ivs.into_iter().collect();
    let dependent_vars: BTreeSet<Var> = dvs.into_iter().collect();

    let sys = ODESystem {
        equations: vec![lorenz1, lorenz2, lorenz3],
        independent_vars,
        dependent_vars,
    };
    // sys.dependent_vars.
    println!("dvs: {:?}", sys.dependent_vars);

    // let y0 = [0., 1.];
    //define the right-hand-side
    fn f(_t: Realtype, y: &[Realtype; 2], ydot: &mut [Realtype; 2], p: &[Realtype]) -> RhsResult {
        // *ydot = [0.into(), -y[1] * k];
        RhsResult::Ok
    }
    //initialize the solver
    // let mut solver = SolverNoSensi::new(
    //     LinearMultistepMethod::Bdf,
    //     f,
    //     0.,
    //     &y0,
    //     1e-4,
    //     AbsTolerance::scalar(1e-4),
    //     1e-2,
    // )
    // .unwrap();
    // //and solve
    // let ts: Vec<_> = (1..1000).collect();
    // println!("0,{},{}", y0[0], y0[1]);
    // for &t in &ts {
    //     let (_tret, &[x, xdot]) = solver.step(t as _, StepKind::Normal).unwrap();
    //     println!("{},{},{}", t, x, xdot);
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expr_variable_extraction() {
        let t = Var {
            name: "t".to_string(),
            independent_vars: BTreeSet::new(),
        };

        let expr = Expr::Variable::<i64, f64>(t.clone());
        let extracted_vars = expr.extract_vars();

        assert_eq!(extracted_vars.len(), 1);
        assert!(extracted_vars.contains(&t));
    }

    #[test]
    fn test_equation_variable_extraction() {
        let t = Var {
            name: "t".to_string(),
            independent_vars: BTreeSet::new(),
        };
        let x = Var {
            name: "x".to_string(),
            independent_vars: BTreeSet::new(),
        };

        let eq = Equation {
            lhs: Expr::Variable::<i64, f64>(t.clone()),
            rhs: Expr::Variable(x.clone()),
        };

        let extracted_vars = eq.extract_vars();

        assert_eq!(extracted_vars.len(), 2);
        assert!(extracted_vars.contains(&t));
        assert!(extracted_vars.contains(&x));
    }

    #[test]
    fn test_odesystem_variable_set() {
        let t = Var {
            name: "t".to_string(),
            independent_vars: BTreeSet::new(),
        };
        let x = Var {
            name: "x".to_string(),
            independent_vars: BTreeSet::new(),
        };

        let eq = Equation {
            lhs: Expr::Variable::<i64, f64>(t.clone()),
            rhs: Expr::Variable(x.clone()),
        };

        let system = ODESystem {
            equations: vec![eq],
            independent_vars: vec![t.clone()].into_iter().collect(),
            dependent_vars: vec![x.clone()].into_iter().collect(),
        };

        assert_eq!(system.independent_vars.len(), 1);
        assert_eq!(system.dependent_vars.len(), 1);
        assert!(system.independent_vars.contains(&t));
        assert!(system.dependent_vars.contains(&x));
    }
}
