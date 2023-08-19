use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use std::collections::HashMap;
use std::fs;
use std::process::Command;

#[derive(Debug, Clone)]
enum Ex {
    Const(f64),
    Var(String), // implicity dependence on time
    Par(String), // independent of t
    Mul(Box<Ex>, Box<Ex>),
    Add(Box<Ex>, Box<Ex>),
    Pow(Box<Ex>, Box<Ex>),
}

#[derive(Debug, Clone)]
struct ODE {
    variable: String,
    expression: Ex,
}

#[derive(Debug, Clone)]
struct ODESystem {
    odes: Vec<ODE>,
    defaults: HashMap<String, f64>,
    tspan: (f64, f64),
}

fn translate_expression_to_code(
    expr: &Ex,
    var_mapping: &HashMap<String, usize>,
    par_mapping: &HashMap<String, usize>,
) -> String {
    match expr {
        Ex::Const(c) => {
            if c.fract() == 0.0 {
                format!("{:.1}", c)
            } else {
                c.to_string()
            }
        }
        Ex::Var(v) => {
            let idx = var_mapping
                .get(v)
                .expect(&format!("Variable {} not found", v));
            format!("y[{}]", idx)
        }
        Ex::Par(p) => {
            let idx = par_mapping
                .get(p)
                .expect(&format!("Parameter {} not found", p));
            format!("p[{}]", idx)
        }
        Ex::Mul(a, b) => format!(
            "({} * {})",
            translate_expression_to_code(a, var_mapping, par_mapping),
            translate_expression_to_code(b, var_mapping, par_mapping)
        ),
        Ex::Add(a, b) => format!(
            "({} + {})",
            translate_expression_to_code(a, var_mapping, par_mapping),
            translate_expression_to_code(b, var_mapping, par_mapping)
        ),
        Ex::Pow(a, b) => format!(
            "({}.powf({}))",
            translate_expression_to_code(a, var_mapping, par_mapping),
            translate_expression_to_code(b, var_mapping, par_mapping)
        ),
    }
}

fn generate_function_from_system(system: &ODESystem) -> String {
    let state_count = system.odes.len();

    // Generate mapping for variables and parameters
    let mut var_mapping = HashMap::new();
    for (idx, ode) in system.odes.iter().enumerate() {
        var_mapping.insert(ode.variable.clone(), idx);
    }
    let mut par_mapping = HashMap::new();
    let mut idx = 0;
    for (k, _) in &system.defaults {
        if !var_mapping.contains_key(k) {
            par_mapping.insert(k.clone(), idx);
            idx += 1;
        }
    }

    let param_count = par_mapping.len();

    let mut result = format!(
        "#[no_mangle]\npub extern \"C\" fn my_ode_function(t: f64, y: &[f64; {}], ydot: &mut [f64; {}], p: &[f64; {}]) -> () {{\n",
        state_count, state_count, param_count
    );

    result += "    let d = [";
    for ode in &system.odes {
        result += &format!(
            "{},",
            translate_expression_to_code(&ode.expression, &var_mapping, &par_mapping)
        );
    }
    result.truncate(result.len() - 1); // Remove the last comma
    result += "];\n";

    result += "    for x in 0..y.len() {\n        ydot[x] = d[x];\n    }\n";
    result += "    }";
    result
}

fn main() {
    // Define sigma, rho, and beta parameters
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;

    let sig = Ex::Par("sigma".to_string());
    let rh = Ex::Par("rho".to_string());
    let bet = Ex::Par("beta".to_string());

    let x = Ex::Var("x".to_string());
    let y = Ex::Var("y".to_string());
    let z = Ex::Var("z".to_string());

    // ODEs for the Lorenz system
    let dx_dt = ODE {
        variable: "x".to_string(),
        expression: Ex::Mul(
            Box::new(sig),
            Box::new(Ex::Add(
                Box::new(y.clone()),
                Box::new(Ex::Mul(Box::new(Ex::Const(-1.)), Box::new(x.clone()))),
            )),
        ),
    };

    let dy_dt = ODE {
        variable: "y".to_string(),
        expression: Ex::Add(
            Box::new(Ex::Mul(
                Box::new(x.clone()),
                Box::new(Ex::Add(
                    Box::new(rh),
                    Box::new(Ex::Mul(Box::new(Ex::Const(-1.)), Box::new(z.clone()))),
                )),
            )),
            Box::new(Ex::Mul(Box::new(Ex::Const(-1.)), Box::new(y.clone()))),
        ),
    };

    let dz_dt = ODE {
        variable: "z".to_string(),
        expression: Ex::Add(
            Box::new(Ex::Mul(Box::new(x), Box::new(y.clone()))),
            Box::new(Ex::Mul(Box::new(bet), Box::new(z.clone()))),
        ),
    };

    // Default values and parameters
    let mut defaults = HashMap::new();
    defaults.insert("x".to_string(), 1.0);
    defaults.insert("y".to_string(), 1.0);
    defaults.insert("z".to_string(), 1.0);
    defaults.insert("sigma".to_string(), sigma);
    defaults.insert("rho".to_string(), rho);
    defaults.insert("beta".to_string(), beta);

    // Constructing the Lorenz ODESystem
    let lorenz_sys = ODESystem {
        odes: vec![dx_dt, dy_dt, dz_dt],
        defaults: defaults,
        tspan: (0.0, 100.0), // As an example, can be adjusted
    };

    println!("{:?}", lorenz_sys);

    let code = generate_function_from_system(&lorenz_sys);

    fs::write("generated_code.rs", code).expect("Unable to write file");

    // Compile the generated code into a shared library
    Command::new("rustc")
        .arg("--crate-type=cdylib")
        .arg("generated_code.rs")
        .status()
        .expect("Failed to compile");

    // Load the shared library

    // Load the function from the library
    unsafe {
        let lib = Library::new("libgenerated_code.dylib").expect("Failed to load the library");
        let func: Symbol<extern "C" fn(f64, &[f64; 3], &mut [f64; 3], &[f64; 3]) -> i32> =
            lib.get(b"my_ode_function").expect("Function not found");
        // Here you'd call the function with appropriate arguments

        let y0 = [1.0, 1.0, 1.0];
        let p = [10.0, 28.0, 8.0 / 3.0];

        //initialize the solver
        let mut solver = SolverNoSensi::new(
            LinearMultistepMethod::Adams,
            |t, y, ydot, k| {
                func(t, y.into(), ydot.into(), k);
                RhsResult::Ok
            },
            0.0,
            &y0,
            1e-4,
            AbsTolerance::scalar(1e-4),
            p,
        )
        .unwrap();

        let ts: Vec<_> = (1..1000).collect();
        println!("0,{:?}", y0);
        for &t in &ts {
            let step = solver.step(t as _, StepKind::Normal).unwrap();
            // let (_tret, &[x, xdot, z]) = 
            println!("{},{:?}", step.0, step.1);
        }
    }
}
