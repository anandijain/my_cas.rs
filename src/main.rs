use cvode_wrap::{AbsTolerance, LinearMultistepMethod, RhsResult, SolverNoSensi, StepKind};
use libloading::{Library, Symbol};
use ndarray;
use std::collections::HashMap;
use std::ffi::c_void;
use std::fs;
use std::process::Command;
use sundials_sys::*;

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
    defaults: Vec<(String, f64)>,
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
    println!("var_mapping: {:?}", var_mapping);
    println!("par_mapping: {:?}", par_mapping);

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
    result += "}";
    result
}

fn extract_initial_conditions_and_parameters(system: &ODESystem) -> (Vec<f64>, Vec<f64>) {
    // Extract y0 from the system
    let y0: Vec<f64> = system
        .odes
        .iter()
        .map(|ode| {
            system
                .defaults
                .iter()
                .find(|&&(ref var, _)| *var == ode.variable)
                .map(|&(_, value)| value)
                .unwrap() // Default to 0.0 if not found (or panic if you prefer)
        })
        .collect();

    // Extract parameters from the system
    let p: Vec<f64> = system
        .defaults
        .iter()
        .filter(|&&(ref var, _)| !system.odes.iter().any(|ode| ode.variable == *var))
        .map(|&(_, value)| value)
        .collect();

    (y0, p)
}

fn simple_ode() {
    unsafe extern "C" fn rhs(
        _t: realtype,
        y: N_Vector,
        dy: N_Vector,
        _user_data: *mut c_void,
    ) -> i32 {
        *N_VGetArrayPointer(dy) = -*N_VGetArrayPointer(y);
        0
    }

    unsafe {
        let y = N_VNew_Serial(1);
        *N_VGetArrayPointer(y) = 1.0;

        let mut cvode_mem = CVodeCreate(CV_ADAMS);

        CVodeInit(cvode_mem, Some(rhs), 0.0, y);
        CVodeSStolerances(cvode_mem, 1e-6, 1e-8);

        let matrix = SUNDenseMatrix(1, 1);
        let solver = SUNDenseLinearSolver(y, matrix);

        CVodeSetLinearSolver(cvode_mem, solver, matrix);
        CVodeSetInitStep(cvode_mem, 0.1);
        let mut t = 0f64;
        CVode(cvode_mem, 1.0, y, &mut t, CV_NORMAL);
        // y[0] is now exp(-1)

        let result = (*N_VGetArrayPointer(y) * 1e6) as i32;
        assert_eq!(result, 367879);

        N_VDestroy(y);
        CVodeFree(&mut cvode_mem);
        SUNLinSolFree(solver);
        SUNMatDestroy(matrix);
    }
}

fn main() {
    simple_ode();
    // Define sigma, rho, and beta parameters

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
            Box::new(sig.clone()),
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
                    Box::new(rh.clone()),
                    Box::new(Ex::Mul(Box::new(Ex::Const(-1.)), Box::new(z.clone()))),
                )),
            )),
            Box::new(Ex::Mul(Box::new(Ex::Const(-1.)), Box::new(y.clone()))),
        ),
    };

    let dz_dt = ODE {
        variable: "z".to_string(),
        expression: Ex::Add(
            Box::new(Ex::Mul(Box::new(x.clone()), Box::new(y.clone()))),
            Box::new(Ex::Mul(
                Box::new(Ex::Mul(Box::new(bet.clone()), Box::new(z.clone()))),
                Box::new(Ex::Const(-1.)),
            )),
        ),
    };

    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;

    // Default values and parameters
    let defaults = vec![
        ("x".to_string(), 1.0),
        ("y".to_string(), 1.0),
        ("z".to_string(), 1.0),
        ("sigma".to_string(), sigma),
        ("rho".to_string(), rho),
        ("beta".to_string(), beta),
    ];

    // Constructing the Lorenz ODESystem
    let lorenz_sys = ODESystem {
        odes: vec![dx_dt, dy_dt, dz_dt],
        defaults: defaults,
        tspan: (0.0, 100.0), // As an example, can be adjusted
    };

    let (y0_vec, p_vec) = extract_initial_conditions_and_parameters(&lorenz_sys);
    let y0: [f64; 3] = [y0_vec[0], y0_vec[1], y0_vec[2]]; // As an example for size 3
    let p: [f64; 3] = [p_vec[0], p_vec[1], p_vec[2]]; // As an example for size 3

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

        //initialize the solver
        let mut solver = SolverNoSensi::new(
            LinearMultistepMethod::Adams,
            |t, y, ydot, k| {
                func(t, y.into(), ydot.into(), k);
                RhsResult::Ok
            },
            0.0,
            &y0,
            1e-8,
            AbsTolerance::scalar(1e-8),
            p,
        )
        .unwrap();

        let ts = ndarray::Array::linspace(0.1, 10.0, 10000);
        println!("0,{:?}", y0);
        for &t in &ts {
            let step = solver.step(t as _, StepKind::Normal).unwrap();
            // let (_tret, &[x, xdot, z]) =
            // println!("{},{:?}", step.0, step.1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use csv::ReaderBuilder;
    use csv::WriterBuilder;
    use libloading::{Library, Symbol};
    use std::fs;
    use std::process::Command;

    #[test]
    fn test_simple_ode_system() {
        let p_val = Ex::Par("p".to_string());
        let x_var = Ex::Var("x".to_string());

        let dx_dt_simple = ODE {
            variable: "x".to_string(),
            expression: Ex::Mul(
                Box::new(Ex::Const(-1.)),
                Box::new(Ex::Mul(Box::new(p_val.clone()), Box::new(x_var.clone()))),
            ),
        };

        let p_param = 2.0;
        let defaults_simple = vec![("x".to_string(), 1.0), ("p".to_string(), p_param)];

        let simple_sys = ODESystem {
            odes: vec![dx_dt_simple],
            defaults: defaults_simple,
            tspan: (0.0, 10.0),
        };

        let (y0_vec, p_vec) = extract_initial_conditions_and_parameters(&simple_sys);
        let y0: [f64; 1] = [y0_vec[0]];
        let p: [f64; 1] = [p_vec[0]];

        let code = generate_function_from_system(&simple_sys);
        fs::write("generated_code_simple.rs", code).expect("Unable to write file");

        Command::new("rustc")
            .arg("--crate-type=cdylib")
            .arg("generated_code_simple.rs")
            .status()
            .expect("Failed to compile");

        let mut results = Vec::new();

        unsafe {
            let lib =
                Library::new("libgenerated_code_simple.dylib").expect("Failed to load the library");
            let func: Symbol<extern "C" fn(f64, &[f64; 1], &mut [f64; 1], &[f64; 1]) -> i32> =
                lib.get(b"my_ode_function").expect("Function not found");

            let mut solver = SolverNoSensi::new(
                LinearMultistepMethod::Bdf,
                |t, y, ydot, k| {
                    func(t, y.into(), ydot.into(), k);
                    RhsResult::Ok
                },
                0.,
                &y0,
                1e-8,
                AbsTolerance::scalar(1e-8),
                p,
            )
            .unwrap();

            let ts = ndarray::Array::linspace(1.0, 10.0, 100);
            for &t in &ts {
                let step = solver.step(t as _, StepKind::Normal).unwrap();
                println!("{:?}", step);
                results.push((step.0, step.1[0]));
            }
        }

        // Saving results to a CSV
        let mut wtr = WriterBuilder::new().from_path("rust_simple.csv").unwrap();
        for &(time, value) in &results {
            wtr.write_record(&[format!("{}", time), format!("{}", value)])
                .unwrap();
        }
        wtr.flush().unwrap();

        // Reading in the simple.csv and comparing
        let mut rdr = ReaderBuilder::new().from_path("simple.csv").unwrap();
        let mut errors = Vec::new();

        for (result, record) in results.iter().zip(rdr.records().skip(1)) {
            let record = record.expect("a record");
            let reference_value: f64 = record[1].parse().expect("a float value");
            let error = (result.1 - reference_value).abs();
            errors.push(error);
        }

        // Printing out the errors
        for error in errors {
            println!("Error: {}", error);
        }
    }

    #[test]
    fn test_time_columns_match() {
        use std::error::Error;

        // Read the time column from the Julia CSV
        let mut rdr = csv::Reader::from_path("simple.csv").unwrap();
        let mut julia_times: Vec<f64> = Vec::new();
        for result in rdr.records() {
            let record = result.unwrap();
            julia_times.push(record[0].parse::<f64>().unwrap());
        }

        // Read the time column from the Rust CSV
        let mut rdr = csv::Reader::from_path("rust_simple.csv").unwrap();
        let mut rust_times: Vec<f64> = Vec::new();
        for result in rdr.records() {
            let record = result.unwrap();
            rust_times.push(record[0].parse::<f64>().unwrap());
        }

        // Ensure the vectors are the same length
        assert_eq!(
            julia_times.len(),
            rust_times.len(),
            "The time columns have different lengths!"
        );

        // Check if every entry in the time column is the same
        for (i, (jt, rt)) in julia_times.iter().zip(rust_times.iter()).enumerate() {
            assert!(
                (jt - rt).abs() < 1e-6,
                "Mismatch at index {}: Julia time is {}, Rust time is {}",
                i,
                jt,
                rt
            );
        }
    }
}
