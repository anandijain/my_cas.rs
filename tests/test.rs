// #[cfg(test)]
// mod tests {
//     use cvode_wrap::AbsTolerance;
//     use cvode_wrap::LinearMultistepMethod;
//     use cvode_wrap::RhsResult;
//     use cvode_wrap::SolverNoSensi;
//     use cvode_wrap::StepKind;
//     use my_cas::*;
//     use super::*;
//     use csv::ReaderBuilder;
//     use csv::WriterBuilder;
//     use libloading::{Library, Symbol};
//     use std::fs;
//     use std::process::Command;

//     #[test]
//     fn test_simple_ode_system() {
//         let p_val = Ex::Par("p".to_string());
//         let x_var = Ex::Var("x".to_string());

//         let dx_dt_simple = ODE {
//             variable: "x".to_string(),
//             expression: Ex::Mul(
//                 Box::new(Ex::Const(-1.)),
//                 Box::new(Ex::Mul(Box::new(p_val.clone()), Box::new(x_var.clone()))),
//             ),
//         };

//         let p_param = 2.0;
//         let defaults_simple = vec![("x".to_string(), 1.0), ("p".to_string(), p_param)];

//         let simple_sys = ODESystem {
//             odes: vec![dx_dt_simple],
//             defaults: defaults_simple,
//             tspan: (0.0, 10.0),
//         };

//         let (y0_vec, p_vec) = extract_initial_conditions_and_parameters(&simple_sys);
//         let y0: [f64; 1] = [y0_vec[0]];
//         let p: [f64; 1] = [p_vec[0]];

//         let code = generate_function_from_system(&simple_sys);
//         fs::write("generated_code_simple.rs", code).expect("Unable to write file");

//         Command::new("rustc")
//             .arg("--crate-type=cdylib")
//             .arg("generated_code_simple.rs")
//             .status()
//             .expect("Failed to compile");

//         let mut results = Vec::new();

//         unsafe {
//             let lib =
//                 Library::new("libgenerated_code_simple.dylib").expect("Failed to load the library");
//             let func: Symbol<extern "C" fn(f64, &[f64; 1], &mut [f64; 1], &[f64; 1]) -> i32> =
//                 lib.get(b"my_ode_function").expect("Function not found");

//             let mut solver = SolverNoSensi::new(
//                 LinearMultistepMethod::Bdf,
//                 |t, y, ydot, k| {
//                     func(t, y.into(), ydot.into(), k);
//                     RhsResult::Ok
//                 },
//                 0.,
//                 &y0,
//                 1e-8,
//                 AbsTolerance::scalar(1e-8),
//                 p,
//             )
//             .unwrap();

//             let ts = ndarray::Array::linspace(1.0, 10.0, 100);
//             for &t in &ts {
//                 let step = solver.step(t as _, StepKind::Normal).unwrap();
//                 println!("{:?}", step);
//                 results.push((step.0, step.1[0]));
//             }
//         }

//         // Saving results to a CSV
//         let mut wtr = WriterBuilder::new().from_path("rust_simple.csv").unwrap();
//         for &(time, value) in &results {
//             wtr.write_record(&[format!("{}", time), format!("{}", value)])
//                 .unwrap();
//         }
//         wtr.flush().unwrap();

//         // Reading in the simple.csv and comparing
//         let mut rdr = ReaderBuilder::new().from_path("simple.csv").unwrap();
//         let mut errors = Vec::new();

//         for (result, record) in results.iter().zip(rdr.records().skip(1)) {
//             let record = record.expect("a record");
//             let reference_value: f64 = record[1].parse().expect("a float value");
//             let error = (result.1 - reference_value).abs();
//             errors.push(error);
//         }

//         // Printing out the errors
//         for error in errors {
//             println!("Error: {}", error);
//         }
//     }

//     #[test]
//     fn test_time_columns_match() {
//         use std::error::Error;

//         // Read the time column from the Julia CSV
//         let mut rdr = csv::Reader::from_path("simple.csv").unwrap();
//         let mut julia_times: Vec<f64> = Vec::new();
//         for result in rdr.records() {
//             let record = result.unwrap();
//             julia_times.push(record[0].parse::<f64>().unwrap());
//         }

//         // Read the time column from the Rust CSV
//         let mut rdr = csv::Reader::from_path("rust_simple.csv").unwrap();
//         let mut rust_times: Vec<f64> = Vec::new();
//         for result in rdr.records() {
//             let record = result.unwrap();
//             rust_times.push(record[0].parse::<f64>().unwrap());
//         }

//         // Ensure the vectors are the same length
//         assert_eq!(
//             julia_times.len(),
//             rust_times.len(),
//             "The time columns have different lengths!"
//         );

//         // Check if every entry in the time column is the same
//         for (i, (jt, rt)) in julia_times.iter().zip(rust_times.iter()).enumerate() {
//             assert!(
//                 (jt - rt).abs() < 1e-6,
//                 "Mismatch at index {}: Julia time is {}, Rust time is {}",
//                 i,
//                 jt,
//                 rt
//             );
//         }
//     }

//     #[test]
//     fn ode_order_lower_test() {
//         let sys = pendulum_sys();
//         println!("{:#?}", sys);
//         println!("{}", sys.equations.len());
//         let transformed_sys = transform_to_first_order(&sys);
//         println!("{:#?}", transformed_sys);
//         println!("{}", transformed_sys.equations.len());
//     }
// }
