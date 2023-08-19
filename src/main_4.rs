use std::collections::HashMap;

enum SymbolicExpression {
    Variable(String),
    Parameter(String),
    Constant(f64),
    Multiplication(Box<SymbolicExpression>, Box<SymbolicExpression>),
    Addition(Box<SymbolicExpression>, Box<SymbolicExpression>),
    // ... add other operations as needed
}

struct ODE {
    variable: String,
    expression: SymbolicExpression,
}

fn evaluate_expression(expr: &SymbolicExpression, values: &HashMap<String, f64>) -> f64 {
    match expr {
        SymbolicExpression::Variable(v) => *values.get(v).expect("Variable value not provided!"),
        SymbolicExpression::Parameter(p) => *values.get(p).expect("Parameter value not provided!"),
        SymbolicExpression::Constant(c) => *c,
        SymbolicExpression::Multiplication(a, b) => evaluate_expression(&**a, values) * evaluate_expression(&**b, values),
        SymbolicExpression::Addition(a, b) => evaluate_expression(&**a, values) + evaluate_expression(&**b, values),
        // ... handle other operations similarly
    }
}

fn generate_ode_fn(ode: ODE, values: HashMap<String, f64>) -> Box<dyn Fn(f64, f64) -> f64> {
    // let mut updated_values = values.clone();
    Box::new(move |_t, y| {
        // updated_values.insert(ode.variable.clone(), y);
        evaluate_expression(&ode.expression, &values)
    })
}

// Simple Euler's method
fn euler_solve(ode_fn: &dyn Fn(f64, f64) -> f64, t0: f64, tf: f64, y0: f64, dt: f64) -> Vec<(f64, f64)> {
    let mut t = t0;
    let mut y = y0;
    let mut results = vec![(t0, y0)];

    while t < tf {
        y = y + dt * ode_fn(t, y);
        t += dt;
        results.push((t, y));
    }

    results
}

fn main() {
    // Defining the ODE dx/dt = -param*x
    let ex =SymbolicExpression::Multiplication(
            Box::new(SymbolicExpression::Variable("x".to_string())),
            Box::new(SymbolicExpression::Parameter("param".to_string())),
        );
        let ex2 = SymbolicExpression::Addition(Box::new(ex), Box::new(SymbolicExpression::Constant(1.0)));
    let my_ode = ODE {
        variable: "x".to_string(),
        expression: ex2,
    };

    let mut values = HashMap::new();
    values.insert("param".to_string(), -2.0);
    values.insert("x".to_string(), 1.0); // Initial condition

    let ode_function = generate_ode_fn(my_ode, values);
    let x2 = ode_function(0.0, 1.0);
    println!("x2: {}", x2);
    let results = euler_solve(&*ode_function, 0.0, 2.0, 1.0, 0.1);

    for (t, y) in results {
        println!("t: {}, x: {}", t, y);
    }
}
