use std::collections::HashMap;

type DerivativeFn = Box<dyn Fn(&HashMap<String, f64>) -> f64>;

struct System {
    derivatives: HashMap<String, DerivativeFn>,
}

impl System {
    pub fn from_dsl(dsl: &str) -> Self {
        let mut derivatives = HashMap::new();
        let ls = dsl.lines().collect::<Vec<_>>();
        println!("ls: {:?}", ls);
        for line in ls {
            let parts: Vec<_> = line.split("~").collect();
            let var = parts[0].trim().split('(').nth(2).unwrap().split(')').next().unwrap();
            // println!("par: {}", var);
            println!("var: {}", var);
            match var {
                "x" => derivatives.insert(
                    "x".to_string(),
                    Box::new(|vars: &HashMap<String, f64>| -vars.get("x").unwrap()) as DerivativeFn,
                ),
                "y" => derivatives.insert(
                    "y".to_string(),
                    Box::new(|vars: &HashMap<String, f64>| *vars.get("y").unwrap()) as DerivativeFn,
                ),
                _ => panic!("Unsupported variable!"),
            };
        }

        System { derivatives }
    }

    pub fn simulate(&self, init: HashMap<String, f64>, dt: f64, steps: usize) {
        let mut values = init;

        for _ in 0..steps {
            let mut diffs = HashMap::new();
            for (var, fun) in &self.derivatives {
                diffs.insert(var.clone(), fun(&values));
            }

            for (var, diff) in diffs.iter() {
                let v = values.get_mut(var).unwrap();
                *v += diff * dt;
            }
        }

        // Print final values
        for (var, val) in values.iter() {
            println!("{}: {}", var, val);
        }
    }
}

fn main() {
    let system = System::from_dsl(
        "Differential(t)(x) ~ -x
         Differential(t)(y) ~ y",
    );

    let mut initial_values = HashMap::new();
    initial_values.insert("x".to_string(), 1.0);
    initial_values.insert("y".to_string(), 1.0);

    system.simulate(initial_values, 0.01, 1000);
}
