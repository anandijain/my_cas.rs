#[no_mangle]
pub extern "C" fn my_ode_function(t: f64, y: &[f64; 1], ydot: &mut [f64; 1], p: &[f64; 1]) -> () {
    let d = [(-1.0 * (p[0] * y[0]))];
    for x in 0..y.len() {
        ydot[x] = d[x];
    }
    }