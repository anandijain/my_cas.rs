#[no_mangle]
pub extern "C" fn my_ode_function(t: f64, y: &[f64; 3], ydot: &mut [f64; 3], p: &[f64; 3]) -> () {
    let d = [
        (p[0] * (y[1] + (-1.0 * y[0]))),
        ((y[0] * (p[1] + (-1.0 * y[2]))) + (-1.0 * y[1])),
        ((y[0] * y[1]) + ((p[2] * y[2]) * -1.0)),
    ];
    for x in 0..y.len() {
        ydot[x] = d[x];
    }
}
