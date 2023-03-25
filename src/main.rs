mod univariate_solvers;

///  x   sin(x)
/// e  + ──────
///        x
fn fct(x : f64) -> f64 {
    if x != 0.0 {
        return f64::sin(x)/x + f64::exp(x);
    } else {
        return 2.0;
    }
}

fn dfct(x : f64) -> f64 {
    if x != 0.0 {
        return f64::exp(x) + f64::cos(x)/x - f64::sin(x)/x.powi(2);
    } else {
        return 1.0;
    }
}

fn main() {
    println!("Testing Rust numerical solvers.");
    let x0 : f64 = 1.0;
    let tol : f64 = 1e-10;
    let max_iter : u32 = 100;
    let x_mathematica = -3.26650043678562449167148755288;
    let x_newton = univariate_solvers::newton_solve(&(fct as fn(f64) -> f64), &(dfct as fn(f64) -> f64), x0, tol, max_iter);
    println!("Mathematica     : x = {}", x_mathematica);
    println!("Newton's method : x = {}", x_newton);
}
