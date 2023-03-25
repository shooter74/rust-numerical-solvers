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

fn ddfct(x : f64) -> f64 {
    // exp(x) - 2*cos(x)/x**2 - (x**2 - 2)*sin(x)/x**3
    if x != 0.0 {
        let tmp0 : f64 = f64::powi(x, 2);
        return (f64::exp(x) - (tmp0 - 2.0)*f64::sin(x)/f64::powi(x, 3) - 2.0*f64::cos(x)/tmp0) as f64;
    } else {
        return 2.0/3.0;
    }
}

fn main() {
    println!("Testing Rust numerical solvers.");
    let x0 : f64 = 1.0;
    let tol : f64 = 1e-10;
    let max_iter : u32 = 100;
    let dx_num : f64 = 1e-6;
    let x_mathematica = -3.26650043678562449167148755288;
    let x_newton = univariate_solvers::newton_solve(&(fct as fn(f64) -> f64), &(dfct as fn(f64) -> f64), x0, tol, max_iter);
    let x_newton_num: f64 = univariate_solvers::newton_solve_num(&(fct as fn(f64) -> f64), x0, tol, dx_num, max_iter);
    let x_halley: f64 = univariate_solvers::halley_solve(&(fct as fn(f64) -> f64), &(dfct as fn(f64) -> f64), &(ddfct as fn(f64) -> f64), x0, tol, max_iter, false).unwrap();
    let x_halley_num: f64 = univariate_solvers::halley_solve_num(&(fct as fn(f64) -> f64), x0, tol, dx_num, max_iter, false).unwrap();
    let x_bisection : f64 = univariate_solvers::bisection_solve(&(fct as fn(f64) -> f64), -5.0, 1.0, tol).unwrap();
    let x_secant : f64 = univariate_solvers::secant_solve(&(fct as fn(f64) -> f64), -1.0, 1.0, tol, max_iter);
    let x_ridder : f64 = univariate_solvers::ridder_solve(&(fct as fn(f64) -> f64), -5.0, 1.0, tol, max_iter).unwrap();
    println!("Mathematica           : x = {}\tf(x) = {}", x_mathematica, fct(x_mathematica));
    println!("Newton's method       : x = {}\tf(x) = {}", x_newton, fct(x_newton));
    println!("Newton's method (num) : x = {}\tf(x) = {}", x_newton_num, fct(x_newton_num));
    println!("Halley's method       : x = {}\tf(x) = {}", x_halley, fct(x_halley));
    println!("Halley's method (num) : x = {}\tf(x) = {}", x_halley_num, fct(x_halley_num));
    println!("Bisection             : x = {}\tf(x) = {}", x_bisection, fct(x_bisection));
    println!("Secant method         : x = {}\tf(x) = {}", x_secant, fct(x_secant));
    println!("Ridder's method       : x = {}\tf(x) = {}", x_ridder, fct(x_ridder));
}
