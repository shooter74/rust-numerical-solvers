extern crate colored;
extern crate nalgebra as na;

mod univariate_solvers;
mod univariate_minimizers;
mod nelder_mead;

use colored::Colorize;

fn rosenbrock(x: &na::DVector<f64>) -> f64 {
    return (1.0-x[0]).powi(2) + 100.0*(x[1] - x[0].powi(2)).powi(2);
}

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

fn check_result(x : f64, x_true : f64, tol : f64, test_name: &str, verbose: bool) -> u32 {
    let diff : f64 = (x - x_true).abs();
    let test_name_padded: String = format!("{:<30}", test_name);
    if diff < tol {
        if verbose {
            println!("{}\t: x = {}\tf(x) = {}\t{}", test_name_padded, format!("{:<20}", x), format!("{:<40}", fct(x)), "passed".green());
        } else {
            println!("{} {}", test_name_padded, "passed".green());
        }
        return 1;
    } else {
        if verbose {
            println!("{}\t: x = {}\tf(x) = {}\t{} (expected {}, delta = {})", test_name_padded, format!("{:<20}", x), format!("{:<40}", fct(x)), "failed".red(), x_true, (x - x_true));
        } else {
            println!("{} {} : expected {}, got {}", test_name_padded, "failed".red(), x_true, x);
        }
        return 0;
    }
}

fn check_result_optim(x: &na::DVector<f64>, f_x: f64, x_true: &na::DVector<f64>, f_x_true: f64, tol_x: f64, tol_f_x: f64, test_name: &str, verbose: bool) -> u32 {
    let diff_x : f64 = (x - x_true).norm();
    let diff_f_x : f64 = (f_x - f_x_true).abs();
    let test_name_padded: String = format!("{:<30}", test_name);
    if diff_x < tol_x && diff_f_x < tol_f_x {
        if verbose {
            println!("{}\t: x = {}\tf(x) = {}\t{}", test_name_padded, format!("{:<20}", x), format!("{:<40}", f_x), "passed".green());
        } else {
            println!("{} {}", test_name_padded, "passed".green());
        }
        return 1;
    } else {
        if verbose {
            println!("{}\t: x = {}\tf(x) = {}\t{} (expected {}, delta = {})", test_name_padded, format!("{:<20}", x), format!("{:<40}", f_x), "failed".red(), x_true, (x - x_true));
        } else {
            println!("{} {} : expected {}, got {}", test_name_padded, "failed".red(), x_true, x);
        }
        return 0;
    }
}

fn print_test_results(num_tests_passed: u32, num_tests_total: u32) {
    let ratio_str:String = format!("{}/{} ({} %)", num_tests_passed, num_tests_total, ((num_tests_passed as f64)/(num_tests_total as f64)*100.0).round());
    if num_tests_passed == num_tests_total {
        println!("{} {}", "Tests passed:", ratio_str.green());
    } else {
        println!("{} {}", "Tests passed:", ratio_str.red());
    }
}

fn test_univariate_solvers(verbose: bool) {
    println!("Testing univariate numerical solvers.");
    let x0 :       f64 = 1.0;
    let tol :      f64 = 1e-10;
    let max_iter : u32 = 100;
    let dx_num :   f64 = 1e-6;
    let mut num_tests_passed : u32 = 0;
    let num_tests_total :      u32 = 7;

    let x_mathematica: f64   = -3.26650043678562449167148755288;// 30 digits of precision
    let x_mathematica_2: f64 = -6.27133405258685307845641527902;// 30 digits of precision
    // let x_mathematica_3: f64 = -9.42553801930504142668603949182;// 30 digits of precision
    let x_newton:     f64 = univariate_solvers::newton_solve(&(fct as fn(f64) -> f64), &(dfct as fn(f64) -> f64), x0, tol, max_iter);
    let x_newton_num: f64 = univariate_solvers::newton_solve_num(&(fct as fn(f64) -> f64), x0, tol, dx_num, max_iter);
    let x_halley:     f64 = univariate_solvers::halley_solve(&(fct as fn(f64) -> f64), &(dfct as fn(f64) -> f64), &(ddfct as fn(f64) -> f64), x0, tol, max_iter, false).unwrap();
    let x_halley_num: f64 = univariate_solvers::halley_solve_num(&(fct as fn(f64) -> f64), x0, tol, dx_num, max_iter, false).unwrap();
    let x_bisection : f64 = univariate_solvers::bisection_solve(&(fct as fn(f64) -> f64), -5.0, 1.0, tol).unwrap();
    let x_secant :    f64 = univariate_solvers::secant_solve(&(fct as fn(f64) -> f64), -1.0, 1.0, tol, max_iter);
    let x_ridder :    f64 = univariate_solvers::ridder_solve(&(fct as fn(f64) -> f64), -5.0, 1.0, tol, max_iter).unwrap();
    num_tests_passed += check_result(x_newton, x_mathematica, tol, "Newton's method", verbose);
    num_tests_passed += check_result(x_newton_num, x_mathematica, tol, "Newton's method (num)", verbose);
    num_tests_passed += check_result(x_halley, x_mathematica_2, tol, "Halley's method", verbose);
    num_tests_passed += check_result(x_halley_num, x_mathematica_2, tol, "Halley's method (num)", verbose);
    num_tests_passed += check_result(x_bisection, x_mathematica, tol, "Bisection method", verbose);
    num_tests_passed += check_result(x_secant, x_mathematica, tol, "Secant method", verbose);
    num_tests_passed += check_result(x_ridder, x_mathematica, tol, "Ridder's method", verbose);

    print_test_results(num_tests_passed, num_tests_total);
}

fn test_univariate_optimizers(verbose: bool) {
    let tol :      f64 = 1e-10;
    let max_iter : u32 = 100;
    let dx_num :   f64 = 1e-6;
    let mut num_tests_passed : u32 = 0;
    let num_tests_total :      u32 = 1;

    let x_mathematica:    f64 = -4.54295618675514754103476876324;// 30 digits of precision
    let y_mathematica:    f64 = -0.206327079359226884630654987440;// 30 digits of precision
    let x_golden_section: f64 = univariate_minimizers::golden_section_minimize(&(fct as fn(f64) -> f64), -7.0, -1.0, tol);
    num_tests_passed += check_result(x_golden_section, x_mathematica, tol*1e2, "Golden section search", verbose);
    print_test_results(num_tests_passed, num_tests_total);
}

fn test_multivariate_optimizers(verbose: bool) {
    let tol :      f64 = 1e-10;
    let max_iter : u32 = 1000;
    let dx_num :   f64 = 1e-6;
    let mut num_tests_passed : u32 = 0;
    let num_tests_total :      u32 = 1;

    let tol_x:      f64 = 1e-4;
    let tol_f_x:    f64 = 1e-6;
    
    let x_true:        na::DVector<f64> = na::DVector::from_vec(vec![1.,1.]);
    let f_x_true:      f64 = rosenbrock(&x_true);
    let sol_nelder_mead: (na::DVector<f64>, f64) = nelder_mead::nelder_mead(&(rosenbrock as fn(&na::DVector<f64>) -> f64), &na::DVector::from_vec(vec![2.0,-1.0]), 0.1, tol, max_iter, false);
    // num_tests_passed += check_result(x_nelder_mead, x_true, tol*1e2, "Nelder-Mead", verbose);
    num_tests_passed += check_result_optim(&sol_nelder_mead.0, sol_nelder_mead.1, &x_true, f_x_true, tol_x, tol_f_x, "Nelder-Mead", verbose);
    print_test_results(num_tests_passed, num_tests_total);
}

fn main() {
    println!("Testing Rust numerical solvers.");
    let verbose : bool = true;
    test_univariate_solvers(verbose);
    test_univariate_optimizers(verbose);
    test_multivariate_optimizers(verbose);
}
