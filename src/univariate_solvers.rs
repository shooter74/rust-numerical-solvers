/// @brief Newton's method for solving a function f(x) = 0
/// @param f function to solve
/// @param df derivative of function f
/// @param x0 initial guess
/// @param tol tolerance
/// @param max_iter maximum number of iterations
/// @return solution
pub fn newton_solve<F, F2>(f : F, df : F2, x0 : f64, tol : f64, max_iter : u32) -> f64
    where F : Fn(f64) -> f64, F2 : Fn(f64) -> f64
{
    let mut x: f64 = x0;
    let mut dx: f64;
    let mut fx: f64;
    let mut dfx: f64;
    for _iter in 0..max_iter {
        fx = f(x);
        dfx = df(x);
        if dfx == 0.0 {
            dx = fx;
        } else {
            dx = fx/dfx;
        }
        x -= dx;
        if f64::abs(dx) < tol {
            break;
        }
    }
    return x;
}

/// @brief Newton's method for solving a function f(x) = 0
/// @param f function to solve
/// @param df derivative of function f
/// @param x0 initial guess
/// @param tol tolerance
/// @param max_iter maximum number of iterations
/// @return solution
pub fn newton_solve_num<F>(f : F, x0 : f64, tol : f64, dx_num : f64, max_iter : u32) -> f64
where F : Fn(f64) -> f64
{
    return newton_solve(&f, |x: f64| {
        (f(x + dx_num) - f(x - dx_num))/(2.0*dx_num)
    }, x0, tol, max_iter);
}

// --------------------------------------------------------------------
// ------------------------ Bracketing methods ------------------------
// --------------------------------------------------------------------

pub fn bisection_solve<F>(f : F, mut a : f64, mut b : f64, tol : f64) -> Result<f64, &'static str>
where F : Fn(f64) -> f64
{
    let mut c: f64;
    let mut fa: f64 = f(a);
    let mut fb: f64 = f(b);
    let mut fc: f64;
    let max_iter: u32 = (f64::log2((b-a)/tol)).ceil() as u32;
    for _iter in 0..max_iter {
        c = (a + b)/2.0;
        fc = f(c);
        if fa*fc < 0.0 {
            b = c;
            fb = fc;
        } else if fb*fc < 0.0 {
            a = c;
            fa = fc;
        } else {
            return Result::Err("The interval does not backet the root.");
        }
    }
    return Result::Ok((a + b)/2.0);
}