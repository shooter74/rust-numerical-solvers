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
/// @param x0 initial guess
/// @param tol tolerance
/// @param dx_num numerical differentiation step size
/// @param max_iter maximum number of iterations
/// @return solution
/// @note This method uses numerical differentiation to compute the first and second derivatives.
pub fn newton_solve_num<F>(f : F, x0 : f64, tol : f64, dx_num : f64, max_iter : u32) -> f64
where F : Fn(f64) -> f64
{
    return newton_solve(&f, |x: f64| {
        (f(x + dx_num) - f(x - dx_num))/(2.0*dx_num)
    }, x0, tol, max_iter);
}

/// Halley's method for solving a function f(x) = 0
/// @param f function to solve
/// @param df derivative of function f
/// @param ddf second derivative of function f
/// @param x0 initial guess
/// @param tol tolerance
/// @param max_iter maximum number of iterations
/// @return solution
/// @note This method is more efficient than Newton's method, but requires the second derivative of f
pub fn halley_solve<F, F2, F3>(f: F, df: F2, ddf: F3, x0: f64, tol: f64, max_iter: u32, verbose: bool) -> Result<f64, &'static str>
where F : Fn(f64) -> f64, F2 : Fn(f64) -> f64, F3 : Fn(f64) -> f64
{
    let mut x: f64 = x0;
    let mut f_x: f64;
    let mut df_x: f64;
    let mut ddf_x: f64;
    for _i in 0..max_iter {
        f_x = f(x);
        df_x = df(x);
        ddf_x = ddf(x);
        if verbose {
            println!("x = {}, f(x) = {}, df(x) = {}, ddf(x) = {}", x, f_x, df_x, ddf_x);
        }
        if f64::abs(f_x) < tol {
            return Ok(x);
        }
        x = x - 2.0*f_x*df_x / (2.0*df_x.powi(2) - f_x*ddf_x);
    }
    return Err("Halley method did not converge after reaching the maximum number of iterations allowed.")
}

/// Halley's method for solving a function f(x) = 0
/// @param f function to solve
/// @param x0 initial guess
/// @param tol tolerance
/// @param dx_num numerical differentiation step size
/// @param max_iter maximum number of iterations
/// @return solution
/// @note This method is more efficient than Newton's method, but requires the second derivative of f.
/// @note This method uses numerical differentiation to compute the first and second derivatives.
pub fn halley_solve_num<F>(f: F, x0: f64, tol: f64, dx_num : f64, max_iter: u32, verbose: bool) -> Result<f64, &'static str>
where F : Fn(f64) -> f64
{
    let mut x: f64 = x0;
    let mut f_x: f64;
    let mut df_x: f64;
    let mut ddf_x: f64;
    let mut fx_m_dx: f64;
    let mut fx_p_dx: f64;
    for _i in 0..max_iter {
        f_x = f(x);
        fx_m_dx = f(x - dx_num);
        fx_p_dx = f(x + dx_num);
        df_x = (fx_p_dx - fx_m_dx)/(2.0*dx_num);
        ddf_x = (fx_p_dx - 2.0*f_x + fx_m_dx) / (dx_num.powi(2));
        if verbose {
            println!("x = {}, f(x) = {}, df(x) = {}, ddf(x) = {}", x, f_x, df_x, ddf_x);
        }
        if f64::abs(f_x) < tol {
            return Ok(x);
        }
        x = x - 2.0*f_x*df_x / (2.0*df_x.powi(2) - f_x*ddf_x);
    }
    return Err("Halley method did not converge after reaching the maximum number of iterations allowed.")
}

// --------------------------------------------------------------------
// ------------------------ Bracketing methods ------------------------
// --------------------------------------------------------------------

/// @brief Bisection method for solving a function f(x) = 0
/// @param f function to solve
/// @param a left bracket
/// @param b right bracket
/// @param tol tolerance
/// @return solution
/// @note The interval [a, b] must bracket the root, meaning f(a) and f(b) must be of a different sign.
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

/// @brief Secant method for solving a function f(x) = 0
/// @param f function to solve
/// @param a left bracket
/// @param b right bracket
/// @param tol tolerance
/// @param max_iter maximum number of iterations
/// @return solution
/// @note The interval [a, b] does not have to bracket the root.
/// @note The secant method is not guaranteed to converge.
pub fn secant_solve<F>(f : F, mut a : f64, mut b : f64, tol : f64, max_iter : u32) -> f64
where F : Fn(f64) -> f64
{
    let mut c: f64 = (a + b)/2.0;
    let mut fa: f64 = f(a);
    let mut fb: f64 = f(b);
    let mut fc: f64;
    for _iter in 0..max_iter {
        // c is x[n], a is x[n-1], b is x[n-2]
        c = a - fa*(a - b)/(fa - fb);
        fc = f(c);
        b = a;
        fb = fa;
        a = c;
        fa = fc;
        if (b - a).abs() < tol {
            break;
        }
    }
    return c;
}

/// @brief Ridder's method for solving a function f(x) = 0
/// @param f function to solve
/// @param a left bracket
/// @param b right bracket
/// @param tol tolerance
/// @return solution
/// @note The interval [a, b] must bracket the root, meaning f(a) and f(b) must be of a different sign.
pub fn ridder_solve<F>(f : F, mut a : f64, mut b : f64, tol : f64, max_iter : u32) -> Result<f64, &'static str>
where F : Fn(f64) -> f64
{
    let mut fa: f64 = f(a);
    let mut fb: f64 = f(b);
    let mut c: f64;
    let mut fc: f64;
    let mut s: f64;
    let mut dx: f64;
    let mut x: f64;
    let mut fx: f64;
    let mut x_old: f64 = (a + b)/2.0;
    if fa == 0.0 { return Ok(a); }
    if fb == 0.0 { return Ok(b); }
    if fa*fb > 0.0 {
        return Err("Root is not bracketed")
    }
    for i in 0..max_iter {
        // Compute the improved root x from Ridder's formula
        c = 0.5*(a + b); fc = f(c);
        s = f64::sqrt(fc.powi(2) - fa*fb);
        if s != 0.0 {
            dx = (c - a)*fc/s;
        } else {
            dx = (c - a)*fc;
        }
        if (fa - fb) < 0.0 { dx = -dx; }
        x = c + dx; fx = f(x);
        // Test for convergence
        if i > 0 {
            if f64::abs(x - x_old) < tol*f64::max(f64::abs(x),1.0) { return Ok(x) }
        }
        x_old = x;
        // Re-bracket the root as tightly as possible
        if fc*fx > 0.0 {
            if fa*fx < 0.0 { b = x; fb = fx; }
            else { a = x; fa = fx; }
        } else {
            a = c; b = x; fa = fc; fb = fx;
        }
    }
    return Err("Maximum number of iterations exceeded.")
}