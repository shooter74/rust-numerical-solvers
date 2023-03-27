/// Golden section search for minimizing a function f(x)
/// @param f function to minimize
/// @param a left bracket
/// @param b right bracket
/// @param tol tolerance
/// @return solution
/// @note The interval [a, b] must bracket the minimum, and the function must have f''(x) > 0 over the interval [a, b] to garantee convergence.
pub fn golden_section_minimize<F>(f : F, mut a: f64, mut b: f64, tol: f64) -> f64
where F : Fn(f64) -> f64
{
    let invphi: f64 = (f64::sqrt(5.0) - 1.0) / 2.0;   // 1 / phi
    let invphi2: f64 = (3.0 - f64::sqrt(5.0)) / 2.0;  // 1 / phi^2

    a = f64::min(a, b);
    b = f64::max(a, b);
    let mut h: f64 = b - a;
    if h <= tol {
        return (a + b)/2.0;
    }

    // Required steps to achieve tolerance
    let n: u32 = (f64::ceil(f64::ln(tol / h) / f64::ln(invphi))) as u32;

    let mut c: f64 = a + invphi2 * h;
    let mut d: f64 = a + invphi * h;
    let mut yc: f64 = f(c);
    let mut yd: f64 = f(d);

    for _ in 0..n {
        if yc < yd {  // yc > yd to find the maximum
            b = d;
            d = c;
            yd = yc;
            h = invphi * h;
            c = a + invphi2 * h;
            yc = f(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi * h;
            d = a + invphi * h;
            yd = f(d);
        }
    }

    if yc < yd { return (a + d)/2.0; }
    else { return (c + b)/2.0; }
}
