extern crate nalgebra as na;

/// Gauss-Newton algorithm to solve a non-linear least squares problem. It minimizes the difference between fct_lsqr(xp, beta) and the data (xp, yp)
/// @param xp: vector of x values of the data points
/// @param yp: vector of y values of the data points
/// @param fct_lsqr: function that computes the least squares function. It takes as input the parameters and the data points and returns the model for the data fit : .
pub fn gauss_newton_lsqr<F: Fn(&na::DVector<f64>, &na::DVector<f64>) -> na::DVector<f64>>(xp: &na::DVector<f64>, yp: &na::DVector<f64>, fct_lsqr: &F, beta0: &na::DVector<f64>, tol: f64, n_iter_max: u32, dx_num: f64, verbose: bool) -> na::DVector<f64> {
    let mut beta: na::DVector<f64> = beta0.clone();

    let n_pts:  usize = xp.len();
    let n_dims: usize = beta.len();

    for iter in 0..n_iter_max {
        let residuals = yp - fct_lsqr(xp, &beta);// Residual vector

        // Compute the Jacobian
        let mut jac: na::DMatrix<f64> = na::DMatrix::zeros(n_pts, n_dims);

        let f_beta: na::DVector<f64> = fct_lsqr(&xp, &beta);
        for j in 0..n_dims {
            let mut beta_dx: na::DVector<f64> = beta.clone();
            beta_dx[j] += dx_num;
            let jac_col = (fct_lsqr(&xp, &beta_dx) - &f_beta) / dx_num;
            for i in 0..n_pts {
                jac[(i, j)] = jac_col[i];
            }
        }

        // Compute the Gauss-Newton step
        let jac_t = jac.transpose();// J^T
        let jac_t_jac = &jac_t*jac;// J^T*J
        let jac_t_res = - &jac_t*&residuals;// J^T*residuals
        let delta_beta = jac_t_jac.qr().solve(&jac_t_res).unwrap();// (J^T*J)^{-1}*J^T*residuals
        
        if verbose {
            println!("iter = {}\tbeta = {}\tresiduals = {}\tdelta_beta = {}", iter, &beta, &residuals, &delta_beta);
        }

        beta = &beta - &delta_beta;

        if delta_beta.norm() < tol {
            break;
        }
    }

    return beta;
}
