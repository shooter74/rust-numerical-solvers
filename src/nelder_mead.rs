extern crate nalgebra as na;

// let mut tuple_list2: Vec<(u16, u16)> = vec![(1, 5), (0, 17), (8, 2)];
// tuple_list2.sort_by(|a, b| a.1.cmp(&b.1));

/// Computes the mean of the points.
/// @param points The points.
/// @return The mean of the points.
fn mean_of_points(points: &Vec<(na::DVector<f64>, f64)>) -> f64 {
    let mut sum: f64 = 0.0;
    for point in points {
        sum += point.1;
    }
    return sum / points.len() as f64;
}

/// Computes the standard deviation of the points.
/// @param points The points.
/// @return The standard deviation of the points.
fn standard_deviation_of_points(points: &Vec<(na::DVector<f64>, f64)>) -> f64 {
    let mean: f64 = mean_of_points(points);
    let mut sum: f64 = 0.0;
    for value in points {
        sum += (value.1 - mean).powi(2);
    }
    return f64::sqrt(sum / points.len() as f64);
}

/// Computes the average distance between each point in the simplex.
/// @param simplex The simplex.
/// @return The average distance between each point in the simplex.
fn compute_simplex_size(simplex: &Vec<(na::DVector<f64>, f64)>) -> f64 {
    let mut sum: f64 = 0.0;
    for i in 0..simplex.len() {
        for j in 0..simplex.len() {
            if &i != &j {
                sum += (&simplex[i].0 - &simplex[j].0).norm();
            }
        }
    }
    return sum / (simplex.len() as f64 * simplex.len() as f64);
}

/// Implements the Nelder-Mead gradient-less optimization algorithm.
/// @param f  The function to optimize.
/// @param x0 The starting point of the algorithm.
/// @param simplex_size The maximum size of the simplex at initialization.
pub fn nelder_mead<F>(f: F, x0: &na::DVector<f64>, simplex_size: f64, tol: f64, max_iter: u32, verbose: bool) -> (na::DVector<f64>, f64)
where F : Fn(&na::DVector<f64>) -> f64
{
    // Parameters
    let alpha: f64 = 1.0; // Reflection coefficient
    let gamma: f64 = 2.0; // Expansion coefficient
    let rho:   f64 = 0.5; // Contraction coefficient
    let sigma: f64 = 0.5; // Shrink coefficient

    // Create simplex around x0
    let mut simplex: Vec<(na::DVector<f64>, f64)> = Vec::new();
    simplex.push((x0.clone(), f(&x0)));// Initial point of the simplex = starting point of the algorithm
    for i in 0..x0.len() {
        simplex.push((x0.clone(), 0.0));
        simplex[i+1].0[i] += simplex_size;  // Initialise each point in the simplex to be simplex_size away from x0 along aech direction
        simplex[i+1].1 = f(&simplex[i+1].0);// Evaluate objective function at each point in the simplex
    }

    let N: &usize = &simplex.len(); // Number of points in the simplex

    if verbose {
        println!("Initial simplex:");
        for v in &simplex {
            println!("{:?}", v);
        }
    }

    for iter in 0..max_iter {
        // Sort the simplex by function value
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        if verbose {
           println!("\nIteration {}\nSorted simplex:", iter);
            for v in &simplex {
                println!("{} -> {}", v.0, v.1);
            }
        }

        // Termination condition 1 : convergence of the simplex values
        let std_fvals: f64 = standard_deviation_of_points(&simplex);
        if std_fvals < tol {
            if verbose {
                println!("Converged on function values after {} iterations", iter);
            }
            return simplex[0].clone();
        }

        // Termination condition 2 : convergence of the simplex size
        let average_simplex_size: f64 = compute_simplex_size(&simplex);
        if average_simplex_size < tol {
            if verbose {
                println!("Converged on simplex size after {} iterations", iter);
            }
            return simplex[0].clone();
        }

        // Compute the centroid of the simplex (excluding the worst point)
        let mut centroid: na::DVector<f64> = na::DVector::zeros(x0.len());
        for i in 0..(N-1) {
            centroid += &simplex[i].0;
        }
        centroid /= (N-1) as f64;

        if verbose {
            println!("Centroid: {}", centroid);
        }
        
        // Compute the reflection point : xr = x0 + alpha*(x0 - x[N+1])
        let reflection: na::DVector<f64> = &centroid + alpha*(&centroid - &simplex[N-1].0);
        let f_reflection: f64 = f(&reflection);

        if verbose {
            println!("Reflection: {} -> {}", reflection, f_reflection);
        }

        // If the reflection point is better than the best point, replace the worst point with the reflection point
        if &simplex[0].1 <= &f_reflection && &f_reflection < &simplex[N-2].1 {
            simplex[N-1] = (reflection, f_reflection);
            continue;// Go to next iteration
        }

        // If the reflection point is better than the best current point, try an expansion
        if &f_reflection < &simplex[0].1 {
            let expansion: na::DVector<f64> = &centroid + gamma*(&reflection - &centroid);
            let f_expansion: f64 = f(&expansion);

            if verbose {
                println!("Expansion: {} -> {}", expansion, f_expansion);
            }

            if &f_expansion <= &f_reflection {
                simplex[N-1] = (expansion, f_expansion);
            } else {
                simplex[N-1] = (reflection, f_reflection);
            }
            continue;// Go to next iteration
        }

        // If the reflection point is worse than the second worst point, try a contraction
        if &f_reflection >= &simplex[N-2].1 {
            let contraction: na::DVector<f64> = &centroid + rho*(&simplex[N-1].0 - &centroid);
            let f_contraction: f64 = f(&contraction);

            if verbose {
                println!("Contraction: {} -> {}", contraction, f_contraction);
            }

            if &f_contraction < &simplex[N-1].1 {
                simplex[N-1] = (contraction, f_contraction);
                continue;// Go to next iteration
            } else {
                // Shrink the simplex
                for i in 1..*N {
                    simplex[i].0 = &simplex[0].0 + sigma*(&simplex[i].0 - &simplex[0].0);
                    simplex[i].1 = f(&simplex[i].0);
                }

                if verbose {
                    println!("Shrinking the whole simplex");
                }
            }
        }
    }

    if verbose {
        println!("Maximum number of iterations reached");
    }

    return simplex[0].clone();
}