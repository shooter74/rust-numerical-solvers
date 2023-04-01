extern crate nalgebra as na;

#[path = "./xorwow.rs"]
mod xorwow;

use xorwow::Xorwow;

/// A particle in the particle swarm optimization algorithm
struct Particle {
    x: na::DVector<f64>,// position
    v: na::DVector<f64>,// velocity
    fx: f64,            // function value at x
}

impl Particle {
    fn new(x: na::DVector<f64>, v: na::DVector<f64>, fx: f64) -> Particle {
        Particle {
            x: x,
            v: v,
            fx: fx,
        }
    }
}

pub fn particle_swarm_minimize<F: Fn(&na::DVector<f64>) -> f64>(f: F, n_particles: u32, lb: &na::DVector<f64>, ub: &na::DVector<f64>, tol: f64, n_iter_max: u32, rng_seed: u32) -> f64 {
    let mut particles: Vec<Particle> = Vec::new();

    let mut rng = Xorwow::new(rng_seed);

    for _ in 0..n_particles {
        let x = lb + (ub - lb) * rng.rand_vec(lb.len());
        let v = (lb - ub) * rng.rand_vec(lb.len());
        let fx = f(&x);
        println!("x = {}\tv = {}\tf(x) = {}", &x, &v, fx);// DEBUG
        particles.push(Particle::new(x, v, fx));
    }

    let x = (lb + ub) / 2.0;
    return f(&x);
}