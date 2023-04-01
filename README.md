# rust-numerical-solvers
A collection of numerical solvers in Rust. This library exists mainly as a coding exercise, but is useful nonetheless.
The solvers that are ~~barred~~ are not yet implemented but are planned to be.

All solvers and optimizers are designed to work with 64-bit floating point numbers (`f64` type in rust).

## Contents of the library

Several types of solvers and optimizers exist in this library:
- Solvers
    - Univariate
        - Gradient-based methods
        - Non gradient-based methods
    - Multivariate
        - Gradient-based methods
- Optimizers
    - Univariate
        - Gradient-based methods
        - Non gradient-based methods
    - Multivariate
        - Gradient-based methods
        - Non gradient-based methods

Some methods that require derivative/gradient information have a fully-numerical version, suffixed `_num`, where derivatives are evaluated using finite-differences. Prefer providing analytical gradients to the methods when possible.

### Solvers

Solvers allow to solve non-linear functions: $$ f(x) = 0 $$ or $$ f(\vec{x}) = \vec{0} $$

#### Univariate solvers

Univariate solvers allow to solve univariate non-linear functions : $$ f(x) = 0 $$

##### Gradient-based solvers

Gradient based solvers use the gradient of the function to be solved in order to help with the convergence of the numerical method. Some methods, such as Halley's or Laguerre's, even require the second order derivative.

Here is a list of the univariate gradient-based solvers implemented in the library :
- Newton's method (`newton_solve`)
- Newton's method with finite-differences derivatives (`newton_solve_num`)
- Halley's method (`halley_solve`)
- Halley's method with finite-differences derivatives (`halley_solve_num`)
- ~~Laguerre's method (`laguerre_solve`)~~
- ~~Laguerre's method with finite-differences derivatives (`laguerre_solve_num`)~~

##### Derivative-free solvers

Derivative-free solvers, as they name suggest, do not require the derivatives of the function in order to solve it. They are usually more robust, but suffer from a lower convergence rate. They are however very useful when the function to solve does not have a closed-form derivative, or when it is too noisy.

Here is a list of the univariate derivative-free solvers implemented in the library :
- Bisection method (`bisection_solve`)
- Secant method (`secant_solve`)
- Ridder's method (`ridder_solve`)

#### Multivariate solvers


### Optimizers

Optimizers allow to minimize an objective function and solve the following class of problems : $$ \min_{x} f(x) $$ or $$ \min_{\vec{x}} f(\vec{x}) $$

#### Univariate optimizers

Univariate optimizers allow to minimize univariate non-linear functions : $$ \min_{x} f(x) $$

##### Gradient-based optimizers

Gradient-based optimizers require the gradient of the function in order to speed-up convergence.

##### Derivative-free optimizers

Derivative-free optimizers do not require the gradient of the function and are generally more robust to noisy objective functions. Some optimizers are local and other are global. Local optimizers tend to converge to the minimum that is closest to the given starting point. Global minimizers are capable of exploring the solution space more extensively and can sometimes return the global minimum of the objective function.

Here is a list of univariate derivative-free optimizers implemented in the library :
- Golden section search (`golden_section_minimize`)
- ~~Cubic Lagrange polynomial optimization (`cubic_lagrange_minimize`)~~

#### Multivariate optimizers

Multivariate optimizers allow to minimize multivariate non-linea functions : $$ \min_{\vec{x}} f(\vec{x}) $$

##### Gradient-based optimizers

Gradient-based optimizers require the gradient of the function in order to speed-up convergence. Some of them even require the Hessian matrix (second derivative) to be computed (Quasi-Newton method for example).

##### Derivative-free optimizers

Derivative-free optimizers do not require the gradient of the function and are generally more robust to noisy objective functions.

Here is a list of multivariate derivative-free optimizers implemented in the library :
- Nelder-Mead (`nelder_mead`)
- ~~Particle Swarm Optimization (`particle_swarm_minimize`)~~
- ~~Differential evolution (`differential_evolution_minimize`)~~

#### Least-squares solvers

The following non-linear least-squares problem can be solved more efficiently using specialised techniques than generic optimizers :$$ \min_{\beta} \sum_{i=0}^{N} (f(x_i, \beta) - y_i)^2 $$

Here is a list of multivariate non-linear least-squares solvers implemented in the library :
- Gauss-Newton (`gauss_newton_lsqr`)