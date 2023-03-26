# rust-numerical-solvers
A collection of numerical solvers in Rust. This library exists mainly as a coding exercise, but is useful nonetheless.

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

Some methods that require derivative/gradient information have a fully-numerical version, where derivatives are evaluated using finite-differences. Prefer providing analytical gradients to the methods.
