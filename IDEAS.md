# Future Feature Ideas for IntU.jl

Based on the current capabilities of `IntU.jl` (symbolic integration over the
Haar measure of $U(d)$ using Weingarten calculus), here are several high-impact
functionalities that could be added:

## 1. Orthogonal ($O(d)$) and Symplectic ($Sp(d)$) Groups
*   **Description**: Extend the package to support integration over random
    orthogonal and symplectic matrices.
*   **Why**: These groups are fundamental in quantum information and random
    matrix theory.
*   **Implementation**: Requires implementing the Weingarten functions for these
    groups, which involve sums over pair partitions (matchings) rather than
    permutations.

~~## 2. Haar Random States (integration over $|\psi\rangle$)~~
~~*   **Description**: Add a dedicated measure `dPsi(psi, d)` for integration over~~
~~    Haar-random pure states.~~
~~*   **Why**: Easier syntax for problems involving state averages rather than~~
~~    full unitaries.~~
~~*   **Implementation**: This is mathematically equivalent to integrating over~~
~~    the first column of a Haar unitary, so it can reuse the existing `U(d)`~~
~~    engine with fixed indices.~~

~~## 3. Large-$d$ Asymptotic Expansions~~
~~*   **Description**: A function `asymptotic(expr, measure, order=1)` that~~
~~    returns the series expansion of the integral in powers of $1/d$.~~
~~*   **Why**: In many physical applications (like holography or many-body~~
~~    physics), the large-$d$ limit is the most relevant regime.~~

## 4. Symbolic Trace Logic (Coordinate-Free Integration)
*   **Description**: Instead of expanding `tr(U*A*U'*B)` into indices
    immediately, use symbolic rewrite rules to perform integration directly at
    the trace level where possible (graphical calculus).
*   **Why**: For complex expressions involving many matrix multiplications,
    expanding to indices $O(d^k)$ is computationally expensive. Graphical
    Weingarten calculus can be much faster.

## 5. Gaussian Unitary Ensemble (GUE)
*   **Description**: Support integration over Gaussian random matrices
    (GUE/GOE/GSE) where entries are i.i.d. Gaussian.
*   **Why**: Useful for comparing Haar-random results with chaotic/Hamiltonian
    dynamics. This is simpler than Haar but nice to have in the same unified
    interface.

~~## 6. Quantum Information Helpers~~
~~*   **Description**: Built-in functions for calculating average Purity,~~
~~    Fidelity, or Entanglement Entropy of subsystems.~~
~~*   **Example**: `average_purity(rho_symbolic, measure)`.~~

## 7. Pre-computed Integral Cache
*   **Description**: Create a cache with all potential interesting integrals that
    will be distributed with the package.
*   **Why**: Instant retrieval of interesting standard results without need for
    re-calculation; serves as a reference library for researchers.
