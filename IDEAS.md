# Future Feature Ideas for IntU.jl

Based on the current capabilities of `IntU.jl` (symbolic integration over the
Haar measure of $U(d)$ using Weingarten calculus), here are several high-impact
functionalities that could be added:

~~## 1. Orthogonal ($O(d)$) and Symplectic ($Sp(d)$) Groups~~
~~*   **Description**: Extend the package to support integration over random~~
~~    orthogonal and symplectic matrices.~~
~~*   **Why**: These groups are fundamental in quantum information and random~~
~~    matrix theory.~~
~~*   **Implementation**: Requires implementing the Weingarten functions for these~~
~~    groups, which involve sums over pair partitions (matchings) rather than~~
~~    permutations.~~

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

~~## 4. Symbolic Trace Logic (Coordinate-Free Integration)~~
~~*   **Description**: Instead of expanding `tr(U*A*U'*B)` into indices~~
~~    immediately, use symbolic rewrite rules to perform integration directly at~~
~~    the trace level where possible (graphical calculus).~~
~~*   **Why**: For complex expressions involving many matrix multiplications,~~
~~    expanding to indices $O(d^k)$ is computationally expensive. Graphical~~
~~    Weingarten calculus can be much faster.~~

~~## 5. Gaussian Unitary Ensemble (GUE)~~
~~*   **Description**: Support integration over Gaussian random matrices~~
~~    (GUE/GOE/GSE) where entries are i.i.d. Gaussian.~~
~~*   **Features**: Includes GUE ($\beta=2$), GOE ($\beta=1$), and GSE ($\beta=4$).~~
~~*   **Why**: Useful for comparing Haar-random results with chaotic/Hamiltonian~~
~~    dynamics. This is simpler than Haar but nice to have in the same unified~~
~~    interface.~~

~~## 6. Quantum Information Helpers~~
~~*   **Description**: Built-in functions for calculating average Purity,~~
~~    Fidelity, or Entanglement Entropy of subsystems.~~
~~*   **Example**: `average_purity(rho_symbolic, measure)`.~~

~~## 7. Pre-computed Integral Cache~~
~~*   **Description**: Create a cache with all potential interesting integrals that~~
~~    will be distributed with the package.~~
~~*   **Why**: Instant retrieval of interesting standard results without need for~~
~~    re-calculation; serves as a reference library for researchers.~~

~~## 8. Circular Ensembles (COE, CSE)~~
~~*   **Description**: Implement integration over the Circular Orthogonal (COE) and Circular Symplectic Operations (CSE).~~
~~*   **Why**: **[Parity with Haarpy]**. These ensembles (unitary matrices with specific symmetries) are distinct from the standard Haar measure on $U(d)$ and are supported by Haarpy.~~
~~*   **Implementation**: Similar to GOE/GSE but for unitary matrices; requires specific Weingarten-like expansions or mappings to the unitary group.~~

~~## 9. Permutation Groups~~
~~*   **Description**: Add explicit support for integration/summation over Permutation Groups (Symmetric Group $S_n$) and Centered Permutation Groups.~~
~~*   **Why**: **[Parity with Haarpy]**. Haarpy provides functionalities for these groups. While `IntU.jl` uses them internally for Weingarten, exposing them as a domain would allow for broader combinatorial applications.~~

~~## 10. ITensors.jl Integration~~
~~*   **Description**: Create a bridge to `ITensors.jl` to allow symbolic integration of tensor network contractions.~~
~~*   **Why**: **[Extension]**. Allows integrating large tensor networks without manually converting them to trace expressions.~~
~~*   **Implementation**: Convert `ITensor` objects into `IntU.jl`'s symbolic representation (or vice versa) and apply the integration engine.~~

## 11. SU(d) Integration
*   **Description**: Support integration over the Special Unitary group $SU(d)$.
*   **Why**: **[Extension]**. While $U(d)$ and $SU(d)$ averages often coincide for "balanced" polynomials, they differ for polynomials with non-zero winding numbers (e.g., $\det(U)$ terms).
*   **Implementation**: Incorporate the $\det(U)=1$ constraint into the Weingarten calculus.

## 12. Advanced Asymptotics & Approximations
*   **Description**: Extend the asymptotic engine to provide higher-order corrections or different limiting regimes.
*   **Why**: **[Extension]**. To study finite-size corrections in detail beyond the leading order.
