# IntU.jl

IntU.jl is a Julia package for the **symbolic** calculation of integrals over
the Haar measure of classical compact groups ($U(d)$, $O(d)$, $Sp(d)$) and
related ensembles. It leverages **Weingarten Calculus** to compute exact results
for polynomial functions of matrix entries, supporting arbitrary symbolic
dimension $d$.

## IntU in action

To introduce the main functionality of IntU, consider the problem of averaging
$|U_{i,j}|^2$ over the unitary group, i.e., computing $\int dU |U_{i,j}|^2 =
\int dU U_{i,j} U_{i,j}^*$.

While numerical approaches (like sampling random matrices) can estimate this,
they are slow and approximate. IntU provides the **exact** analytic result
instantly, even for symbolic dimensions.

```julia
using IntU, Symbolics

# Define symbolic dimension 'd' and a unitary matrix 'U'
@variables d
@variables U[1:d, 1:d]::Complex

# Define the Haar measure
measure = dU(U, d)

# Compute the integral of |U_{i,j}|^2
# Note: IntU handles symbolic indices automatically if defined, 
# but here we use concrete 1,1 for simplicity which yields the same result by symmetry.
integrate(abs(U[1,1])^2, measure)
# Output: 1 / d
```

For more complex moments, such as $\int dU |U_{1,1}|^2 |U_{1,2}|^2$, IntU handles the combinatorics (Weingarten functions) automatically:

```julia
integrate(abs(U[1,1])^2 * abs(U[1,2])^2, measure)
# Output: 1 / (d * (1 + d))
```

## IntU functionality

IntU implements Weingarten calculus for the Unitary, Orthogonal, and Symplectic
groups, as well as Haar-random pure states.

### Unitary group
Unitary matrices $U$ are complex matrices satisfying $U U^\dagger = I_d$. One
can calculate averages over the unitary Haar measure using `dU` and `integrate`.

```julia
# 4-th moment of a diagonal entry
integrate(abs(U[1,1])^4, dU(U, d))
# Output: 2 / (d * (1 + d))
```

### Orthogonal group
Orthogonal matrices $O$ are real matrices satisfying $O O^T = I_d$. Averages are
computed using the `dO` measure.

```julia
@variables O_mat[1:d, 1:d]::Real
integrate(O_mat[1,1]^4, dO(O_mat, d))
# Output: 3 / (d * (2 + d))
```

### Symplectic group
Symplectic matrices $S$ are unitary matrices of even dimension $2n$ that
preserve the symplectic form, $S \Omega S^T = \Omega$. Use `dSp`.

```julia
@variables S_mat[1:d, 1:d]::Complex
integrate(abs(S_mat[1,1])^2, dSp(S_mat, d))
# Output: 1 / d
```

### Random Pure States
IntU can integrate polynomial functions of the components of a Haar-random pure
state vector $|\psi\rangle$ of dimension $d$.

```julia
@variables dim
@variables psi[1:dim]::Complex
measure_psi = dPsi(psi, dim)

# Average of |ψ_1|^2
integrate(abs(psi[1])^2, measure_psi)
# Output: 1 / dim
```

### Symbolic Traces
IntU supports index-free notation for integrating traces of products of random
matrices, which is often more convenient for quantum information tasks.

```julia
using IntU: tr
# Define symbolic matrices A, B (constant) and U (random)
A = SymbolicMatrix(:A)
B = SymbolicMatrix(:B)
U_sym = SymbolicMatrix(:U, false, :U) # unitary

# Compute ∫ tr(U A U† B) dU
expr = tr(U_sym * A * U_sym' * B)
integrate(expr, dU(U, d))
# Output: tr(A)*tr(B) / d
```

## Installation

IntU is tested with Julia 1.11 or later. Installation can be done through the
Pkg REPL:

```julia
import Pkg; Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## How to cite this work

If you use IntU.jl in your research, please cite:

```bibtex
@misc{intu2024,
  author = {Pawela, Łukasz and Krawiec, Adam},
  title = {IntU.jl: Symbolic integration over the Haar measure of classical compact groups},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/iitis/IntU.jl}}
}
```
