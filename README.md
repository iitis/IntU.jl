# IntU.jl

IntU.jl is a Julia package for the **symbolic** calculation of integrals over
the Haar measure of classical compact groups ($U(d)$, $O(d)$, $Sp(d)$) and
related ensembles. It leverages **Weingarten Calculus** to compute exact results
for polynomial functions of matrix entries, supporting arbitrary symbolic
dimension $d$.

For detailed documentation, please visit [iitis.github.io/IntU.jl](https://iitis.github.io/IntU.jl).

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

# New: Convenient measure constructors
measure = dU(d) # No matrix variable required
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

### Special Unitary Group
The Special Unitary group $SU(d)$ consists of unitary matrices with determinant 1. Use `dSU`.

```julia
integrate(abs(U[1,1])^2, dSU(U, d))
# Output: 1/d
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

### Circular Ensembles
IntU also supports Circular Ensembles (CUE, COE, CSE) which are commonly used in random matrix theory.
- **CUE (Circular Unitary Ensemble)**: Equivalent to the Haar measure on $U(d)$. Use `dCUE`.
- **COE (Circular Orthogonal Ensemble)**: Ensemble of symmetric unitary matrices. Use `dCOE`.
- **CSE (Circular Symplectic Ensemble)**: Ensemble of self-dual unitary matrices of even dimension $2n$. Use `dCSE`.

```julia
@variables S_coe[1:d, 1:d]::Complex
# COE moment E[|S_{1,1}|^2]
integrate(abs(S_coe[1,1])^2, dCOE(S_coe, d))
# Output: 2 / (d + 1)
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

### Permutation Groups
IntU supports integration over the Symmetric Group $S_d$ (permutation matrices) and centered permutation matrices $Y = P - J/d$.

```julia
@variables P[1:d, 1:d]
measure = dPerm(P, d)
# E[P_11 * P_22]
integrate(P[1,1] * P[2,2], measure)
# Output: 1 / (d * (d - 1))

@variables Y[1:d, 1:d]
m_centered = dCPerm(Y, d)
# E[Y_11^2]
integrate(Y[1,1]^2, m_centered)
# Output: (d - 1) / d^2
```

### Diagonal Unitary Matrices (Torus group)
IntU supports integration over the group of diagonal unitary matrices, which
corresponds to independent phase averaging for each diagonal entry.

```julia
@variables V[1:d, 1:d]::Complex
measure = dDiagUnitary(V, d)

# E[|V_11|^2]
integrate(abs(V[1,1])^2, measure)
# Output: 1
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
integrate(expr, dU(d))
# Output: tr(A)*tr(B) / d
```

## ITensors.jl Integration

IntU.jl provides a bridge to [ITensors.jl](https://github.com/ITensors/ITensors.jl) for symbolic integration of tensor networks.

```julia
using IntU, ITensors
i, j = Index(2), Index(2)
U_it = randomITensor(i, j)
U = ITensorUnitary(U_it; out_indices=[i], in_indices=[j])

A = randomITensor(j, i)
res = integrate([U, A], dU(2))
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
