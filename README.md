# IntU.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://iitis.github.io/IntU.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://iitis.github.io/IntU.jl/dev/)

IntU.jl is a Julia package for the exact symbolic and numeric calculation of
integrals over the Haar measure of classical compact groups ($U(d)$, $O(d)$,
$Sp(d)$) and related ensembles. It leverages **Weingarten Calculus** to compute
exact results for polynomial functions of matrix entries, with broad support
for symbolic dimension $d$ in entry-wise and trace-polynomial workflows.

For detailed documentation, please visit [iitis.github.io/IntU.jl](https://iitis.github.io/IntU.jl).

## IntU in action

To introduce the main functionality of IntU, consider the problem of averaging
$|U_{i,j}|^2$ over the unitary group, i.e., computing
$\int dU\, |U_{i,j}|^2 = \int dU\, U_{i,j}\,\overline{U_{i,j}}$.
This is the diagonal special case of the general second-moment correlator
$\int dU\, U_{i,j}\,\overline{U_{k,l}}$ (obtained by setting $k=i$, $l=j$).

IntU provides exact analytic results through a unified interface:
`integrate(expr, measure)`. It supports matrix-valued expressions (for concrete
integer result dimensions) and provides the `@integrate` macro for intuitive
symbolic integration. Symbolic dimensions are supported for a broad class of
workflows; selected paths currently require concrete integer dimensions
(notably higher pure trace moments `|tr(U)|^(2k), k > 1`, `hciz` on
`SymbolicMatrix` inputs, and direct matrix-valued integration of
`SymbolicMatrix`/`SymbolicMatrixProduct` expressions).

The `@integrate` macro implicitly identifies random matrices based on the measure:
- `dU`, `dSU`, `dCUE`, `dDesign` $\rightarrow$ `U` (Unitary / CUE / t-design)
- `dO` $\rightarrow$ `O` (Orthogonal)
- `dSp` $\rightarrow$ `Sp` (Symplectic)
- `dPerm` $\rightarrow$ `P` (Permutation)
- `dCPerm` $\rightarrow$ `Y` (Centered Permutation)
- `dCOE`, `dCSE` $\rightarrow$ `S` (Circular Orthogonal/Symplectic)
- `dGUE`, `dGOE`, `dGSE` $\rightarrow$ `H` (Gaussian Ensembles)
- `dGinUE`, `dGinOE`, `dGinSE` $\rightarrow$ `G` (Ginibre Ensembles)
- `dPsi` $\rightarrow$ `psi` (Pure State)
- `dDiagUnitary` $\rightarrow$ `D` (Diagonal Unitary)
- `dStiefel` $\rightarrow$ `V` (Stiefel Manifold)

Unknown symbols (like `A`, `B`, `d`) are automatically treated as constants or dimensions.

> [!NOTE]
> **Symbol Scope and Redefinition**: The `@integrate` macro manages a persistent symbolic state in your session. If you use a symbol (e.g., `U`) as a random matrix in one call and then use it in a different context (e.g., as a constant in an orthogonal integral), the macro automatically re-declares and re-binds the symbol to the correct type and measure context. This "Safety Rebind" mechanism prevents silent math errors when running examples in sequence.

```julia
using IntU, Symbolics, LinearAlgebra

# Integrate the matrix expression U * U'
# The @integrate macro knows 'U' is the random matrix for measure dU(2)
res = @integrate U * U' dU(2)
# Output: [1.0 0.0; 0.0 1.0] (2x2 Identity Matrix)
```

```julia
using IntU
# Compute the integral of |U_{1,1}|^2
@integrate abs(U[1, 1])^2 dU(d)
# Output: 1 / d
```

For more complex moments, such as $\int dU |U_{1,1}|^2 |U_{1,2}|^2$, IntU handles the combinatorics (Weingarten functions) automatically:

```julia
@integrate abs(U[1, 1])^2 * abs(U[1, 2])^2 dU(d)
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
@integrate abs(U[1, 1])^4 dU(d)
# Output: 2 / (d * (1 + d))
```

### Special Unitary Group
The Special Unitary group $SU(d)$ consists of unitary matrices with determinant 1. Use `dSU`.

```julia
@integrate abs(U[1, 1])^2 dSU(d)
# Output: 1/d
```

For the currently supported stable-range workflow, `dSU` delegates to the same
balanced Weingarten rules as `dU`. As a result, unbalanced monomials evaluate
to zero:

```julia
@integrate U[1, 1] dSU(d)
# Output: 0
```

### Orthogonal group
Orthogonal matrices $O$ are real matrices satisfying $O O^T = I_d$. Averages are
computed using the `dO` measure.

```julia
@integrate O[1, 1]^4 dO(d)
# Output: 3 / (d * (2 + d))
```

### Symplectic group
Symplectic matrices $Sp$ are unitary matrices of even dimension $2n$ that
preserve the symplectic form, $Sp \Omega Sp^T = \Omega$. Use `dSp`.

```julia
@integrate abs(Sp[1, 1])^2 dSp(d)
# Output: 1 / d
```

### Ginibre Ensembles
Ginibre ensembles consist of non-Hermitian matrices with i.i.d. Gaussian entries. 
- **GinUE (Complex Ginibre Ensemble)**: i.i.d. complex Gaussian entries. Use `dGinUE`.
- **GinOE (Real Ginibre Ensemble)**: i.i.d. real Gaussian entries. Use `dGinOE`.
- **GinSE (Symplectic Ginibre Ensemble)**: i.i.d. quaternionic Gaussian entries. Use `dGinSE`.

```julia
# E[Tr(G G')] = d^2
@integrate tr(G * G') dGinUE(d)
```

### Circular Ensembles
IntU also supports Circular Ensembles (CUE, COE, CSE) which are commonly used in random matrix theory.
- **CUE (Circular Unitary Ensemble)**: Equivalent to the Haar measure on $U(d)$. Use `dCUE`.
- **COE (Circular Orthogonal Ensemble)**: Ensemble of symmetric unitary matrices. Use `dCOE`.
- **CSE (Circular Symplectic Ensemble)**: Ensemble of self-dual unitary matrices of even dimension $2n$. Use `dCSE`.

```julia
# COE moment E[|S_{1,1}|^2]
@integrate abs(S[1, 1])^2 dCOE(d)
# Output: 2 / (d + 1)
```

### Random Pure States
IntU can integrate polynomial functions of the components of a Haar-random pure
state vector $|\psi\rangle$ of dimension $d$.

```julia
# Average of |ψ_1|^2
# The @integrate macro treats 'psi' as a (d, 1) column vector
@integrate abs(psi[1, 1])^2 dPsi(d)
# Output: 1 / d
```

### Permutation Group
IntU supports integration over the Symmetric group $S_d$ (permutation matrices).
- **dPerm**: Integration over the set of $d \times d$ permutation matrices.

```julia
# E[P_11]
@integrate P[1, 1] dPerm(d)
# Output: 1 / d
```

IntU also supports centered permutation matrices $Y = P - J/d$.
- **dCPerm**: Integration over centered permutation matrices.

```julia
# E[Y_11^2]
@integrate Y[1, 1]^2 dCPerm(d)
# Output: (d - 1) / d^2
```

### Diagonal Unitary Matrices (Torus group)
IntU supports integration over the group of diagonal unitary matrices, which
corresponds to independent phase averaging for each diagonal entry.

```julia
# E[|D_11|^2]
@integrate abs(D[1, 1])^2 dDiagUnitary(d)
# Output: 1
```

### Symbolic Traces
IntU supports index-free notation for integrating traces of products of random
matrices, which is often more convenient for quantum information tasks.

```julia
# Compute ∫ tr(U A U† B) dU
@integrate tr(U * A * U' * B) dU(d)
# Output: (tr(A)*tr(B)) / d
```

### Harish-Chandra-Itzykson-Zuber (HCIZ) Integrals
IntU supports calculating HCIZ integrals of the form $\int_{U(d)} dU \exp(\text{Tr}(A U B U^\dagger))$.

> [!NOTE]
> HCIZ has two user paths: eigenvalue vectors and matrix inputs. For
> `SymbolicMatrix` inputs, `hciz` requires a concrete integer dimension (symbolic
> `d` is not supported). For numeric matrix inputs with degenerate spectra, IntU
> sorts both spectra and applies tiny independent perturbations to both before
> evaluating the determinant formula.

```julia
using IntU, LinearAlgebra
A = diagm([1.0, 2.0])
B = diagm([0.5, 1.5])
# Compute ∫ dU exp(Tr(A U B U'))
hciz(A, B)
# Output: 20.9329...
```

### Stiefel Manifolds
IntU supports integration over the Stiefel manifold $V_k(\mathbb{C}^d)$, which represents the set of $d \times k$ matrices with orthonormal columns. This generalizes Haar-random pure states ($k=1$).

```julia
# E[|V_{1,1}|^2]
@integrate abs(V[1, 1])^2 dStiefel(d, 2)
# Output: 1 / d
```

## ITensors.jl Integration

IntU.jl provides a bridge to [ITensors.jl](https://github.com/ITensors/ITensors.jl) for symbolic integration of tensor networks.

```julia
using IntU, ITensors
i = Index(2, "Out")
j = Index(2, "In")
i2 = Index(2, "Out2")
j2 = Index(2, "In2")

# Mark a Haar-random unitary U and its adjoint U_dag
U = ITensorUnitary(out_indices=[i], in_indices=[j])
U_dag = ITensorUnitary(out_indices=[j2], in_indices=[i2], is_adj=true)

# Integrate a balanced network E[Tr(U A U_dag B)] over U(2)
A = randomITensor(j, j2)
B = randomITensor(i2, i)
res = integrate([U, A, U_dag, B], dU(2))
```

## Running Examples and Benchmarks

IntU.jl comes with a comprehensive set of examples and benchmarks.

### Examples

You can run all examples at once using the provided script:

```bash
sh examples/runexamples.sh
```

To run a single example, pass the example number to the script. For instance, to run `01_basics.jl`:

```bash
sh examples/runexamples.sh 01
```

### Benchmarks

To run all benchmarks:

```bash
sh benchmarks/runbenchmarks.sh
```

To run a single benchmark, pass the number to the script:

```bash
sh benchmarks/runbenchmarks.sh 01
```

## Installation

IntU is tested with Julia 1.11 or 1.12. Installation can be done through the
Pkg REPL:

```julia
import Pkg; Pkg.add(url="https://github.com/iitis/IntU.jl")
```

## How to cite this work

If you use IntU.jl in your research, please cite:

```bibtex
@misc{intu2024,
  author = {Pawela, Łukasz and Puchała, Zbigniew},
  title = {IntU.jl: Symbolic integration over the Haar measure of classical compact groups},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/iitis/IntU.jl}}
}
```
