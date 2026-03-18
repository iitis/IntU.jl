# Diagonal Unitary Matrices (Torus Group)

Integration over the group of diagonal unitary matrices, also known as the
**torus group** $T^d$.

## Mathematical Background

A diagonal unitary matrix $D \in T^d$ has the form

$$D = \mathrm{diag}(e^{i\theta_1}, e^{i\theta_2}, \dots, e^{i\theta_d}),$$

where each phase $\theta_k \in [0, 2\pi)$ is drawn independently and uniformly.
Integration over $T^d$ is a product of independent phase averages:

$$\int_{T^d} f(D)\, dD = \prod_{k=1}^d \left( \frac{1}{2\pi} \int_0^{2\pi} f(e^{i\theta_k})\, d\theta_k \right).$$

### Integration rule for monomials

A monomial in the matrix entries $D_{ij}$ integrates to a non-zero value only
if both conditions hold simultaneously:

1. **All entries are diagonal**: $i = j$ for every factor.
2. **Phase balance**: for each index $k$, the number of $D_{kk}$ factors equals
   the number of $\bar{D}_{kk}$ factors.

$$\int_{T^d} D_{i_1 i_1} \cdots D_{i_n i_n}\, \bar{D}_{j_1 j_1} \cdots \bar{D}_{j_n j_n}\, dD = \begin{cases} 1 & \{i_1,\dots,i_n\} = \{j_1,\dots,j_n\} \text{ as multisets,} \\ 0 & \text{otherwise.} \end{cases}$$

### Physical context

Diagonal unitaries arise naturally in:
- **Dephasing channels**: random independent phase errors applied to each
  computational basis state.
- **Quantum coherence theory**: measuring robustness of off-diagonal density
  matrix elements to local phase noise.
- **Randomised benchmarking**: certain protocol families use torus-group
  averages as a computationally cheap proxy for full Haar randomness.

## Usage

Use the `dDiagUnitary(d)` measure. The `@integrate` macro automatically
identifies `D` as the random diagonal unitary.

### Basic examples

```julia
using IntU, Symbolics
@variables d

# |D_{11}|^2 = 1 deterministically (unit modulus)
@integrate abs(D[1, 1])^2 dDiagUnitary(d)
# Output: 1

# Off-diagonal entries are structurally zero
@integrate abs(D[1, 2])^2 dDiagUnitary(d)
# Output: 0
```

### Independence of distinct phases

Different diagonal entries carry independent phases, so their cross-moments
factorise:

```julia
# E[D_{11} * conj(D_{22})] = 0  (independent mean-zero phases)
@integrate D[1, 1] * conj(D[2, 2]) dDiagUnitary(d)
# Output: 0

# E[|D_{11}|^2 * |D_{22}|^2] = 1  (product of unit-modulus entries)
@integrate abs(D[1, 1])^2 * abs(D[2, 2])^2 dDiagUnitary(d)
# Output: 1
```

### Phase-balance condition

Higher moments vanish unless the multisets of $D$ and $\bar{D}$ indices match:

```julia
# Balanced: two D_{11} and two conj(D_{11})
@integrate D[1, 1]^2 * conj(D[1, 1])^2 dDiagUnitary(d)
# Output: 1

# Unbalanced: two D_{11}, only one conj(D_{11})
@integrate D[1, 1]^2 * conj(D[1, 1]) dDiagUnitary(d)
# Output: 0
```

### Manual integration

```julia
using IntU, Symbolics
@variables d
D = SymbolicMatrix(:D, :DiagUnitary)

integrate(abs(D[1, 1])^2, dDiagUnitary(d))         # 1
integrate(D[1, 1] * conj(D[2, 2]), dDiagUnitary(d)) # 0
integrate(abs(D[1, 2])^2, dDiagUnitary(d))           # 0
```

## Performance Note

Integration over $T^d$ is significantly faster than over the full unitary group
$U(d)$. The engine bypasses Weingarten functions entirely and reduces to a
multiset equality check on diagonal indices — an $\mathcal{O}(n \log n)$
operation in the number of factors, with no symbolic summation.

## See Also

- [Unitary Integration](unitary_integration.md) — full $U(d)$ Haar measure
- [Asymptotic Expansions](asymptotic.md) — large-$d$ behaviour

See [`dDiagUnitary`](@ref) in the [API Reference](api.md).
