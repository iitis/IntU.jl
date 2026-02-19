# Stiefel Manifold Integration

IntU supports integration over the **Stiefel manifold** $V_k(\mathbb{C}^d)$,
denoted as $V(k, d)$ or $St(k, d)$. This manifold consists of all $d \times k$
matrices $V$ with orthonormal columns, i.e., satisfying the condition:

$$V^\dagger V = I_k$$

where $I_k$ is the $k \times k$ identity matrix.

## Mathematical Definition

The uniform measure on the Stiefel manifold is induced by the Haar measure on
the unitary group $U(d)$. Specifically, a random matrix $V \in
V_k(\mathbb{C}^d)$ can be constructed by taking the first $k$ columns of a
Haar-random unitary matrix $U \in U(d)$.

This implies that for any polynomial function $f(V)$, the integral over the
Stiefel manifold is equivalent to:

$$\int_{V_k(\mathbb{C}^d)} f(V) dV = \int_{U(d)} f(U_{1..d, 1..k}) dU$$

When $k=1$, the Stiefel manifold corresponds to the set of unit vectors in
$\mathbb{C}^d$, and the measure reduces to the Fubini-Study measure on pure
quantum states $|\psi\rangle$.

## Usage

To perform integration over the Stiefel manifold, use the `dStiefel(V, d, k)` measure.

```julia
using IntU, Symbolics

@variables d
V = SymbolicMatrix(:V)
# Integration over V_2(C^d)
integrate(abs(V[1, 1])^2, dStiefel(d, 2))
```

The system automatically handles the mapping to the unitary group and applies
Weingarten calculus.

### Asymptotic Expansions

Large-$d$ expansions are fully supported:

```@example stiefel
using IntU, Symbolics
@variables d
V = SymbolicMatrix(:V)
measure = dStiefel(d, 2)
expr = abs(V[1, 1])^2 * abs(V[1, 2])^2
asymptotic(expr, measure, 2)
```

## API Reference

```@docs
IntU.StiefelMeasure
dStiefel
```
