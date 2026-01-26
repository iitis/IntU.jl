# Pure State Integration

Integration over the measure of random pure states $|\psi\rangle$ in
$\mathbb{C}^d$ (Fubini-Study measure).

## Theory

A random pure state $|\psi\rangle$ distributed according to the Haar measure can
be simulated by taking the first column of a Haar-random unitary matrix $U$.

```math
\int f(\psi) d\psi = \int_{U(d)} f(U |0\rangle) dU
```

Where $|0\rangle = (1, 0, \dots, 0)^T$. This implies $\psi_i = U_{i,1}$.

IntU.jl implements `dPsi` by internally mapping the state variable $\psi$ to the
first column of an underlying symbolic Unitary matrix and invoking the standard
$U(d)$ integration engine.

## Usage

```julia
using IntU, Symbolics

@variables d
@variables psi[1:d]::Complex
measure = dPsi(psi, d)

# Calculate average overlap with a fixed state phi
@variables phi[1:d]::Complex
expr = abs(sum(conj(psi[i]) * phi[i] for i in 1:d))^2
integrate(expr, measure)
# Result: sum(|phi_i|^2) / d = 1/d if phi is normalized
```

## Pitfalls

- **Normalization**: The measure is normalized such that $\int d\psi = 1$.
- **Vector indexing**: Ensure you index `psi` as a vector `psi[i]` rather than a
  matrix `psi[i,1]`.

## References

- Życzkowski, K., & Sommers, H. J. (2001). Induced measures in the space of
  mixed quantum states. *Journal of Physics A: Mathematical and General*,
  34(35), 7111.
