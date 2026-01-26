# Theory

## Symbolic Integration over the Haar Measure

This package performs symbolic integration of polynomial functions of unitary
matrix elements over the Haar measure of the unitary group $U(d)$. The core
method relies on the **Weingarten calculus** [CoSn06], which reduces integrals
of matrix coefficients to combinatorial sums involving characters of the
symmetric group.

### The Problem

We are interested in computing integrals of the form:

```math
I = \int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU
```
where $dU$ is the normalized Haar measure on $U(d)$. If the number of $U$ and
$\bar{U}$ factors differs, the integral vanishes due to phase invariance.

### Schur-Weyl Duality

The theoretical foundation for the Weingarten formula lies in **Schur-Weyl Duality**. Consider the tensor product space $(\mathbb{C}^d)^{\otimes n}$. There are two natural group actions on this space:
1. The unitary group $U(d)$ acting diagonally: $U \cdot (v_1 \otimes \dots \otimes v_n) = (Uv_1) \otimes \dots \otimes (Uv_n)$.
2. The symmetric group $S_n$ permuting the tensor factors: $\sigma \cdot (v_1 \otimes \dots \otimes v_n) = v_{\sigma^{-1}(1)} \otimes \dots \otimes v_{\sigma^{-1}(n)}$.

Schur-Weyl duality states that the actions of $U(d)$ and $S_n$ commute and span each other's commutants. Consequently, any operator $A$ acting on $(\mathbb{C}^d)^{\otimes n}$ that commutes with all unitary matrices must be a linear combination of permutation operators $P_\sigma$.

The integral we wish to compute can be viewed as the projection onto the $U(d)$-invariant subspace of $(\mathbb{C}^d \otimes (\mathbb{C}^d)^*)^{\otimes n} \cong \text{End}((\mathbb{C}^d)^{\otimes n})$. The result is an element of the commutant of $U(d)$, and thus lies in the algebra spanned by permutations.

### The Weingarten Formula

The general formula for the integral of a monomial in entries of a
Haar-distributed unitary matrix $U \in U(d)$ is given by [CoSn06]:

```math
\int_{U(d)} U_{i_1 j_1} \dots U_{i_n j_n} \bar{U}_{k_1 l_1} \dots \bar{U}_{k_n l_n} dU = \sum_{\sigma, \tau \in S_n} \delta_{i, k_\sigma} \delta_{j, l_\tau} \text{Wg}(\sigma \tau^{-1}, d)
```

where:
- $S_n$ is the symmetric group of degree $n$.
- $\delta_{i, k_\sigma} = \prod_{m=1}^n \delta_{i_m, k_{\sigma(m)}}$.
- $\text{Wg}(\pi, d)$ is the **Weingarten function** on $S_n$.

### The Weingarten Function

The Weingarten function $\text{Wg}(\pi, d)$ depends only on the cycle type of
the permutation $\pi$ and the dimension $d$. It can be computed via the theory
of Schur functions:

```math
\text{Wg}(\pi, d) = \frac{1}{(n!)^2} \sum_{\lambda \vdash n, l(\lambda) \le d} \frac{\chi^\lambda(\pi) \chi^\lambda(e)^2}{s_{\lambda}(1^d)}
```

where:
- $\lambda \vdash n$ sums over partitions of $n$.
- $\chi^\lambda$ is the irreducible character of $S_n$ indexed by $\lambda$.
- $s_{\lambda}(1^d)$ is the Schur polynomial evaluated at identity (dimension of
  the irrep of $U(d)$).

For $d < n$, the sum is restricted to partitions with length $\le d$. However,
`IntU.jl` treats $d$ symbolically and assumes $d \ge n$ for the generic rational
function form, which is valid for symbolic manipulation.

### Explicit Values

For small values of $n$, the Weingarten functions have simple rational forms.

**$n=1$**:
```math
\text{Wg}([1], d) = \frac{1}{d}
```

**$n=2$**:
```math
\text{Wg}([1,1], d) = \frac{1}{d^2-1}, \quad \text{Wg}([2], d) = \frac{-1}{d(d^2-1)}
```

**$n=3$**:
```math
\text{Wg}([1,1,1], d) = \frac{d^2-2}{d(d^2-1)(d^2-4)}, \quad \text{Wg}([2,1], d) = \frac{-1}{(d^2-1)(d^2-4)}, \quad \text{Wg}([3], d) = \frac{2}{d(d^2-1)(d^2-4)}
```

Note that singularities occur at integers $d < n$, reflecting that the integral
over $U(d)$ is different when the "length" of the partition exceeds the
dimension. `IntU.jl` assumes $d \ge n$ for the generic form.

### Large-$d$ Asymptotics

For large dimension $d$, the Weingarten function behaves asymptotically as:

```math
\text{Wg}(\pi, d) = d^{-n - |\pi|} \prod_{c \in \text{cycles}(\pi)} (-1)^{|c|-1} C_{|c|-1} + O(d^{-n-|\pi|-2})
```
where $|\pi|$ is the minimal number of transpositions to generate $\pi$ ($n$
minus number of cycles), and $C_k$ are the Catalan numbers.

$1/d$ for large $d$.

### Integration over Pure States

Integration over Haar-random pure states $|\psi\rangle$ is based on the
observation that the first column of a Haar-random unitary matrix $U$ is a
Haar-random pure state. Specifically, if $U \sim \text{Haar}(U(d))$, then:
```math
|\psi\rangle = U |0\rangle = \sum_{i=1}^d U_{i,1} |i\rangle
```
is a random state vector distributed uniformly on the complex unit sphere in
$\mathbb{C}^d$.

Thus, integrals over $|\psi\rangle$ can be mapped to integrals over $U$:
```math
\int f(\psi) d\psi = \int_{U(d)} f(U e_1) dU
```
In `IntU.jl`, this is handled by substituting $\psi_i \to U_{i,1}$ and then
invoking the standard Weingarten engine.

### Quantum Information Applications

The package provides helpers for common quantities in Quantum Information
Theory.

#### Average Purity
The purity of a quantum state $\rho$ is defined as $\text{Tr}(\rho^2)$. If
$\rho$ transforms under a random unitary channel $\rho \to U \rho_{in}
U^\dagger$, we often want to compute:
```math
\mathbb{E}_U [\text{Tr}((U \rho_{in} U^\dagger)^2)]
```
This is a polynomial of degree 2 in $U$ and $\bar{U}$, which can be evaluated
exactly using the Weingarten formula for $n=2$.

#### Average Fidelity
The fidelity between a fixed state $\sigma$ and a randomly rotated state $\rho =
U \rho_{in} U^\dagger$ is:
```math
F(\rho, \sigma) = \text{Tr}(\rho \sigma) = \text{Tr}(U \rho_{in} U^\dagger \sigma)
```
Its average is:
```math
\mathbb{E}_U [F(\rho, \sigma)] = \frac{\text{Tr}(\rho_{in}) \text{Tr}(\sigma)}{d}
```
assuming trace-1 normalization, this simplifies to $1/d$.

## Symbolic Trace Logic

While the standard Weingarten formula (Eq. 37) works with explicit indices
$U_{ij}$, it is often more convenient to work with coordinate-free expressions
involving traces of matrix products, such as $\operatorname{Tr}(U A U^\dagger
B)$. This approach is known as the **Graphical Weingarten Calculus**.

### Theoretical Basis

The integral $\int dU U_{i_1 j_1} \dots \bar{U}_{k_n l_n}$ can be interpreted as
a sum over pairings of the input and output indices of $U$ and $\bar{U}$. In the
graphical notation:
- Each $U$ is a box with an input wire (index $j$) and an output wire (index
  $i$).
- Each $\bar{U}$ (or $U^\dagger$) is a box with an input wire (index $l$) and an
  output wire (index $k$).
- The Weingarten formula sums over permutations $\sigma, \tau \in S_n$ that
  connect these wires:
  - $\sigma$ connects the output of the $m$-th $U$ to the output of the
    $\sigma(m)$-th $U^\dagger$ (indices $i$ and $k$).
  - $\tau$ connects the input of the $m$-th $U$ to the input of the $\tau(m)$-th
    $U^\dagger$ (indices $j$ and $l$).

When we integrate a trace like $\operatorname{Tr}(U A U^\dagger B)$, we are
essentially closing the loops of indices with the constant matrices $A$ and $B$.

### Algorithm

The symbolic trace integration algorithm in `IntU.jl` operates as follows:
1.  **Identify $U$ and $U^\dagger$ sites**: The expression is treated as a
    cyclic word of symbolic matrices.
2.  **Generate Wires**: Constant matrices between unitary factors form "wires"
    or partial traces.
3.  **Sum over Permutations**: For each pair $(\sigma, \tau)$:
    - Connect $U$ and $U^\dagger$ indices according to $\sigma$ and $\tau$.
    - Traverse the resulting graph. Cycles formed by traversing the graph
      correspond to full traces of products of the constant matrices encountered
      along the path.
    - Each cycle contributes a factor $\operatorname{Tr}(\dots)$.
    - If a cycle contains no constant matrices, it contributes a factor of $d$.
4.  **Weighting**: The contribution is weighted by
    $\operatorname{Wg}(\sigma\tau^{-1}, d)$.

This allows computing results like:
```math
\int_{U(d)} \operatorname{Tr}(U A U^\dagger B) dU = \frac{\operatorname{Tr}(A)\operatorname{Tr}(B)}{d}
```
symbolically, without expanding into $d^4$ tensor indices.

### References

- **[CoSn06]** Collins, B., & Śniady, P. (2006). Integration with respect to the
  Haar measure on unitary, orthogonal and symplectic groups. *Communications in
  Mathematical Physics*, 264(3), 773-795.
  [arXiv:math-ph/0402073](https://arxiv.org/abs/math-ph/0402073)
- **[Broud16]** Brouder, C. (2016). A Mathematica package for symbolic
  integration over the unitary group. *Journal of Symbolic Computation*.
- **[DiaSha98]** Diaconis, P., & Shahshahani, M. (1998). On the eigenvalues of
  random matrices. *Journal of Applied Probability*, 31A, 49-62.
- **[PuMi17]** Puchała, Z., & Miszczak, J. A. (2017). Symbolic integration with
  respect to the Haar measure on the unitary groups. *Bulletin of the Polish
  Academy of Sciences Technical Sciences*, 65(1), 21-27.
  [arXiv:1109.4244](https://arxiv.org/abs/1109.4244)
- **[Novak14]** Novak, J. (2014). Three lectures on free probability. *Random
  Matrix Theory, Interacting Particle Systems, and Integrable Systems*, 65,
  309-383. [arXiv:1205.2097](https://arxiv.org/abs/1205.2097)
