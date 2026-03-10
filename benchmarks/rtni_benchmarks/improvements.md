# Benchmark Results Comparison: IntU.jl vs RTNI

## 1. Result Representation Discrepancy

| Feature | IntU.jl (Element API) | RTNI (Multinomial API) |
| :--- | :--- | :--- |
| **Integrand Type** | Scalar element/Trace (e.g., `abs(U[1,1])^2`) | Matrix product (e.g., $E[U \otimes \bar{U}]$) |
| **Result Type** | Scalar number/symbolic fraction (e.g., `1/d`) | Matrix-valued expression/integrated graph |
| **Simplification** | Automatic simplification to lowest scalar form | Symbolic representation of tensor contraction |

### Why they look different:
The `IntU.jl` benchmarks use the high-level scalar API (`integrate(expr, measure)`), which is designed to return simplified analytical results. RTNI's `MultinomialexpectationvalueHaar` returns the "integrated graph," which includes symbolic projectors (like `x` for $E_{11}$). 

If you apply `Tr[]` to the RTNI result and simplify, it matches the `IntU.jl` scalar result exactly.

---

## 2. Graphical Integration Capabilities

`IntU.jl` is fully capable of the same graph-based integration as RTNI through its **Graphical Engine** (located in `src/itensors_integration.jl`).

### How to use it:
To obtain matrix-valued or tensor results similar to RTNI, use the `GraphicalUnitary` interface (requires `ITensors.jl`):

```julia
using ITensors, IntU

# Define open indices for a matrix-valued result
i, j, i_dag, j_dag = Index(d), Index(d), Index(d), Index(d)

# Define U and U_dagger
U = GraphicalUnitary([i], [j], false)
Ub = GraphicalUnitary([i_dag], [j_dag], true)

# Integrate for a matrix/tensor result
# Empty constants [] means we just integrate the unitaries themselves
res = integrate_graphical([], [U, Ub], dU(d))
# Result: An ITensor representing (1/d) * δ(i, i_dag) * δ(j, j_dag)
```

---

## 3. Calculating $|U_{11}|^{2k}$ as a Tensor Network

In a tensor network setting, $|U_{11}|^{2k}$ is represented as a **closed network** (vacuum diagram):

1.  **Nodes**: $k$ copies of $U$ and $k$ copies of $\bar{U}$.
2.  **Connections**: Plug output/input legs of $U$ and $\bar{U}$ into rank-1 projectors $E_{11} = |1\rangle\langle 1|$. 
3.  **Result**: Because the network has no open indices, the integration result is a rank-0 tensor—a **scalar number**.

RTNI's output for $|U_{11}|^2$ appears matrix-like only because it leaves the projector $E_{11}$ symbolic (represented by the variable `x`).
