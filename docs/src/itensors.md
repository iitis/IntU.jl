# ITensors.jl Integration

IntU.jl provides a seamless bridge to [ITensors.jl](https://github.com/ITensors/ITensors.jl), allowing you to integrate entire tensor networks symbolically. This is particularly useful for studying random tensor networks, quantum circuits, and holographic models.

## Basic Usage

The primary way to use the integration is to wrap the tensors you wish to integrate in the `ITensorUnitary` struct. For pure symbolic work, you don't even need a "dummy" ITensor object; you can just specify the indices.

```julia
using IntU, ITensors

# Define indices
i = Index(2, "Out")
j = Index(2, "In")
i2 = Index(2, "Out2")
j2 = Index(2, "In2")

# Mark a Haar-random unitary U and its adjoint U_dag
U = ITensorUnitary(out_indices=[i], in_indices=[j])
U_dag = ITensorUnitary(out_indices=[j2], in_indices=[i2], is_adj=true)

# Integrate a balanced network [U, A, U_dag, B] over U(2)
# Here A and B connect the dangling indices
A = randomITensor(j, j2) 
B = randomITensor(i2, i)
res = integrate([U, A, U_dag, B], dU(2))
# The result is a scalar ITensor from Weingarten delta contractions
# (for d=2 the overall prefactor is 1/2)
```

## Defining Random Unitaries

The `ITensorUnitary` struct is used to mark tensors for integration. There are two primary ways to create it, depending on whether you are working purely symbolically or with existing data.

### 1. Symbolic Placeholders (Recommended)
For pure symbolic integration, you only need to specify the indices. No "dummy" ITensor object is required.

```julia
# Mark a Haar-random unitary by its indices
U = ITensorUnitary(out_indices=[i], in_indices=[j])
```

> [!TIP]
> This is usually the cleanest approach for deriving analytical formulas or performing benchmarks where the specific tensor values don't matter.

### 2. Wrapping Existing Tensors
If you already have an `ITensor` (e.g., from a numeric simulation) and want to integrate over its values as if it were random, you can wrap the existing object.

```julia
U_it = randomITensor(i, j)
U = ITensorUnitary(U_it; out_indices=[i], in_indices=[j])
```

In this case, the `integrate` function will ignore the current numerical contents of `U_it` and treat it as a random placeholder. After integration, the result will correctly contract with any other tensors connected to those indices.

## Supported Measures

The ITensors integration supports the following measure types:

| Measure | Usage |
| :--- | :--- |
| **Haar Unitary** | `dU(dim)` |
| **Orthogonal Group** | `dO(dim)` |
| **Symplectic Group** | `dSp(dim)` |
| **Unitary $t$-designs** | `dDesign(dim, t)` |

### Example: Orthogonal Integration
```julia
using IntU, ITensors

# Define indices for a 3x3 orthogonal matrix O
o1 = Index(3, "Out")
i1 = Index(3, "In")

# Mark O as an orthogonal random matrix
O = ITensorUnitary(out_indices=[o1], in_indices=[i1])

# Constant tensor A
A = randomITensor(o1, i1)

# Integrate Tr(O A) - should be zero as E[O_ij] = 0
res = integrate([O, A], dO(3))
```

## Symbolic Dimensions

You can use symbolic dimensions from `Symbolics.jl` even when working with ITensors. While ITensor objects require a concrete size for construction, IntU will interpret the integration dimension symbolically if a symbolic variable is passed to the measure.

```julia
using IntU, ITensors, Symbolics

@variables d_sym
i = Index(2, "Out"); j = Index(2, "In")
i2 = Index(2, "Out2"); j2 = Index(2, "In2")

U = ITensorUnitary(out_indices=[i], in_indices=[j])
U_dag = ITensorUnitary(out_indices=[j2], in_indices=[i2], is_adj=true)

# Integrate E[U ⊗ U†]
res = integrate([U, U_dag], dU(d_sym))
# Result is an ITensor with indices {i, j, i2, j2} involving d_sym
```

## Nested and Sequential Integration

For networks with multiple independent random unitaries, you can perform integration in steps.

```julia
using IntU, ITensors

i = Index(2); j = Index(2); k = Index(2)
U_it = randomITensor(i, j); V_it = randomITensor(j, k)
A = randomITensor(k, i) # Constant tensor to close the loop

# Wrap U and V
U_wrap = ITensorUnitary(U_it, out_indices=[i], in_indices=[j])
V_wrap = ITensorUnitary(V_it, out_indices=[j], in_indices=[k])

# Network: Tr(U * V * A)
# Step 1: Integrate over U
res_partial = integrate([U_wrap, V_it, A], dU(2))

# Step 2: Integrate the result over V
res_final = integrate([V_wrap, res_partial], dU(2))
```

## Performance & Algorithms

The ITensors integration uses a **Graphical Weingarten Engine**. Instead of
expanding the full matrix, it:
1.  Identifies the network topology (which indices are connected to which tensors).
2.  Generates the possible Weingarten contractions as deltas over original ITensor indices.
3.  Returns a sum of contracted ITensor networks.

This is much more efficient than traditional matrix methods for large sparse networks.

## API Reference

```@docs
ITensorUnitary
GraphicalUnitary
IntU.integrate_graphical
```
