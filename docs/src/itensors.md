# ITensors.jl Integration

IntU.jl provides a seamless bridge to [ITensors.jl](https://github.com/ITensors/ITensors.jl), allowing you to integrate entire tensor networks symbolically. This is particularly useful for studying random tensor networks, quantum circuits, and holographic models.

## Basic Usage

The primary way to use the integration is to wrap the tensors you wish to integrate in the `ITensorUnitary` struct.

```julia
using IntU, ITensors

# Define indices
i = Index(2, "Out")
j = Index(2, "In")

# Create a random ITensor
U_it = randomITensor(i, j)

# Wrap it to mark as a Haar-random unitary
U = ITensorUnitary(U_it; out_indices=[i], in_indices=[j])

# Integrate Tr(U A) over U(2)
A = randomITensor(j, i)
res = integrate([U, A], dU(2))
```

## Supported Measures

The ITensors integration supports all measure types provided by IntU:

| Measure | Usage |
| :--- | :--- |
| **Haar Unitary** | `dU(dim)` |
| **Orthogonal Group** | `dO(dim)` |
| **Symplectic Group** | `dSp(dim)` |
| **Unitary $t$-designs** | `dDesign(t, dim)` |

### Example: Orthogonal Integration
```julia
res = integrate([O1, O2, A], dO(3))
```

## Symbolic Dimensions

You can use symbolic dimensions from `Symbolics.jl` even when working with ITensors. While ITensor objects require a concrete size for construction, IntU will interpret the integration dimension symbolically if a symbolic variable is passed to the measure.

```julia
using Symbolics
@variables d_sym
res = integrate([U, U_dag], dU(d_sym))
# Result will be an ITensor with scalar value involving d_sym (e.g., 1/d_sym)
```

## Nested and Sequential Integration

For networks with multiple independent random unitaries, you can perform integration in steps.

```julia
# Network: U * V
# Step 1: Integrate over U
res_partial = integrate([U_wrap, V_it, A], dU(2))

# Step 2: Integrate the result over V
res_final = integrate([V_wrap, res_partial], dU(2))
```

## Performance & Algorithms

The ITensors integration uses a **Graphical Weingarten Engine**. Instead of
expanding the Full matrix, it:
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
