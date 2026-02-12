# src/macros.jl

"""
    @integrate expr measure

Convenience macro for performing integration.
Equivalent to `integrate(expr, measure)`.

Example:
```julia
@variables d
U = SymbolicMatrix(:U)
@integrate abs(U[1,1])^2 dU(d)
```
"""
macro integrate(expr, measure)
    return quote
        integrate($(esc(expr)), $(esc(measure)))
    end
end

"""
    @symbolic_dimension expr

Legacy macro for defining a symbolic matrix with a symbolic dimension.
Example: `@symbolic_dimension V[1:d, 1:d]`
"""
macro symbolic_dimension(expr)
    if Meta.isexpr(expr, :ref)
        name = expr.args[1]
        indices = expr.args[2:end]
        
        # Extract dimension from something like 1:d
        dim = nothing
        if length(indices) >= 1 && Meta.isexpr(indices[1], :call) && indices[1].args[1] == :(:)
            dim = indices[1].args[3]
        end
        
        return quote
            $(esc(name)) = SymbolicMatrix($(QuoteNode(name)), :U, $(esc(dim)))
        end
    end
    error("Invalid usage of @symbolic_dimension. Expected format: @symbolic_dimension V[1:d, 1:d]")
end
