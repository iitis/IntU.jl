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
