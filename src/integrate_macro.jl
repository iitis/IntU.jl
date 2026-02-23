


const EXCLUDED_MACRO_SYMS = Set{Symbol}([
    :abs,
    :abs2,
    :real,
    :imag,
    :conj,
    :tr,
    :+,
    :-,
    :*,
    :/,
    :^,
    :getindex,
    :setindex!,
    :adjoint,
    :transpose,
    :I,
    :Number,
    :Complex,
    :Int,
    :Float64,
    :sum,
    :prod,
    :sqrt,
    :exp,
    :log,
    :sin,
    :cos,
    :tan,
    Symbol(":"),
])

"""
    @integrate expr measure

Symbolically integrate an expression over a measure, automatically declaring variables.
Heuristics are used to identify which symbols represent random matrices and dimensions.
"""
macro integrate(expr, measure)
    # 1. Parse measure to get dimension and type
    # measure is something like dU(d) or dO(2)
    if !Meta.isexpr(measure, :call)
        error("Measure must be a function call, e.g., dU(d)")
    end

    m_name = measure.args[1]
    m_dim = measure.args[2]

    # 2. Identify "natural" symbol and tag based on measure name
    natural_sym, tag = if m_name == :dU || m_name == :dCUE || m_name == :dSU
        :U, :U
    elseif m_name == :dO
        :O, :O
    elseif m_name == :dCOE
        :S, :COE
    elseif m_name == :dSp
        :Sp, :Sp
    elseif m_name == :dCSE
        :S, :CSE
    elseif m_name == :dPsi
        :psi, :psi
    elseif m_name == :dPerm
        :P, :Perm
    elseif m_name == :dCPerm
        :Y, :CPerm
    elseif m_name == :dDiagUnitary
        :D, :DiagUnitary
    elseif m_name == :dStiefel
        :V, :U
    elseif m_name == :dGUE || m_name == :dGOE || m_name == :dGSE
        :H, Symbol(string(m_name)[2:end])
    elseif m_name == :dGinUE || m_name == :dGinOE || m_name == :dGinSE
        :G, Symbol(string(m_name)[2:end])
    else
        :U, :U # default
    end

    # 3. Find all symbols in expr
    integrand_syms = Symbol[]
    MacroTools.postwalk(expr) do x
        if x isa Symbol
            push!(integrand_syms, x)
        end
        return x
    end

    # Exclude common functions and keywords
    integrand_syms = filter(s -> !(s in EXCLUDED_MACRO_SYMS), unique(integrand_syms))

    # 4. Generate declarations
    decls = []

    # Handle dimension
    if m_dim isa Symbol
        push!(decls, :(@variables $m_dim))
    end

    # Handle integrand symbols
    for s in unique(integrand_syms)
        if s == natural_sym
            push!(
                decls,
                :(
                    if !@isdefined($s) ||
                       !($s isa SymbolicMatrix) ||
                       $s.special_type !== $(QuoteNode(tag)) ||
                       !isequal($s.dim, $m_dim)
                        ;$s = SymbolicMatrix($(QuoteNode(s)), $(QuoteNode(tag)), $m_dim);
                    end
                ),
            )
        else
            push!(decls, :(
                if !@isdefined($s) || ($s isa SymbolicMatrix && $s.special_type !== :Constant)
                    ;$s = SymbolicMatrix($(QuoteNode(s)), :Constant, nothing);
                end
            ))
        end
    end

    return esc(quote
        $(decls...)
        integrate($expr, $measure)
    end)
end
