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
    :kron,
    Symbol(":"),
])

"""
    @integrate expr measure

Symbolically integrate an expression over a measure, automatically declaring variables.
Heuristics are used to identify which symbols represent random matrices and dimensions.
"""
macro integrate(expr, measure)
    if !Meta.isexpr(measure, :call)
        throw(ArgumentError("Measure must be a function call, e.g., dU(d)"))
    end

    m_name = measure.args[1]
    m_dim = measure.args[2]

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
        :V, :V
    elseif m_name == :dGUE || m_name == :dGOE || m_name == :dGSE
        :H, Symbol(string(m_name)[2:end])
    elseif m_name == :dGinUE || m_name == :dGinOE || m_name == :dGinSE
        :G, Symbol(string(m_name)[2:end])
    else
        :U, :U # default
    end

    integrand_syms = Symbol[]
    MacroTools.postwalk(expr) do x
        if x isa Symbol
            push!(integrand_syms, x)
        end
        return x
    end

    integrand_syms = filter(s -> !(s in EXCLUDED_MACRO_SYMS), unique(integrand_syms))

    decls = []

    for arg in measure.args[2:end]
        if arg isa Symbol
            push!(decls, :(
                if !@isdefined($arg)
                    @variables $arg
                end
            ))
        end
    end

    for s in unique(integrand_syms)
        if s == natural_sym
            push!(
                decls,
                :(
                    if !@isdefined($s) ||
                       !($s isa SymbolicMatrix) ||
                       $s.special_type !== $(QuoteNode(tag)) ||
                       !isequal($s.dim, $m_dim)
                        dim = if $(QuoteNode(tag)) === :psi
                            (IntegrateUnitary._ensure_symbolic_dim($m_dim), 1)
                        elseif $(QuoteNode(tag)) === :V
                            m_args = $(measure.args)
                            (
                                IntegrateUnitary._ensure_symbolic_dim(m_args[2]),
                                IntegrateUnitary._ensure_symbolic_dim(m_args[3]),
                            )
                        else
                            IntegrateUnitary._ensure_symbolic_dim($m_dim)
                        end
                        $s = SymbolicMatrix($(QuoteNode(s)), $(QuoteNode(tag)), dim);
                    end
                ),
            )
        else
            push!(
                decls,
                :(
                    if !@isdefined($s) ||
                       ($s isa SymbolicMatrix && $s.special_type !== :Constant)
                        ;$s = SymbolicMatrix($(QuoteNode(s)), :Constant, nothing);
                    end
                ),
            )
        end
    end

    return esc(quote
        $(decls...)
        integrate($expr, $measure)
    end)
end
