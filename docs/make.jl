using IntegrateUnitary
using Documenter

DocMeta.setdocmeta!(IntegrateUnitary, :DocTestSetup, :(using IntegrateUnitary); recursive = true)

makedocs(;
    modules = [IntegrateUnitary],
    authors = "Łukasz Pawela and Zbigniew Puchała",
    repo = Documenter.Remotes.GitHub("iitis", "IntegrateUnitary.jl"),
    sitename = "IntegrateUnitary.jl",
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://iitis.github.io/IntegrateUnitary.jl",
        assets = String[],
    ),
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Unitary Integration" => "unitary_integration.md",
            "Diagonal Unitaries" => "diagonal_unitary.md",
            "Orthogonal & Symplectic" => "orthogonal_integration.md",
            "Permutation Groups" => "permutation_integration.md",
            "Circular Ensembles" => "circular_ensembles.md",
            "Gaussian Ensembles" => "gaussian_integration.md",
            "Pure States" => "pure_states.md",
            "Stiefel Manifolds" => "stiefel_manifold.md",
            "Symbolic Trace Integration" => "symbolic_trace.md",
            "Asymptotics" => "asymptotic.md",
            "Integral Library" => "integral_library.md",
            "QI Helpers" => "qi_helpers.md",
            "ITensors Integration" => "itensors.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/iitis/IntegrateUnitary.jl", devbranch = "main")
