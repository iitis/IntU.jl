using IntU
using Documenter

DocMeta.setdocmeta!(IntU, :DocTestSetup, :(using IntU); recursive = true)

makedocs(;
    modules = [IntU],
    authors = "USER",
    repo = Documenter.Remotes.GitHub("iitis", "IntU.jl"),
    sitename = "IntU.jl",
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://iitis.github.io/IntU.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Unitary Integration" => "unitary_integration.md",
            "Special Unitary Group" => "su_integration.md",
            "Diagonal Unitaries" => "diagonal_unitary.md",
            "Unitary Designs" => "unitary_designs.md",
            "Orthogonal & Symplectic" => "orthogonal_integration.md",
            "Permutation Groups" => "permutation_integration.md",
            "Circular Ensembles" => "circular_ensembles.md",
            "Gaussian Ensembles" => "gaussian_integration.md",
            "Pure States" => "pure_states.md",
            "Stiefel Manifolds" => "stiefel_manifold.md",
            "Symbolic Traces" => "symbolic_trace.md",
            "Asymptotics" => "asymptotic.md",
            "QI Helpers" => "qi_helpers.md",
            "ITensors Integration" => "itensors.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/iitis/IntU.jl", devbranch = "main")
