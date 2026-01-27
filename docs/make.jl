using IntU
using Documenter

DocMeta.setdocmeta!(IntU, :DocTestSetup, :(using IntU); recursive=true)

makedocs(;
    modules=[IntU],
    authors="USER",
    repo=Documenter.Remotes.GitHub("iitis", "IntU.jl"),
    sitename="IntU.jl",
    remotes=nothing,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iitis.github.io/IntU.jl",
        assets=String[],
    ),
    pages=[

        "Home" => "index.md",
        "Manual" => [
            "Unitary Integration" => "unitary_integration.md",
            "Unitary Designs" => "unitary_designs.md",
            "Orthogonal & Symplectic" => "orthogonal_integration.md",

            "Gaussian Ensembles" => "gaussian_integration.md",
            "Pure States" => "pure_states.md",
            "Symbolic Traces" => "symbolic_trace.md",
            "Asymptotics" => "asymptotic.md",
            "QI Helpers" => "qi_helpers.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/iitis/IntU.jl",
    devbranch="main",
)
