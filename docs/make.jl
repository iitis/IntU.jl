using IntU
using Documenter

DocMeta.setdocmeta!(IntU, :DocTestSetup, :(using IntU); recursive=true)

makedocs(;
    modules=[IntU],
    authors="USER",
    repo="https://github.com/iitis/IntU.jl/blob/{commit}{path}#{line}",
    sitename="IntU.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iitis.github.io/IntU.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/iitis/IntU.jl",
    devbranch="main",
)
