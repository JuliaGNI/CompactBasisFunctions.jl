using CompactBasisFunctions
using Documenter

makedocs(;
    modules=[CompactBasisFunctions],
    authors="Michael Kraus",
    repo="https://github.com/JuliaGNI/CompactBasisFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="CompactBasisFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGNI.github.io/CompactBasisFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/CompactBasisFunctions.jl",
    devurl = "latest",
    devbranch = "main",
)
