using CompactBasisFunctions
using Documenter
using DocumenterCitations

# What the doctests need in scope. Shared with the `doctest` job of `.github/workflows/CI.yml`,
# which includes the same file, so a build and a doctest run cannot disagree.
include(joinpath(@__DIR__, "doctestsetup.jl"))

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    plugins = [bib],
    modules = [CompactBasisFunctions],
    authors = "Michael Kraus",
    repo = Remotes.GitHub("JuliaGNI", "CompactBasisFunctions.jl"),
    sitename = "CompactBasisFunctions.jl",
    # the internal helpers carry docstrings explaining the recurrences and the reasons
    # behind them; those belong at the call site, not in the manual
    checkdocs = :exports,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaGNI.github.io/CompactBasisFunctions.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Polynomial Approximation" => "approximation.md",
        "Basis Functions" => [
            "Lagrange" => "lagrange.md",
            "Chebyshev" => "chebyshev.md",
            "Legendre" => "legendre.md",
            "Bernstein" => "bernstein.md"
        ],
        "Usage" => "usage.md",
        "Library" => "library.md",
        "References" => "references.md"
    ]
)

deploydocs(;
    repo = "github.com/JuliaGNI/CompactBasisFunctions.jl",
    devurl = "latest",
    devbranch = "main"
)
