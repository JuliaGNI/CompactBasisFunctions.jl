using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Vandermonde matrix" include("vandermonde_matrix.jl")
    @safetestset "Basis" include("basis.jl")
    @safetestset "Bernstein" include("bernstein.jl")
    @safetestset "Chebyshev" include("chebyshev.jl")
    @safetestset "Lagrange" include("lagrange.jl")
    @safetestset "Legendre" include("legendre.jl")
end
if "slow" in GROUPS
    @safetestset "Doctests" include("quality/doctests.jl")
end
