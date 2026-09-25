using Documenter
using CompactBasisFunctions

# the same doctest setup as the documentation build and the Doctests job of CI
include(joinpath(@__DIR__, "..", "..", "docs", "doctestsetup.jl"))

doctest(CompactBasisFunctions; manual = false)
