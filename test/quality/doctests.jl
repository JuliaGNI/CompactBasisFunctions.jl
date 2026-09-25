using Documenter
using CompactBasisFunctions

# the same doctest setup as the documentation build and the Doctests job of CI
include(joinpath(@__DIR__, "..", "..", "docs", "doctestsetup.jl"))

# the manual pages resolve `CurrentModule = CompactBasisFunctions` in `Main`, and a
# `@safetestset` runs this file in a module of its own
@eval Main import CompactBasisFunctions

doctest(CompactBasisFunctions)
