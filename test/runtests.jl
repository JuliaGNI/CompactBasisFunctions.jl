using CompactBasisFunctions
using Random
using Test

# several tests draw their evaluation points with rand(); seed so that a failure is
# reproducible, and so that the Vandermonde tests cannot hit an ill-conditioned draw
# on one run and pass on the next
Random.seed!(0x6a1d3f2e)

include("vandermonde_tests.jl")
include("basis_tests.jl")
include("bernstein_tests.jl")
include("chebyshev_tests.jl")
include("lagrange_tests.jl")
include("legendre_tests.jl")
