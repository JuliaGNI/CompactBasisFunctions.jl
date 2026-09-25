using CompactBasisFunctions
using Random
using Test

# a fixed seed, so that a failure is reproducible and an ill-conditioned draw cannot pass on
# one run and fail on the next
Random.seed!(0x6a1d3f2e)

import CompactBasisFunctions: vandermonde_matrix, vandermonde_matrix_inverse
import LinearAlgebra: I
import OffsetArrays: OffsetArray

@testset "$(rpad("Vandermonde Matrix",80))" begin
    n = 5
    x = rand(n)
    A = vandermonde_matrix(x)
    B = vandermonde_matrix_inverse(x)

    @test A * B ≈ Matrix(I, n, n)
    @test inv(A) ≈ B
    @test inv(B) ≈ A

    y = OffsetArray(x, 0:(n - 1))
    C = vandermonde_matrix(y)
    D = vandermonde_matrix_inverse(y)

    @test C * D ≈ Matrix(I, n, n)
    @test inv(C) ≈ D
    @test inv(D) ≈ C

    @test A == C
    @test B == D

    # the entries are the ones the docstrings pin down
    z = [0.0, 0.5, 1.0]

    @test vandermonde_matrix(z) == [1.0 0.0 0.0; 1.0 0.5 0.25; 1.0 1.0 1.0]
    @test vandermonde_matrix_inverse(z) == [1.0 0.0 0.0; -3.0 4.0 -1.0; 2.0 -4.0 2.0]

    # V maps monomial coefficients to values at the nodes, and its inverse maps back. This
    # is what the two matrices are for.
    c = [1.0, 0.0, 1.0]                                  # 1 + x²
    v = [1.0, 1.25, 2.0]                                 # its values at z

    @test vandermonde_matrix(z) * c == v
    @test vandermonde_matrix_inverse(z) * v == c
    @test vandermonde_matrix(z) \ v ≈ c

    # ... for any node set and any polynomial of degree < n
    coeffs = rand(n)
    values = [sum(coeffs[j] * xi^(j-1) for j in 1:n) for xi in x]

    @test A * coeffs ≈ values
    @test B * values ≈ coeffs

    # the element type is carried through: an exact arithmetic stays exact, and a narrow or
    # an extended one is neither widened nor silently computed in Float64
    r = [0 // 1, 1 // 2, 1 // 1]

    @test vandermonde_matrix(r) isa Matrix{Rational{Int}}
    @test vandermonde_matrix_inverse(r) isa Matrix{Rational{Int}}
    @test vandermonde_matrix(r) * vandermonde_matrix_inverse(r) == Matrix(I, 3, 3)

    @test eltype(vandermonde_matrix(Float32.(z))) == Float32
    @test eltype(vandermonde_matrix_inverse(Float32.(z))) == Float32

    setprecision(BigFloat, 256) do
        xb = BigFloat[0, 1 // 4, 1 // 2, 1]
        @test maximum(abs, vandermonde_matrix(xb) * vandermonde_matrix_inverse(xb) - I) <
              1e-70
    end

    # a single node interpolates by the constant polynomial, so both matrices are [1]
    @test vandermonde_matrix([0.7]) == ones(1, 1)
    @test vandermonde_matrix_inverse([0.7]) == ones(1, 1)
end
