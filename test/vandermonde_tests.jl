
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

    y = OffsetArray(x, 0:n-1)
    C = vandermonde_matrix(y)
    D = vandermonde_matrix_inverse(y)

    @test C * D ≈ Matrix(I, n, n)
    @test inv(C) ≈ D
    @test inv(D) ≈ C

    @test A == C
    @test B == D

end
