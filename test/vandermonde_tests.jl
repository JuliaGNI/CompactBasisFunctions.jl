
import CompactBasisFunctions: vandermonde_matrix, vandermonde_matrix_inverse
import LinearAlgebra: I

@testset "$(rpad("Vandermonde Matrix",80))" begin

    n = 5
    x = rand(n)
    A = vandermonde_matrix(x)
    B = vandermonde_matrix_inverse(x)

    @test A * B ≈ Matrix(I, n, n)
    @test inv(A) ≈ B
    @test inv(B) ≈ A

end
