
import CompactBasisFunctions: Basis, nodes, nnodes, eachnode

@testset "$(rpad("Basis Tests",80))" begin

    struct BasisTest{T} <: Basis{T} end

    test_basis = BasisTest{Float64}()

    @test_throws ErrorException basis(test_basis)
    @test_throws ErrorException nodes(test_basis)

    @test_throws ErrorException nbasis(test_basis)
    @test_throws ErrorException nnodes(test_basis)

    @test_throws ErrorException eachbasis(test_basis)
    @test_throws ErrorException eachnode(test_basis)

    @test_throws ErrorException order(test_basis)
    @test_throws ErrorException degree(test_basis)

end
