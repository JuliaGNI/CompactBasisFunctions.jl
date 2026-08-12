
import CompactBasisFunctions: Basis, nodes, nnodes
import ContinuumArrays: grid

@testset "$(rpad("Basis Tests",80))" begin

    struct BasisTest{T} <: Basis{T} end

    test_basis = BasisTest{Float64}()

    @test_throws ErrorException basis(test_basis)
    @test_throws ErrorException nodes(test_basis)

    @test_throws ErrorException nbasis(test_basis)
    @test_throws ErrorException nnodes(test_basis)

    @test_throws ErrorException order(test_basis)
    @test_throws ErrorException degree(test_basis)

    @test_throws ErrorException eachindex(test_basis)

    # the message names the offending type, rather than just "Not implemented!"
    @test occursin("BasisTest", sprint(showerror, try basis(test_basis) catch e; e end))


    # the modal bases have no nodes, and say so. `grid` in particular used to fall through
    # to ContinuumArrays and raise a MethodError about `grid_axis`.
    for b in (Bernstein(4), Legendre(4))
        for f in (nodes, nnodes, grid)
            @test_throws ErrorException f(b)
            @test occursin("modal basis", sprint(showerror, try f(b) catch e; e end))
        end
    end

    # the nodal bases answer all three
    for b in (ChebyshevT(4), ChebyshevU(4), LagrangeGauß(4), LagrangeLobatto(4))
        @test nnodes(b) == 4
        @test length(nodes(b)) == 4
        @test grid(b) == nodes(b)
    end

end
