
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


    # (d*b)[x,j] must really be the derivative of b[x,j], for every basis. This is checked
    # against a central difference in BigFloat, so the step can be small enough that the
    # comparison is sharp. The Legendre derivative used to drop the sqrt(2j+1) factor that
    # its basis functions carry, and was wrong by exactly that factor for every j > 0;
    # nothing related the two, so every direct assertion agreed with itself.
    setprecision(BigFloat, 256) do
        h = BigFloat(1) / 10^30

        for b in (Bernstein(BigFloat, 5), Legendre(BigFloat, 5),
                  Chebyshev{1}(BigFloat, 5), Chebyshev{2}(BigFloat, 5),
                  Lagrange(BigFloat[0, 1//5, 1//2, 4//5, 1]))
            d = Derivative(axes(b, 1))

            for x in (BigFloat(1)/4, BigFloat(1)/3, BigFloat(7)/10), j in eachindex(b)
                fd = (b[x+h, j] - b[x-h, j]) / 2h
                @test isapprox(fd, (d*b)[x, j], atol=1e-25)
            end
        end
    end

end
