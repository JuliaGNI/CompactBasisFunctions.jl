
import CompactBasisFunctions: Basis, nodes, nnodes
import ContinuumArrays
import ContinuumArrays: grid
import GeometricBase

@testset "$(rpad("Accessor bindings",80))" begin

    # the accessors must be methods on the shared generics, not functions of our own, or
    # loading CompactBasisFunctions together with another package of the ecosystem makes the
    # exported names resolve to nothing. Everything else in this file passes just as well
    # with five functions of our own, which is how the collision survived for so long.

    @test CompactBasisFunctions.basis === GeometricBase.basis
    @test CompactBasisFunctions.degree === GeometricBase.degree
    @test CompactBasisFunctions.nodes === GeometricBase.nodes
    @test CompactBasisFunctions.nnodes === GeometricBase.nnodes
    @test CompactBasisFunctions.order === GeometricBase.order

    @test CompactBasisFunctions.grid === ContinuumArrays.grid

    # nbasis is ours: nothing else in the ecosystem declares it
    @test parentmodule(CompactBasisFunctions.nbasis) == CompactBasisFunctions
end

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
    @test occursin("BasisTest", sprint(showerror, try
        basis(test_basis)
    catch e
        e
    end))

    # the modal bases have no nodes, and say so. `grid` in particular used to fall through
    # to ContinuumArrays and raise a MethodError about `grid_axis`.
    for b in (Bernstein(4), Legendre(4))
        for f in (nodes, nnodes, grid)
            @test_throws ErrorException f(b)
            @test occursin("modal basis", sprint(showerror, try
                f(b)
            catch e
                e
            end))
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
            Lagrange(BigFloat[0, 1 // 5, 1 // 2, 4 // 5, 1]))
            d = Derivative(axes(b, 1))

            for x in (BigFloat(1)/4, BigFloat(1)/3, BigFloat(7)/10), j in eachindex(b)

                fd = (b[x + h, j] - b[x - h, j]) / 2h
                @test isapprox(fd, (d * b)[x, j], atol = 1e-25)
            end
        end
    end

    # An evaluation runs in the wider of the basis's element type and the argument's, so a
    # basis of extended precision is not silently held to the precision of the point it is
    # asked about. `Legendre{BigFloat}` used to run Bonnet's recurrence in Float64 and widen
    # only the trailing sqrt(2j+1), returning a BigFloat accurate to Float64 and no further;
    # `Bernstein` and `Chebyshev` returned a Float64 outright, although their derivatives
    # already promoted. The assertion is that a Float64 argument gives exactly what its own
    # exact BigFloat value gives — i.e. that the argument's type limits nothing.
    setprecision(BigFloat, 256) do
        for b in (Bernstein(BigFloat, 6), Legendre(BigFloat, 6),
            Chebyshev{1}(BigFloat, 6), Chebyshev{2}(BigFloat, 6),
            Lagrange(BigFloat[0, 1 // 5, 1 // 2, 4 // 5, 1]))
            d = Derivative(axes(b, 1))

            # 0.1 is here because the shift onto [-1,1] must happen after the conversion and
            # not before: 2y-1 is exact in Float64 for y ≥ 0.25 by Sterbenz, so a test that
            # only asks about such points passes even when the shift still runs at the
            # argument's precision. At 0.1 it rounds, and the rounding is then frozen into
            # the widened value.
            for x in (0.1, 0.3, 0.75, 1.0), j in eachindex(b)

                @test b[x, j] isa BigFloat
                @test (d * b)[x, j] isa BigFloat

                @test b[x, j] == b[BigFloat(x), j]
                @test (d * b)[x, j] == (d * b)[BigFloat(x), j]

                # basis(b)[j](x) == b[x,j] is documented, so the promotion has to live in
                # the closures that `basis` hands out, not in `getindex`
                @test basis(b)[j](x) == b[x, j]
            end
        end
    end

    # ... and the promotion only ever widens: an argument more precise than the basis keeps
    # its own precision, as it did before
    setprecision(BigFloat, 256) do
        for b in (Bernstein(6), Legendre(6), ChebyshevT(6), ChebyshevU(6),
            Lagrange([0.0, 0.2, 0.5, 0.8, 1.0]))
            d = Derivative(axes(b, 1))

            @test b[BigFloat(3) / 10, 1] isa BigFloat
            @test b[0.3, 1] isa Float64

            @test (d * b)[BigFloat(3) / 10, 1] isa BigFloat
            @test (d * b)[0.3, 1] isa Float64
        end
    end
end
