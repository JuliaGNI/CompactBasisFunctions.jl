
import CompactBasisFunctions: Basis, nodes, nnodes
import ContinuumArrays
import ContinuumArrays: grid
import GeometricBase

"""
The allocations of the scalar surface of a basis: one evaluation, one derivative evaluation,
one application of `Derivative`, the counting accessors together and the comparisons together.

Measured inside a function so that the basis arrives with a concrete type. A closure over a
loop variable drawn from a heterogeneous tuple captures it as `Any`, and the dynamic call
then boxes its own result — sixteen bytes per call that the package never allocated.
"""
function scalar_allocations(b)
    d = Derivative(axes(b, 1))
    j = first(eachindex(b)) + 1

    evaluate() = b[0.3, j]
    evaluate_derivative() = (d * b)[0.3, j]
    differentiate() = Derivative(axes(b, 1)) * b
    counts() = nbasis(b) + order(b) + degree(b) + length(eachindex(b))
    compare() = hash(b) + UInt(b == b) + UInt(isequal(b, b)) + UInt(isapprox(b, b))

    # warm up: the first call in a process measures compilation
    evaluate(), evaluate_derivative(), differentiate(), counts(), compare()

    (evaluate = @allocated(evaluate()), derivative = @allocated(evaluate_derivative()),
        differentiate = @allocated(differentiate()), counts = @allocated(counts()),
        compare = @allocated(compare()))
end

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

    bases = (Bernstein(4), Legendre(4), ChebyshevT(4), ChebyshevU(4),
        LagrangeGauß(4), LagrangeLobatto(4))

    # the hierarchy is what carries the shared implementation: the accessors, the four
    # indexing forms, equality and the derivative product are each defined once, on
    # PolynomialBasis or on the nodal/modal split beneath it
    @test PolynomialBasis <: Basis
    @test NodalBasis <: PolynomialBasis
    @test ModalBasis <: PolynomialBasis

    @test Bernstein <: ModalBasis
    @test Legendre <: ModalBasis
    @test Chebyshev <: NodalBasis
    @test Lagrange <: NodalBasis

    for b in bases
        @test b isa PolynomialBasis
        @test b' isa PolynomialBasisDerivative
        @test Derivative(axes(b, 1)) * b isa PolynomialBasisDerivative
    end

    # ... and nothing is defined for a `Basis` this package does not own. Answering for one
    # would be type piracy — neither the generic, which is GeometricBase's, nor the argument
    # type, which is ContinuumArrays', belongs to this package — and it would claim to speak
    # for bases it knows nothing about.
    test_basis = BasisTest{Float64}()

    for f in (basis, nbasis, nodes, nnodes, order, degree)
        @test_throws MethodError f(test_basis)
    end

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

    # the three counts are one number under three names, for every basis
    for b in bases
        @test nbasis(b) == length(basis(b)) == length(eachindex(b)) == 4
        @test order(b) == nbasis(b)
        @test degree(b) == nbasis(b) - 1
        @test eltype(b) == Float64
    end

    # a basis is identified by the family it belongs to and the data that family is built
    # from, so two bases of different families are distinct even where that data agrees.
    # Nothing asserted this while each family carried its own `==`, which could only ever
    # be reached by two bases of the same family.
    same_nodes = Lagrange(nodes(ChebyshevU(3)))

    for (b1, b2) in ((Bernstein(4), Legendre(4)), (ChebyshevT(3), ChebyshevU(3)),
        (same_nodes, ChebyshevU(3)))
        @test b1 != b2
        @test !isequal(b1, b2)
        @test !isapprox(b1, b2)
        @test hash(b1) != hash(b2)
    end

    # a tolerance is a statement about nodes. A modal basis is built from a number of
    # functions, which nothing makes approximate, so it compares exactly however loose the
    # tolerance — forwarding the keywords there would make Bernstein(3) ≈ Bernstein(5).
    @test isapprox(ChebyshevT(Float32, 3), ChebyshevT(3))
    @test isapprox(Lagrange([0.0, 0.5, 1.0]), Lagrange([0.0, 0.5 + 1e-9, 1.0]), atol = 1e-6)
    @test !isapprox(Lagrange([0.0, 0.5, 1.0]), Lagrange([0.0, 0.5 + 1e-3, 1.0]), atol = 1e-6)

    @test isapprox(Bernstein(3), Bernstein(3), atol = 3)
    @test !isapprox(Bernstein(3), Bernstein(5), atol = 3)
    @test !isapprox(Legendre(3), Legendre(4), rtol = 0.5)

    # a function in a basis is the lazy product of the basis with its coefficients, which
    # materialises on indexing. This is the interface the README leads with, and the
    # coefficients are indexed as the basis is — a vector whose axis does not match is a
    # DimensionMismatch rather than a silent off-by-one.
    for b in bases
        c = [1 / (1 + j^2) for j in eachindex(b)]

        @test (b * c)[0.3] ≈ sum(c[j] * b[0.3, j] for j in eachindex(b))

        if first(eachindex(b)) == 0
            @test_throws DimensionMismatch b * collect(parent(c))
        end
    end

    # an index outside the basis is a BoundsError, for the basis and for its derivative
    # alike, in every family and in both of the indexing forms that take an index
    for b in bases
        d = Derivative(axes(b, 1))
        lo = first(eachindex(b)) - 1
        hi = last(eachindex(b)) + 1

        for j in (lo, hi)
            @test_throws BoundsError b[0.5, j]
            @test_throws BoundsError (d * b)[0.5, j]
            @test_throws BoundsError (d * b)[[0.25, 0.75], j]
        end
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

    # evaluation is type stable, which is what the promotion above must not cost: the
    # element type follows from the basis and the argument alone
    for b in bases
        d = Derivative(axes(b, 1))
        j = first(eachindex(b))

        @test @inferred(b[0.3, j]) isa Float64
        @test @inferred((d * b)[0.3, j]) isa Float64
    end

    # ... and it allocates nothing. `b[x,j]` and `(d*b)[x,j]` sit in the innermost loop of
    # every method that expands in a basis, and `Derivative(axes(b,1)) * b` is rebuilt there
    # too. An instability, or a closure that boxes what it captures, shows up here first and
    # changes nothing else. The array-returning forms `b[x,:]`, `b[X,j]` and `b[X,:]`
    # allocate their result and are not covered.
    #
    # Asserted only where bounds checking is left at `auto`. A run that forces
    # `--check-bounds=yes`, as `Pkg.test()` does on the older Julia versions of the CI
    # matrix, inflates allocation counts and would make any figure here meaningless.
    if Base.JLOptions().check_bounds == 0
        for b in bases
            a = scalar_allocations(b)

            @test a.evaluate == 0
            @test a.derivative == 0
            @test a.differentiate == 0
            @test a.counts == 0
            @test a.compare == 0
        end
    end
end
