import ContinuumArrays: apply, MulQuasiMatrix
import OffsetArrays: OffsetArray
import QuadratureRules: gauss_legendre_nodes

@testset "$(rpad("Lagrange Basis Tests",80))" begin
    x = [0.0, 1.0]
    l = Lagrange(x)
    d = Derivative(axes(l, 1))

    @test apply(*, d, l) isa LagrangeDerivative
    @test d*l isa LagrangeDerivative
    @test l' isa LagrangeDerivative

    @test grid(l) == x
    @test basis(l) == l.b
    @test nodes(l) == x
    @test nbasis(l) == 2
    @test nnodes(l) == 2
    @test eachindex(l) == 1:2
    @test order(l) == 2
    @test degree(l) == 1

    l0 = Lagrange(OffsetArray(x, 0:1))
    l1 = Lagrange([0.0, 1.0])
    l2 = Lagrange([0, 1])
    l3 = Lagrange([0.25, 0.75])
    l4 = Lagrange{Float32}(x)

    @test hash(l) == hash(l0)
    @test hash(l) == hash(l1)
    @test hash(l) == hash(l2)
    @test hash(l) != hash(l3)
    @test hash(l) == hash(l4)

    @test l == l0
    @test l == l1
    @test l == l2
    @test l != l3
    @test l == l4

    @test isequal(l, l0)
    @test isequal(l, l1)
    @test !isequal(l, l2)
    @test !isequal(l, l3)
    @test !isequal(l, l4)

    @test isapprox(l, l0)
    @test isapprox(l, l1)
    @test isapprox(l, l2)
    @test !isapprox(l, l3)
    @test isapprox(l, l4)

    @test l != LagrangeGauß(2)
    @test l == LagrangeLobatto(2)

    y = rand(5)

    z1 = [l[y[i], j] for i in eachindex(y), j in eachindex(l)]
    z2 = hcat([l[y[i], :] for i in eachindex(y)]...)'
    z3 = hcat([l[y, j] for j in eachindex(l)]...)
    z4 = l[y, :]

    @test z1 == z2 == z3 == z4

    z1 = [(d * l)[y[i], j] for i in eachindex(y), j in eachindex(l)]
    z2 = hcat([(d * l)[y[i], :] for i in eachindex(y)]...)'
    z3 = hcat([(d * l)[y, j] for j in eachindex(l)]...)
    z4 = (d * l)[y, :]

    @test z1 == z2 == z3 == z4

    @test l(0.0, 1) == 1.0
    @test l(0.5, 1) == 0.5
    @test l(1.0, 1) == 0.0

    @test l(0.0, 2) == 0.0
    @test l(0.5, 2) == 0.5
    @test l(1.0, 2) == 1.0

    @test l[0.0, 1] == 1.0
    @test l[0.5, 1] == 0.5
    @test l[1.0, 1] == 0.0

    @test l[0.0, 2] == 0.0
    @test l[0.5, 2] == 0.5
    @test l[1.0, 2] == 1.0

    @test (d * l)[0.0, 1] == -1.0
    @test (d * l)[0.5, 1] == -1.0
    @test (d * l)[1.0, 1] == -1.0

    @test (d * l)[0.0, 2] == 1.0
    @test (d * l)[0.5, 2] == 1.0
    @test (d * l)[1.0, 2] == 1.0

    @test (d * l)[1, 1] == -1.0
    @test (d * l)[2, 1] == -1.0

    @test (d * l)[1, 2] == 1.0
    @test (d * l)[2, 2] == 1.0

    # the internal buffers must be built in T, not in Float64, otherwise an
    # arbitrary-precision basis silently carries only Float64 precision
    @test eltype(Lagrange{Float32}([0.0f0, 0.25f0, 1.0f0]).diffs) == Float32

    setprecision(BigFloat, 256) do
        lb = Lagrange(gauss_legendre_nodes(BigFloat, 4))
        db = Derivative(axes(lb, 1))

        @test eltype(lb.denom) == BigFloat
        @test eltype(lb.diffs) == BigFloat

        z = BigFloat(1) / 7

        # partition of unity, and its derivative, hold to full BigFloat precision;
        # with Float64 buffers both deviate by ~1e-16
        @test abs(sum(lb[z, j] for j in eachindex(lb)) - 1) < 1e-70
        @test abs(sum((db * lb)[z, j] for j in eachindex(lb))) < 1e-70
    end

    # duplicate nodes used to give denom = [-4.0, -Inf, -Inf, 4.0] and evaluate to NaN
    @test_throws ArgumentError Lagrange([0.0, 0.5, 0.5, 1.0])
    @test_throws ArgumentError Lagrange([1.0, 1.0])
    @test_throws ArgumentError Lagrange{Float64}([0.0, 0.25, 0.25])

    # The check tests the product of the node differences, not the distinctness of the node
    # list. `allunique` compares with `isequal`, under which 0.0 and -0.0 are distinct — yet
    # their difference is zero, which is the whole point of the check, and the basis came out
    # with an Inf denominator. Conversely a lone NaN or Inf node is `isequal` to nothing at
    # all and passed straight through, poisoning every difference.
    @test_throws ArgumentError Lagrange([0.0, -0.0])
    @test_throws ArgumentError Lagrange([-0.0, 0.5, 0.0])
    @test_throws ArgumentError Lagrange([0.0, NaN, 1.0])
    @test_throws ArgumentError Lagrange([0.0, Inf])
    @test_throws ArgumentError Lagrange([0.0, 0.5, -Inf, 1.0])

    # repeated NaN was caught before, by `isequal(NaN, NaN)`, and still is
    @test_throws ArgumentError Lagrange([NaN, NaN])

    # ... while a legitimately signed zero among distinct nodes is fine
    @test nodes(Lagrange([-0.0, 0.5, 1.0])) == [-0.0, 0.5, 1.0]

    # A degenerate product does not on its own say which fault produced it: nodes that are
    # distinct and finite can still multiply out to zero or to an infinity, by underflow
    # when they are packed into a narrow range and by overflow when they are spread over a
    # wide one. Both used to give a silent Inf or 0 denominator; each is now rejected, and
    # the message says the product is at fault rather than sending the reader to look for a
    # repeated node that is not there.
    @test_throws ArgumentError Lagrange(Float32.(range(0, 0.01, length = 25)))
    @test_throws ArgumentError Lagrange([0.0, 1e160, 2e160, 3e160])

    for f in (() -> Lagrange(Float32.(range(0, 0.01, length = 25))),
        () -> Lagrange([0.0, 1e160, 2e160, 3e160]))
        @test !occursin("must be distinct", (try
            f()
        catch e
            e.msg
        end))
        @test !occursin("must be finite", (try
            f()
        catch e
            e.msg
        end))
    end

    # ... and a repeated node is still reported as one even when the surviving differences
    # overflow around it, so the diagnosis does not depend on which fault is noticed first
    @test occursin("must be distinct",
        (try
            Lagrange([0.0, 1e200, -1e200, 0.0])
        catch e
            e.msg
        end))

    # the same node sets are fine once they are scaled to where their product is
    # representable, which is what the message suggests
    @test nnodes(Lagrange(collect(range(0, 0.01, length = 25)))) == 25
    @test nnodes(Lagrange([0.0, 1.0, 2.0, 3.0])) == 4

    # ... and the cardinal property Lᵢ(xⱼ) = δᵢⱼ holds for the node sets that are accepted.
    # Written against positions rather than index values, since Lagrange numbers its basis
    # functions from 1 where the other three bases number theirs from 0.
    for l in (LagrangeGauß(4), LagrangeLobatto(4), Lagrange([0.0, 0.1, 0.7, 1.0]))
        xs = collect(nodes(l))
        idx = collect(eachindex(l))
        @test length(xs) == length(idx)

        for i in eachindex(idx), j in eachindex(xs)

            @test l[xs[j], idx[i]] ≈ (i == j ? 1 : 0) atol=1e-14
        end
    end
end
