import ContinuumArrays: apply, MulQuasiMatrix
import LinearAlgebra: I
import QuadratureRules: GaussLegendreQuadrature, weights

@testset "$(rpad("Legendre Basis Tests",80))" begin

    l = Legendre(2)
    d = Derivative(axes(l,1))

    @test apply(*,d,l) isa LegendreDerivative
    @test d*l isa LegendreDerivative
    @test l' isa LegendreDerivative

    @test basis(l) == l.b
    @test nbasis(l) == 2
    @test eachindex(l) == 0:1
    @test order(l) == 2
    @test degree(l) == 1


    l1 = Legendre(Float64, 2)
    l2 = Legendre(Integer, 2)
    l3 = Legendre(Float64, 3)

    @test hash(l) == hash(l1)
    @test hash(l) == hash(l2)
    @test hash(l) != hash(l3)

    @test l == l1
    @test l == l2
    @test l != l3

    @test  isequal(l, l1)
    @test !isequal(l, l2)
    @test !isequal(l, l3)


    y = rand(5)

    z1 = parent([ l[y[i], j] for i in eachindex(y), j in eachindex(l)])
    z2 =   hcat([ l[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ l[y,    j] for j in eachindex(l)]...)
    z4 = parent(  l[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*l)[y[i], j] for i in eachindex(y), j in eachindex(l)])
    z2 =   hcat([ (d*l)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*l)[y,    j] for j in eachindex(l)]...)
    z4 = parent(  (d*l)[y,    :] )

    @test z1 == z2 == z3 == z4


    @test l(0.0, 0) ==  1.0
    @test l(0.5, 0) ==  1.0
    @test l(1.0, 0) ==  1.0

    @test l(0.0, 1) == -sqrt(3)
    @test l(0.5, 1) ==  0.0
    @test l(1.0, 1) == +sqrt(3)


    @test l[0.0, 0] ==  1.0
    @test l[0.5, 0] ==  1.0
    @test l[1.0, 0] ==  1.0

    @test l[0.0, 1] == -sqrt(3)
    @test l[0.5, 1] ==  0.0
    @test l[1.0, 1] == +sqrt(3)

    @test l[0, 0] ==  1.0
    @test l[1, 0] ==  1.0
    @test l[2, 0] ==  1.0

    @test l[0, 1] == -sqrt(3)
    @test l[1, 1] == +sqrt(3)
    @test l[2, 1] == +sqrt(27)


    @test (d*l)[0.0, 0] == 0.0
    @test (d*l)[0.5, 0] == 0.0
    @test (d*l)[1.0, 0] == 0.0
    @test (d*l)[2.0, 0] == 0.0

    @test (d*l)[0.0, 1] ==  2sqrt(3)
    @test (d*l)[0.5, 1] ==  2sqrt(3)
    @test (d*l)[1.0, 1] ==  2sqrt(3)
    @test (d*l)[2.0, 1] ==  2sqrt(3)

    @test (d*l)[0, 0] == 0.0
    @test (d*l)[1, 0] == 0.0
    @test (d*l)[2, 0] == 0.0

    @test (d*l)[0, 1] ==  2sqrt(3)
    @test (d*l)[1, 1] ==  2sqrt(3)
    @test (d*l)[2, 1] ==  2sqrt(3)


    l = Legendre(3)
    d = Derivative(axes(l,1))

    @test l[0.0, 0] == 1.0
    @test l[0.5, 0] == 1.0
    @test l[1.0, 0] == 1.0

    @test l(0.0, 1) == -sqrt(3)
    @test l(0.5, 1) ==  0.0
    @test l(1.0, 1) == +sqrt(3)

    @test l[0.0, 2] == +1.0 * sqrt(5)
    @test l[0.5, 2] == -0.5 * sqrt(5)
    @test l[1.0, 2] == +1.0 * sqrt(5)

    @test (d*l)[0.0, 0] ==  0.0
    @test (d*l)[0.5, 0] ==  0.0
    @test (d*l)[1.0, 0] ==  0.0

    @test (d*l)[0.0, 1] ==   2sqrt(3)
    @test (d*l)[0.5, 1] ==   2sqrt(3)
    @test (d*l)[1.0, 1] ==   2sqrt(3)

    @test (d*l)[0.0, 2] == -6sqrt(5)
    @test (d*l)[0.5, 2] ==  0.0
    @test (d*l)[1.0, 2] == +6sqrt(5)


    # the sqrt(2i+1) scaling makes the basis orthonormal on [0,1], i.e. the mass matrix
    # is the identity. Nothing asserted this, although it is the point of the scaling.
    quad = GaussLegendreQuadrature(24)

    for n in 1:8
        b = Legendre(n)
        M = [sum(weights(quad)[k] * b[nodes(quad)[k], i] * b[nodes(quad)[k], j]
                 for k in eachindex(nodes(quad))) for i in eachindex(b), j in eachindex(b)]
        @test parent(M) ≈ Matrix(I, n, n) atol=1e-13
    end

    # the recurrences used to descend into two subproblems per step, so a single value cost
    # O(φʲ): 571 µs at n=25 and unusable beyond. The bound is far above what iteration
    # needs (~470 ns at n=80) and far below what recursion would take. A liveness check
    # rather than a benchmark, cf. the same guard in bernstein_tests.jl.
    b = Legendre(80)
    db = Derivative(axes(b,1))
    b[0.3, 79]; (db*b)[0.3, 79]
    @test (@elapsed for _ in 1:100; b[0.3, 79]; end) < 1.0
    @test (@elapsed for _ in 1:100; (db*b)[0.3, 79]; end) < 1.0

end
