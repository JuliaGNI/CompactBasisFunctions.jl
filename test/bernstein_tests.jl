import ContinuumArrays: apply, MulQuasiMatrix

@testset "$(rpad("Bernstein Basis Tests",80))" begin

    l = Bernstein(2)
    d = Derivative(axes(l,1))

    @test apply(*,d,l) isa BernsteinDerivative
    @test d*l isa BernsteinDerivative
    @test l' isa BernsteinDerivative

    @test basis(l) == l.b
    @test nbasis(l) == 2
    @test eachindex(l) == 0:1
    @test order(l) == 2
    @test degree(l) == 1


    l1 = Bernstein(Float64, 2)
    l2 = Bernstein(Integer, 2)
    l3 = Bernstein(Float64, 3)

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


    @test l(0.0, 0) == 1.0
    @test l(0.5, 0) == 0.5
    @test l(1.0, 0) == 0.0

    @test l(0.0, 1) == 0.0
    @test l(0.5, 1) == 0.5
    @test l(1.0, 1) == 1.0


    @test l[0.0, 0] == 1.0
    @test l[0.5, 0] == 0.5
    @test l[1.0, 0] == 0.0

    @test l[0.0, 1] == 0.0
    @test l[0.5, 1] == 0.5
    @test l[1.0, 1] == 1.0

    @test l[0, 0] == +1.0
    @test l[1, 0] ==  0.0
    @test l[2, 0] == -1.0

    @test l[0, 1] == 0.0
    @test l[1, 1] == 1.0
    @test l[2, 1] == 2.0


    @test (d*l)[0.0, 0] == -1.0
    @test (d*l)[0.5, 0] == -1.0
    @test (d*l)[1.0, 0] == -1.0

    @test (d*l)[0.0, 1] == +1.0
    @test (d*l)[0.5, 1] == +1.0
    @test (d*l)[1.0, 1] == +1.0

    @test (d*l)[0, 0] == -1.0
    @test (d*l)[1, 0] == -1.0
    @test (d*l)[2, 0] == -1.0

    @test (d*l)[0, 1] == +1.0
    @test (d*l)[1, 1] == +1.0
    @test (d*l)[2, 1] == +1.0


    # the closed form agrees with the two-index recurrence it replaced, evaluated in
    # BigFloat, and the out-of-range branch that the derivative relies on still holds
    function _bernstein_ref(j, p, x::BigFloat)
        (j < 0 || j > p) && return zero(BigFloat)
        p == 0 && return one(BigFloat)
        return _bernstein_ref(j, p-1, x) * (1-x) + _bernstein_ref(j-1, p-1, x) * x
    end

    setprecision(BigFloat, 256) do
        for p in 0:12, j in -1:p+1, x in (-0.7, 0.0, 1//3, 0.5, 1.0, 1.7)
            @test CompactBasisFunctions._bernstein(j, p, float(x)) ≈
                  Float64(_bernstein_ref(j, p, BigFloat(x))) rtol=1e-13 atol=1e-14
        end
    end

    # partition of unity, and its derivative, on every basis
    for n in 1:12
        b = Bernstein(n)
        db = Derivative(axes(b,1))
        for x in (0.0, 0.1, 0.5, 1//3, 1.0)
            @test sum(b[float(x), j] for j in eachindex(b)) ≈ 1
            @test abs(sum((db*b)[float(x), j] for j in eachindex(b))) < 1e-12
        end
    end

    # the recurrence used to descend into two subproblems per step, so a single value cost
    # O(2^p): 2.3 µs at n=10, 45 ms at n=25, and unusable beyond. The bound here is far
    # above what the closed form needs (~150 ns) and far below what recursion would take.
    b = Bernstein(60)
    b[0.3, 30]
    @test (@elapsed for _ in 1:100; b[0.3, 30]; end) < 1.0

end
