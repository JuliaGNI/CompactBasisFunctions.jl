import ContinuumArrays: apply, MulQuasiMatrix
import OffsetArrays: OffsetArray

@testset "$(rpad("Chebyshev Basis Tests",80))" begin

    t = ChebyshevT(2)
    u = ChebyshevU(2)
    d = Derivative(axes(t,1))

    @test t isa Chebyshev
    @test u isa Chebyshev

    @test ChebyshevTDerivative <: ChebyshevDerivative
    @test ChebyshevUDerivative <: ChebyshevDerivative

    @test apply(*,d,t) isa ChebyshevTDerivative
    @test apply(*,d,t) isa ChebyshevDerivative
    @test d*t isa ChebyshevTDerivative
    @test d*t isa ChebyshevDerivative
    @test t' isa ChebyshevTDerivative
    @test t' isa ChebyshevDerivative

    @test apply(*,d,u) isa ChebyshevUDerivative
    @test apply(*,d,u) isa ChebyshevDerivative
    @test d*u isa ChebyshevUDerivative
    @test d*u isa ChebyshevDerivative
    @test u' isa ChebyshevUDerivative
    @test u' isa ChebyshevDerivative

    @test grid(t) == t.x
    @test basis(t) == t.b
    @test nodes(t) == t.x
    @test nbasis(t) == 2
    @test nnodes(t) == 2
    @test eachindex(t) == 0:1
    @test order(t)  == 2
    @test degree(t) == 1

    @test grid(u) == u.x
    @test basis(u) == u.b
    @test nodes(u) == u.x
    @test nbasis(u) == 2
    @test nnodes(u) == 2
    @test eachindex(u) == 0:1
    @test order(u)  == 2
    @test degree(u) == 1


    t0 = Chebyshev(Float64, 2, Val(1))
    t1 = ChebyshevT(Float64, 2)
    t2 = ChebyshevT(Float32, 2)
    t3 = ChebyshevT(Float64, 3)
    t4 = Chebyshev(2, Val(1))

    @test hash(t) == hash(t0)
    @test hash(t) == hash(t1)
    @test hash(t) != hash(t2)
    @test hash(t) != hash(t3)
    @test hash(t) == hash(t4)

    @test t == t0
    @test t == t1
    @test t ≈  t2
    @test t ≠  t3
    @test t == t4

    @test  isequal(t, t0)
    @test  isequal(t, t1)
    @test !isequal(t, t2)
    @test !isequal(t, t3)
    @test  isequal(t, t4)


    u0 = Chebyshev(Float64, 2, Val(2))
    u1 = ChebyshevU(Float64, 2)
    u2 = ChebyshevU(Integer, 2)
    u3 = ChebyshevU(Float64, 3)
    u4 = Chebyshev(2, Val(2))

    @test hash(u) == hash(u0)
    @test hash(u) == hash(u1)
    @test hash(u) == hash(u2)
    @test hash(u) != hash(u3)
    @test hash(u) == hash(u4)

    @test u == u0
    @test u == u1
    @test u ≈  u2
    @test u ≠  u3
    @test u == u4

    @test  isequal(u, u0)
    @test  isequal(u, u1)
    @test !isequal(u, u2)
    @test !isequal(u, u3)
    @test  isequal(u, u4)


    @test hash(t) != hash(u)


    y = rand(5)

    z1 = parent([ t[y[i], j] for i in eachindex(y), j in eachindex(t)])
    z2 =   hcat([ t[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ t[y,    j] for j in eachindex(t)]...)
    z4 = parent(  t[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ u[y[i], j] for i in eachindex(y), j in eachindex(u)] )
    z2 =   hcat([ u[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ u[y,    j] for j in eachindex(u)]...)
    z4 = parent(  u[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*t)[y[i], j] for i in eachindex(y), j in eachindex(t)])
    z2 =   hcat([ (d*t)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*t)[y,    j] for j in eachindex(t)]...)
    z4 = parent(  (d*t)[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*u)[y[i], j] for i in eachindex(y), j in eachindex(u)])
    z2 =   hcat([ (d*u)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*u)[y,    j] for j in eachindex(u)]...)
    z4 = parent(  (d*u)[y,    :] )

    @test z1 == z2 == z3 == z4


    # the nodes live on [0,1], are ascending, and lie inside axes(C,1)
    # (ChebyshevU is not defined for a single node, hence n ≥ 2)
    for kind in (1, 2), n in 2:8
        C = Chebyshev{kind}(n)
        @test nnodes(C) == n
        @test all(x -> x ∈ axes(C, 1), grid(C))
        @test issorted(nodes(C))
    end

    @test nodes(ChebyshevT(2)) ≈ [(1 - sqrt(2)/2) / 2, (1 + sqrt(2)/2) / 2]
    @test nodes(ChebyshevU(2)) == [0.0, 1.0]


    # basis and derivative values in closed form, with x̃ = 2x-1
    t = ChebyshevT(3)
    u = ChebyshevU(3)
    d = Derivative(axes(t,1))

    @test t[0.0, 0] == 1.0
    @test t[0.5, 0] == 1.0
    @test t[1.0, 0] == 1.0

    @test t[0.0, 1] == -1.0
    @test t[0.5, 1] ==  0.0
    @test t[1.0, 1] == +1.0

    @test t[0.0, 2] == +1.0
    @test t[0.5, 2] == -1.0
    @test t[1.0, 2] == +1.0

    # d/dx Tᵢ(2x-1) = 2i U_{i-1}(2x-1), so the i=2 branch is 16x-8
    @test (d*t)[0.0, 0] == 0.0
    @test (d*t)[0.5, 0] == 0.0
    @test (d*t)[1.0, 0] == 0.0

    @test (d*t)[0.0, 1] == 2.0
    @test (d*t)[0.5, 1] == 2.0
    @test (d*t)[1.0, 1] == 2.0

    @test (d*t)[0.0, 2] == -8.0
    @test (d*t)[0.5, 2] ==  0.0
    @test (d*t)[1.0, 2] == +8.0

    @test u[0.0, 0] == 1.0
    @test u[0.5, 0] == 1.0
    @test u[1.0, 0] == 1.0

    @test u[0.0, 1] == -2.0
    @test u[0.5, 1] ==  0.0
    @test u[1.0, 1] == +2.0

    @test u[0.0, 2] == +3.0
    @test u[0.5, 2] == -1.0
    @test u[1.0, 2] == +3.0

    @test (d*u)[0.25, 0] == 0.0
    @test (d*u)[0.50, 0] == 0.0
    @test (d*u)[0.75, 0] == 0.0

    @test (d*u)[0.25, 1] == 4.0
    @test (d*u)[0.50, 1] == 4.0
    @test (d*u)[0.75, 1] == 4.0

    @test (d*u)[0.25, 2] == -8.0
    @test (d*u)[0.50, 2] ==  0.0
    @test (d*u)[0.75, 2] == +8.0

    # the endpoints x = 0 and x = 1 are where x̃ = ∓1, at which the closed form
    # U'ᵢ = ((i+1)Tᵢ₊₁ - x̃ Uᵢ) / (x̃² - 1) is 0/0 and used to return NaN. The
    # derivative is a polynomial and finite there, with U'ᵢ(±1) known in closed form.
    u = ChebyshevU(6)
    d = Derivative(axes(u,1))

    @test [(d*u)[1.0, i] for i in eachindex(u)] == [2 * i*(i+1)*(i+2)//3 for i in eachindex(u)]
    @test parent([(d*u)[1.0, i] for i in eachindex(u)]) == [0, 4, 16, 40, 80, 140]

    @test [(d*u)[0.0, i] for i in eachindex(u)] == [2 * (-1)^(i+1) * i*(i+1)*(i+2)//3 for i in eachindex(u)]
    @test parent([(d*u)[0.0, i] for i in eachindex(u)]) == [0, 4, -16, 40, -80, 140]

    # ... and they are the limits of the interior values, not a separate branch
    for i in eachindex(u)
        @test (d*u)[1.0 - 1e-8, i] ≈ (d*u)[1.0, i] atol=1e-4
        @test (d*u)[0.0 + 1e-8, i] ≈ (d*u)[0.0, i] atol=1e-4
    end

    # first-kind derivatives are finite at the endpoints too: T'ᵢ(±1) = ±... i²
    t = ChebyshevT(6)
    dt = Derivative(axes(t,1))

    @test [(dt*t)[1.0, i] for i in eachindex(t)] == [2 * i^2 for i in eachindex(t)]
    @test [(dt*t)[0.0, i] for i in eachindex(t)] == [2 * (-1)^(i+1) * i^2 for i in eachindex(t)]

    # out-of-range basis indices are a BoundsError, as for the basis itself
    @test_throws BoundsError (d*u)[0.5, -1]
    @test_throws BoundsError (d*u)[0.5, nbasis(u)]


    # the recurrences used to descend into two subproblems per step, so a single value cost
    # O(φʲ): 323 µs at n=25 and unusable beyond. The bound is far above what iteration
    # needs (~200 ns at n=80) and far below what recursion would take. A liveness check
    # rather than a benchmark, cf. the same guard in bernstein_tests.jl.
    for kind in (1, 2)
        b = Chebyshev{kind}(80)
        db = Derivative(axes(b,1))
        b[0.3, 79]; (db*b)[0.3, 79]
        @test (@elapsed for _ in 1:100; b[0.3, 79]; end) < 1.0
        @test (@elapsed for _ in 1:100; (db*b)[0.3, 79]; end) < 1.0
    end

end
