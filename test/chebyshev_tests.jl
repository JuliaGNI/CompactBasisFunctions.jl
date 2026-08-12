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

    # d/dx Uᵢ(2x-1) is evaluated via a formula that is singular at x̃ = ±1,
    # i.e. at x = 0 and x = 1, so only interior points are checked here
    @test (d*u)[0.25, 0] == 0.0
    @test (d*u)[0.50, 0] == 0.0
    @test (d*u)[0.75, 0] == 0.0

    @test (d*u)[0.25, 1] == 4.0
    @test (d*u)[0.50, 1] == 4.0
    @test (d*u)[0.75, 1] == 4.0

    @test (d*u)[0.25, 2] == -8.0
    @test (d*u)[0.50, 2] ==  0.0
    @test (d*u)[0.75, 2] == +8.0

end
