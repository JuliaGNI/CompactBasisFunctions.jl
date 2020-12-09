import ContinuumArrays: apply, MulQuasiMatrix
import FastTransforms: chebyshevpoints
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
    @test eachnode(t) == 1:2
    @test eachbasis(t) == 0:1
    @test order(t)  == 2
    @test degree(t) == 1

    @test grid(u) == u.x
    @test basis(u) == u.b
    @test nodes(u) == u.x
    @test nbasis(u) == 2
    @test nnodes(u) == 2
    @test eachnode(u) == 1:2
    @test eachbasis(u) == 0:1
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

    z1 = parent([ t[y[i], j] for i in eachindex(y), j in eachbasis(t)])
    z2 =   hcat([ t[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ t[y,    j] for j in eachbasis(t)]...)
    z4 = parent(  t[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ u[y[i], j] for i in eachindex(y), j in eachbasis(u)] )
    z2 =   hcat([ u[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ u[y,    j] for j in eachbasis(u)]...)
    z4 = parent(  u[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*t)[y[i], j] for i in eachindex(y), j in eachbasis(t)])
    z2 =   hcat([ (d*t)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*t)[y,    j] for j in eachbasis(t)]...)
    z4 = parent(  (d*t)[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*u)[y[i], j] for i in eachindex(y), j in eachbasis(u)])
    z2 =   hcat([ (d*u)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*u)[y,    j] for j in eachbasis(u)]...)
    z4 = parent(  (d*u)[y,    :] )

    @test z1 == z2 == z3 == z4

end
