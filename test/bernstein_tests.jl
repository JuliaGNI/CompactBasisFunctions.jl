import ContinuumArrays: apply, MulQuasiMatrix

@testset "$(rpad("Bernstein Basis Tests",80))" begin

    l = Bernstein(2)
    d = Derivative(axes(l,1))

    @test apply(*,d,l) isa BernsteinDerivative
    @test d*l isa BernsteinDerivative

    @test basis(l) == l.b
    @test nbasis(l) == 2
    @test eachbasis(l) == 0:1
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

    z1 = parent([ l[y[i], j] for i in eachindex(y), j in eachbasis(l)])
    z2 =   hcat([ l[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ l[y,    j] for j in eachbasis(l)]...)
    z4 = parent(  l[y,    :] )

    @test z1 == z2 == z3 == z4

    z1 = parent([ (d*l)[y[i], j] for i in eachindex(y), j in eachbasis(l)])
    z2 =   hcat([ (d*l)[y[i], :] for i in eachindex(y)]...)'
    z3 =   hcat([ (d*l)[y,    j] for j in eachbasis(l)]...)
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


    @test (d*l)[0.0, 0] == -1.0
    @test (d*l)[0.5, 0] == -1.0
    @test (d*l)[1.0, 0] == -1.0

    @test (d*l)[0.0, 1] == 1.0
    @test (d*l)[0.5, 1] == 1.0
    @test (d*l)[1.0, 1] == 1.0

    @test (d*l)[1, 0] == -1.0
    @test (d*l)[2, 0] == -1.0

    @test (d*l)[1, 1] == 1.0
    @test (d*l)[2, 1] == 1.0

end
