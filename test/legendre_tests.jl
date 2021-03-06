import ContinuumArrays: apply, MulQuasiMatrix

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

    @test (d*l)[0.0, 1] == 2.0
    @test (d*l)[0.5, 1] == 2.0
    @test (d*l)[1.0, 1] == 2.0
    @test (d*l)[2.0, 1] == 2.0

    @test (d*l)[0, 0] == 0.0
    @test (d*l)[1, 0] == 0.0
    @test (d*l)[2, 0] == 0.0

    @test (d*l)[0, 1] == 2.0
    @test (d*l)[1, 1] == 2.0
    @test (d*l)[2, 1] == 2.0


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

    @test (d*l)[0.0, 1] ==  2.0
    @test (d*l)[0.5, 1] ==  2.0
    @test (d*l)[1.0, 1] ==  2.0

    @test (d*l)[0.0, 2] == -6.0
    @test (d*l)[0.5, 2] ==  0.0
    @test (d*l)[1.0, 2] == +6.0

end
