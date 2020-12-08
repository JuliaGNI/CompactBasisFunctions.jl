import ContinuumArrays: apply, MulQuasiMatrix

@testset "$(rpad("Lagrange Basis Tests",80))" begin

    x = [0.0, 1.0]
    l = Lagrange(x)
    d = Derivative(axes(l,1))

    @test apply(*,d,l) isa LagrangeDerivative
    @test d*l isa LagrangeDerivative

    @test grid(l) == x
    @test nodes(l) == x
    @test nnodes(l) == 2
    @test nbasis(l) == 2
    @test degree(l) == 1
    @test eachnode(l) == 1:2


    l1 = Lagrange([0.0,  1.0 ])
    l2 = Lagrange([0,    1   ])
    l3 = Lagrange([0.25, 0.75])

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

    z1 =      [ l[y[i], j] for i in eachindex(y), j in eachnode(l) ]
    z2 = hcat([ l[y[i], :] for i in eachindex(y)]...)'
    z3 = hcat([ l[y,    j] for j in eachnode(l) ]...)
    z4 =        l[y,    :]

    @test z1 == z2 == z3 == z4

    z1 =      [ (d*l)[y[i], j] for i in eachindex(y), j in eachnode(l) ]
    z2 = hcat([ (d*l)[y[i], :] for i in eachindex(y)]...)'
    z3 = hcat([ (d*l)[y,    j] for j in eachnode(l) ]...)
    z4 =        (d*l)[y,    :]

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


    @test (d*l)[0.0, 1] == -1.0
    @test (d*l)[0.5, 1] == -1.0
    @test (d*l)[1.0, 1] == -1.0

    @test (d*l)[0.0, 2] == 1.0
    @test (d*l)[0.5, 2] == 1.0
    @test (d*l)[1.0, 2] == 1.0

    @test (d*l)[1, 1] == -1.0
    @test (d*l)[2, 1] == -1.0

    @test (d*l)[1, 2] == 1.0
    @test (d*l)[2, 2] == 1.0

end
