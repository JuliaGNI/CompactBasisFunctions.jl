import ContinuumArrays: apply, MulQuasiMatrix

@testset "$(rpad("Lagrange Basis Tests",80))" begin

    x = [0.0, 1.0]
    l = Lagrange(x)
    d = Derivative(axes(l,1))

    @test apply(*,d,l) isa LagrangeDerivative
    @test d*l isa LagrangeDerivative

    @test nodes(l) == x
    @test nnodes(l) == 2
    @test nbasis(l) == 2
    @test degree(l) == 1
    @test eachnode(l) == 1:2


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
