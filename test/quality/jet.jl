# JET.jl optimisation analysis of the hot path. See https://github.com/aviatesk/JET.jl.
#
# The entry points are the functions that `test/basis.jl` asserts with `@allocated`, at the
# concrete argument types that test passes: the six bases of `Float64`, a `Float64` point and
# an `Int` index.

using CompactBasisFunctions
using JET
using Test

import ContinuumArrays: Derivative

const JET_WORKS = isdefined(JET, :JET_AVAILABLE) ? JET.JET_AVAILABLE : JET.JET_LOADABLE

const BASES = ("Bernstein" => Bernstein(4), "Legendre" => Legendre(4),
    "ChebyshevT" => ChebyshevT(4), "ChebyshevU" => ChebyshevU(4),
    "LagrangeGauß" => LagrangeGauß(4), "LagrangeLobatto" => LagrangeLobatto(4))

if JET_WORKS
    @testset "$(rpad("JET report_opt $name", 80))" for (name, b) in BASES
        B = typeof(b)
        D = typeof(Derivative(axes(b, 1)))
        DB = typeof(Derivative(axes(b, 1)) * b)
        m = (CompactBasisFunctions,)

        @test isempty(JET.get_reports(JET.report_opt(getindex, (B, Float64, Int); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(getindex, (DB, Float64, Int); target_modules = m)))
        # the `@simplify` method is one call to `ArrayLayouts.Mul`, so a runtime dispatch it
        # causes is reported in that constructor's frame, which only `AnyFrameModule` keeps
        @test isempty(JET.get_reports(JET.report_opt(
            *, (D, B); target_modules = (JET.AnyFrameModule(CompactBasisFunctions),))))
        @test isempty(JET.get_reports(JET.report_opt(nbasis, (B,); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(order, (B,); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(degree, (B,); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(eachindex, (B,); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(hash, (B, UInt); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(==, (B, B); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(isequal, (B, B); target_modules = m)))
        @test isempty(JET.get_reports(JET.report_opt(isapprox, (B, B); target_modules = m)))
    end
else
    @test_skip "JET does not work on Julia $VERSION" # aviatesk/JET.jl#681
end
