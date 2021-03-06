module CompactBasisFunctions

    using ContinuumArrays
    using StaticArrays

    import Base: *
    import ContinuumArrays: Mul, QMul2, QuasiAdjoint, (..), @simplify

    export Derivative, ℵ₁

    include("vandermonde_matrix.jl")

    export Basis, basis, nbasis, grid, degree, order

    include("basis.jl")

    export Bernstein, BernsteinDerivative,
           Chebyshev, ChebyshevDerivative,
           ChebyshevT, ChebyshevTDerivative,
           ChebyshevU, ChebyshevUDerivative,
           Lagrange, LagrangeDerivative,
           LagrangeGauß, LagrangeLobatto,
           Legendre, LegendreDerivative

    include("bernstein.jl")
    include("chebyshev.jl")
    include("lagrange.jl")
    include("legendre.jl")

end
