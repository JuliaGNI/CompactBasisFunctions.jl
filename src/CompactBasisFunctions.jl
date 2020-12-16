module CompactBasisFunctions

    using ContinuumArrays
    using StaticArrays

    import Base: *
    import ContinuumArrays: Mul, QMul2, QuasiAdjoint, (..), @simplify

    export Derivative, ℵ₁

    include("vandermonde_matrix.jl")

    export Basis, basis, eachbasis, nbasis, grid, degree, order

    include("basis.jl")

    export Bernstein, BernsteinDerivative,
           Chebyshev, ChebyshevDerivative,
           ChebyshevT, ChebyshevTDerivative,
           ChebyshevU, ChebyshevUDerivative,
           Lagrange, LagrangeDerivative,
           LagrangeGauß, LagrangeLobatto

    include("bernstein.jl")
    include("chebyshev.jl")
    include("lagrange.jl")

end
