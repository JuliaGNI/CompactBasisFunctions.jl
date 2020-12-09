module CompactBasisFunctions

    using ContinuumArrays
    using StaticArrays

    import Base: *
    import ContinuumArrays: Mul, QMul2, QuasiAdjoint, (..), @simplify

    export Derivative, ℵ₁

    include("vandermonde_matrix.jl")

    export basis, eachbasis, nbasis,
           nodes, eachnode, nnodes,
           grid, degree, order

    include("basis.jl")

    export Lagrange, LagrangeDerivative,
           LagrangeGauß, LagrangeLobatto

    include("lagrange.jl")

end
