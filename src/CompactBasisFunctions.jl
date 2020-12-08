module CompactBasisFunctions

    using ContinuumArrays
    using StaticArrays

    import Base: *
    import ContinuumArrays: Mul, QMul2, QuasiAdjoint, (..), @simplify

    export Derivative, ℵ₁

    include("vandermonde_matrix.jl")

    export degree, eachnode, grid, nbasis, nnodes, nodes, order

    include("basis.jl")

    export Lagrange, LagrangeDerivative

    include("lagrange.jl")

end
