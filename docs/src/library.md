```@meta
CurrentModule = CompactBasisFunctions
```

# Library

The complete API of *CompactBasisFunctions.jl*. See [Usage](@ref) for how these fit together
and [Polynomial Approximation](@ref) for the underlying theory.

```@index
```

## Generic API

Every basis is a subtype of ContinuumArrays' `Basis`, and answers the accessors below. The
methods shown here are the fallbacks on the abstract type, which report that a concrete basis
has not implemented them; the per-family pages document what each one returns.

`basis`, `degree`, `nodes`, `nnodes` and `order` are imported from `GeometricBase`, and `grid`
from `ContinuumArrays`, so that the packages of the ecosystem extend one generic function per
accessor rather than defining one each. `nbasis` is this package's own.

```@docs
basis
nbasis
nodes
nnodes
order
degree
```

## Basis functions

The four families, in the order the manual introduces them.

### Lagrange

```@docs
Lagrange
LagrangeGauß
LagrangeLobatto
```

### Chebyshev

```@docs
Chebyshev
```

### Legendre

```@docs
Legendre
```

### Bernstein

```@docs
Bernstein
```

## Derivative types

Applying `Derivative(axes(b,1))` to a basis produces a lazy product, one type per family, that
evaluates the derivative on indexing. See [Derivatives](@ref) for how they are used.

```@docs
BernsteinDerivative
ChebyshevDerivative
ChebyshevTDerivative
ChebyshevUDerivative
LagrangeDerivative
LegendreDerivative
```

## Vandermonde matrices

Not exported; reach them as `CompactBasisFunctions.vandermonde_matrix` and
`CompactBasisFunctions.vandermonde_matrix_inverse`.

```@docs
vandermonde_matrix
vandermonde_matrix_inverse
```
