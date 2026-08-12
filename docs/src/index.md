```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# CompactBasisFunctions

This package provides a set of basis functions, mostly compactly supported, which are
implemented as [ContinuumArrays](https://github.com/JuliaApproximation/ContinuumArrays.jl).
Bases are accessed like arrays with continuous dimensions, e.g. as `b[0.1, 2]` to evaluate the
second basis function in the point `0.1`. Operations such as derivatives, inner products and
mass matrices are implemented as [LazyArray](https://github.com/JuliaArrays/LazyArrays.jl)
operations, providing a high-level linear algebra interface. Functions in a basis are
represented by a lazy multiplication of the basis and a vector of coefficients, which
materialises only upon evaluation of the product.

## Installation

*CompactBasisFunctions.jl* and all of its dependencies can be installed via the Julia REPL by
typing

```julia
]add CompactBasisFunctions
```

## The four bases

All four span the polynomials of a given degree on the reference interval ``[0,1]``; they
differ in what their coefficients mean and in what properties they guarantee.

| basis | coefficients are | carries nodes | notable for |
|---|---|---|---|
| [`Lagrange`](@ref) | values at the nodes | yes | interpolation needs no solve |
| [`Chebyshev`](@ref) | expansion coefficients | yes | node sets that defeat the Runge phenomenon |
| [`Legendre`](@ref) | expansion coefficients | no | orthonormal, so the mass matrix is `I` |
| [`Bernstein`](@ref) | control values | no | non-negative, partition of unity, convex hull |

[Polynomial Approximation](@ref) explains the distinctions — interpolation versus projection,
nodal versus modal, the reference interval, the choice of nodes — and each basis then has its
own page.

## Basic usage

Construct a basis and index it:

```jldoctest intro
julia> using CompactBasisFunctions

julia> b = Legendre(3);

julia> nbasis(b), degree(b)
(3, 2)

julia> b[0.5, 0], b[0.5, 1], b[0.5, 2]
(1.0, 0.0, -1.118033988749895)
```

Note that these three bases index their functions from `0`; [`Lagrange`](@ref) indexes from
`1`. Writing `eachindex(b)` avoids having to remember which:

```jldoctest intro
julia> b[0.5, :]
3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:
  1.0
  0.0
 -1.118033988749895
```

Differentiate by multiplying with a `Derivative` over the basis's own axis:

```jldoctest intro
julia> d = Derivative(axes(b, 1));

julia> (d*b)[0.5, 1], (d*b)[1.0, 2]
(3.4641016151377544, 13.416407864998739)
```

A nodal basis makes interpolation immediate, since its coefficients are the function values:

```jldoctest
julia> b = LagrangeLobatto(3);

julia> f(x) = 1 + x^2;

julia> c = [f(x) for x in nodes(b)];

julia> sum(c[j] * b[0.3, j] for j in eachindex(b)) ≈ f(0.3)
true
```

See [Usage](@ref) for the full interface, including the array-valued indexing forms,
expansions, element types and how the accessors relate to those of
[QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl).

## References

If you use CompactBasisFunctions.jl in your work, please consider citing it by

```
@misc{Kraus:2020:CompactBasisFunctions,
  title={CompactBasisFunctions.jl: Compactly supported basis functions in Julia},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/CompactBasisFunctions.jl}},
  doi={10.5281/zenodo.4317806}
}
```
