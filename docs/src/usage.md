```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# Usage

## Constructing a basis

Every basis is built from the number of basis functions, optionally preceded by an element
type, except [`Lagrange`](@ref), which is built from its nodes:

```jldoctest
julia> Bernstein(4), Legendre(4), ChebyshevT(4), ChebyshevU(4);

julia> Bernstein(Float32, 4), Legendre(BigFloat, 4), Chebyshev{1}(Float32, 4);

julia> Lagrange([0.0, 0.25, 1.0]), LagrangeGauß(4), LagrangeLobatto(4);
```

The element type defaults to `Float64`.

## Evaluation

A basis behaves like a matrix with one continuous axis, so evaluation is indexing: the first
index is the point, the second the basis function.

```jldoctest usage
julia> b = Legendre(3);

julia> b[0.5, 1]                     # one point, one basis function
0.0

julia> b(0.5, 1)                     # the same, as a function call
0.0

julia> b[0.5, :]                     # one point, all basis functions
3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:
  1.0
  0.0
 -1.118033988749895

julia> b[[0.25, 0.75], 1]            # several points, one basis function
2-element Vector{Float64}:
 -0.8660254037844386
  0.8660254037844386
```

Indexing with two collections gives the full matrix of values, points along the rows and
basis functions along the columns:

```jldoctest usage
julia> b[[0.0, 0.5, 1.0], :]
3×3 OffsetArray(::Matrix{Float64}, 1:3, 0:2) with eltype Float64 with indices 1:3×0:2:
 1.0  -1.73205   2.23607
 1.0   0.0      -1.11803
 1.0   1.73205   2.23607
```

### Indexing conventions

[`Bernstein`](@ref), [`Chebyshev`](@ref) and [`Legendre`](@ref) index their basis functions
from **0**; [`Lagrange`](@ref) indexes from **1**. Use `eachindex` and this never matters:

```jldoctest
julia> eachindex(Legendre(3))
OffsetArrays.IdOffsetRange(values=0:2, indices=0:2)

julia> eachindex(Lagrange([0.0, 0.5, 1.0]))
SOneTo(3)
```

### Evaluation outside the domain

`axes(b, 1)` is `Inclusion(0..1)` for every basis, but indexing does **not** check it.
Outside the interval the polynomial is extrapolated:

```jldoctest
julia> b = Bernstein(2);

julia> b[2.0, 0], b[2.0, 1]          # B₀(x) = 1-x and B₁(x) = x, evaluated at x = 2
(-1.0, 2.0)
```

This is deliberate — it is what makes the bases usable as local building blocks on an element
of a larger mesh — but note that guarantees which hold on ``[0,1]``, such as Bernstein
positivity and the convex-hull property, do not survive outside it.

## Accessors

```jldoctest
julia> b = ChebyshevT(4);

julia> nbasis(b), order(b), degree(b)
(4, 4, 3)

julia> nnodes(b)
4

julia> length(nodes(b)) == length(grid(b)) == nnodes(b)
true
```

| accessor | meaning | defined for |
|---|---|---|
| [`nbasis`](@ref) | number of basis functions | all |
| [`order`](@ref) | `== nbasis`; polynomials of degree `< order` are exact | all |
| [`degree`](@ref) | `== nbasis - 1`, the highest degree spanned | all |
| [`basis`](@ref) | the basis functions as callables | all |
| [`nodes`](@ref) | the nodes | nodal bases only |
| [`nnodes`](@ref) | number of nodes | nodal bases only |
| `grid` | same as `nodes`, under ContinuumArrays' name | nodal bases only |

The modal bases have no nodes, and say so:

```jldoctest
julia> nodes(Bernstein(4))
ERROR: Bernstein is a modal basis and has no nodes, so nodes is not defined for it.
[...]
```

See [Nodal and modal bases](@ref).

## Derivatives

The derivative operator comes from ContinuumArrays and is applied to the basis by
multiplication. The product is lazy: it stores the basis and evaluates on indexing.

```jldoctest usage2
julia> b = Legendre(3);

julia> d = Derivative(axes(b, 1));

julia> d * b isa LegendreDerivative
true

julia> b' isa LegendreDerivative        # adjoint is the same thing
true

julia> (d*b)[0.5, 2]
0.0

julia> (d*b)[1.0, 2]           # 6√5: the √(2j+1) of the basis carries to its derivative
13.416407864998739
```

Indexing a derivative supports the same four forms as the basis itself: `[x, j]`, `[x, :]`,
`[X, j]` and `[X, :]`.

Each basis has its own derivative type — [`BernsteinDerivative`](@ref),
[`ChebyshevDerivative`](@ref), [`LagrangeDerivative`](@ref), [`LegendreDerivative`](@ref) —
which is what makes the formula dispatch on the family.

## Expansions

A function in a basis is the product of the basis with a coefficient vector, which is again
lazy and materialises on evaluation:

```jldoctest
julia> b = LagrangeLobatto(3);

julia> f(x) = 1 + x^2;

julia> c = [f(x) for x in nodes(b)];            # nodal basis: coefficients are values

julia> sum(c[j] * b[0.3, j] for j in eachindex(b)) ≈ f(0.3)
true
```

For a modal basis the coefficients have to be computed. With [`Legendre`](@ref) that is one
quadrature sum per coefficient, since the basis is orthonormal:

```jldoctest
julia> b = Legendre(4);  quad = GaussLegendreQuadrature(10);

julia> f(x) = x^3 - x;

julia> c = [sum(weights(quad)[k] * f(nodes(quad)[k]) * b[nodes(quad)[k], j]
                for k in eachindex(nodes(quad))) for j in eachindex(b)];

julia> sum(c[j] * b[0.3, j] for j in eachindex(b)) ≈ f(0.3)   # c is indexed like b
true
```

## Element types and precision

The element type propagates through every internal quantity, so an arbitrary-precision basis
really carries its precision:

```jldoctest
julia> setprecision(BigFloat, 256) do
           b = Lagrange(gauss_legendre_nodes(BigFloat, 4))
           eltype(b.denom), eltype(b.diffs)
       end
(BigFloat, BigFloat)
```

For [`Chebyshev`](@ref) the element type must be able to *represent* the nodes, so an integer
type generally fails at construction; see [Element type](@ref).

## Working alongside QuadratureRules and ContinuumArrays

[`nodes`](@ref), [`nnodes`](@ref) and [`order`](@ref) are the same generic functions that
`QuadratureRules` extends — both packages import them from `GeometricBase` — so a basis and a
quadrature rule can be used together without qualifying anything:

```jldoctest
julia> using QuadratureRules

julia> order(Legendre(3)), order(GaussLegendreQuadrature(3))
(3, 6)

julia> nnodes(ChebyshevT(4)), nnodes(GaussLegendreQuadrature(3))
(4, 3)
```

`grid` is likewise ContinuumArrays' own function, not a second one of the same name.

!!! warning "`basis` is ambiguous with ContinuumArrays"
    ContinuumArrays exports its own `basis`, which means something different: for a basis
    object it returns *the basis itself*, where this package returns the collection of basis
    functions. Loading both with `using` therefore leaves the name ambiguous:

    ```julia
    using ContinuumArrays, CompactBasisFunctions
    basis(Legendre(3))   # UndefVarError: `basis` not defined
    ```

    Qualify it as `CompactBasisFunctions.basis`, or import the one that is wanted with
    `import CompactBasisFunctions: basis`.
