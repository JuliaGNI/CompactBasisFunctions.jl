```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# Legendre basis

## Definition

The Legendre polynomials ``P_j`` on ``[-1,+1]`` satisfy Bonnet's recurrence

```math
j \, P_j(\tilde{x}) = (2j-1) \, \tilde{x} \, P_{j-1}(\tilde{x}) - (j-1) \, P_{j-2}(\tilde{x}) ,
\qquad P_0 = 1 , \quad P_1 = \tilde{x} ,
```

which is how they are evaluated here, iterated upwards; see
[Evaluation by recurrence](@ref). The basis provided is the shifted and **normalised**
version,

```math
L_j(x) = \sqrt{2j+1} \, P_j(2x-1) , \qquad j = 0, \dots, p .
```

## Domain and the shift

The argument is shifted by ``\tilde{x} = 2x-1`` onto ``[0,1]``; see
[The reference interval](@ref).

```jldoctest
julia> l = Legendre(3);

julia> axes(l, 1)
Inclusion(0 .. 1)

julia> l[0.5, 1]                            # x = 0.5 is x̃ = 0, and P₁(0) = 0
0.0

julia> l[1.0, 1] ≈ sqrt(3)                  # L₁(1) = √3 P₁(1)
true
```

## No nodes

Legendre is a **modal** basis: its coefficients are the amplitudes of an expansion, not
values at points. It has no nodes, and the accessors say so rather than failing obscurely:

```jldoctest
julia> nbasis(Legendre(4)), order(Legendre(4)), degree(Legendre(4))
(4, 4, 3)

julia> nnodes(Legendre(4))
ERROR: Legendre is a modal basis and has no nodes, so nnodes is not defined for it.
[...]
```

See [Nodal and modal bases](@ref). To recover the coefficients of a given function, project
it — which is what the orthonormality below makes cheap.

## Orthonormality

This is the reason to choose this basis. Since

```math
\int_0^1 P_i(2x-1) \, P_j(2x-1) \, dx = \frac{\delta_{ij}}{2j+1} ,
```

the factor ``\sqrt{2j+1}`` in the definition makes the basis orthonormal on ``[0,1]``:

```math
\int_0^1 L_i(x) \, L_j(x) \, dx = \delta_{ij} .
```

The mass matrix is therefore the identity, and the coefficients of an ``L^2`` projection are
simply the inner products against the basis functions — no linear system.

```jldoctest
julia> l = Legendre(5);  quad = GaussLegendreQuadrature(10);

julia> M = [sum(weights(quad)[k] * l[nodes(quad)[k], i] * l[nodes(quad)[k], j]
                for k in eachindex(nodes(quad))) for i in 0:4, j in 0:4];

julia> maximum(abs, M - [i == j for i in 0:4, j in 0:4]) < 1e-14
true
```

Projecting a function is then a matter of one quadrature sum per coefficient:

```jldoctest
julia> l = Legendre(4);  quad = GaussLegendreQuadrature(10);

julia> f(x) = x^3 - x;

julia> c = [sum(weights(quad)[k] * f(nodes(quad)[k]) * l[nodes(quad)[k], j]
                for k in eachindex(nodes(quad))) for j in 0:3];

julia> sum(c[j+1] * l[0.3, j] for j in 0:3) ≈ f(0.3)   # exact: f has degree 3
true
```

Unlike [Chebyshev](@ref "Chebyshev basis"), whose orthogonality needs a weight, this holds in
the plain unweighted inner product.

## Derivative

Differentiating Bonnet's recurrence term by term gives

```math
j \, P_j' = (2j-1) \, P_{j-1} + (2j-1) \, \tilde{x} \, P_{j-1}' - (j-1) \, P_{j-2}' ,
```

iterated upwards alongside ``P_j`` itself, which the right-hand side needs. The chain-rule
factor `2` of the shift is applied at the end, and the ``\sqrt{2j+1}`` normalisation carries
through.

```jldoctest
julia> l = Legendre(3);  d = Derivative(axes(l, 1));

julia> [(d*l)[1.0, j] for j in eachindex(l)]     # 0, 2√3, 6√5
3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:
  0.0
  3.4641016151377544
 13.416407864998739
```

!!! warning "Fixed in 0.3.0"
    Earlier versions omitted the ``\sqrt{2j+1}`` factor here, although the basis functions
    carry it, so the derivative was wrong by exactly that factor for every ``j > 0``.

## Pitfalls

- **No nodes**: `nodes`, `nnodes` and `grid` all throw. Use [`Lagrange`](@ref) or
  [`Chebyshev`](@ref) if a nodal basis is wanted.
- **The normalisation is part of the basis.** These are not the plain ``P_j``; every basis
  function is scaled by ``\sqrt{2j+1}``, so `l[1.0, 1]` is ``\sqrt{3}`` rather than ``1``.
- **Indexing starts at 0**, as for Bernstein and Chebyshev, but not Lagrange.
- Evaluation outside ``[0,1]`` extrapolates rather than erroring; see
  [Evaluation outside the domain](@ref).
