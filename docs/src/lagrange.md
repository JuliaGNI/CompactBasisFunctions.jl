```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# Lagrange basis

## Definition

Given ``n`` distinct nodes ``x_1, \dots, x_n``, the Lagrange basis is

```math
\ell_j(x) = \prod_{i \neq j} \frac{x - x_i}{x_j - x_i} , \qquad j = 1, \dots, n ,
```

a basis of the polynomials of degree ``\le n-1``. This is the *cardinal* basis of the nodes,

```math
\ell_j(x_i) = \delta_{ij} ,
```

which is the property everything else follows from.

## Domain and nodes

Nodal, and the only basis here whose nodes the user supplies. The declared domain is
``[0,1]``, but the nodes are whatever they are given as: they need not be sorted, and nothing
forces them into the interval.

Note that [`Lagrange`](@ref) indexes its basis functions from **1**, where the other three
bases index from **0**. Write loops over `eachindex(b)` rather than over `0:degree(b)` and
this never comes up.

```jldoctest
julia> l = Lagrange([0.0, 0.5, 1.0]);

julia> eachindex(l)
SOneTo(3)

julia> nnodes(l), nbasis(l), order(l), degree(l)
(3, 3, 3, 2)
```

## Interpolation is free

Because the basis is cardinal, the coefficients of an interpolant *are* the values at the
nodes — there is no system to solve. This is the reason to choose a Lagrange basis.

```jldoctest
julia> l = LagrangeLobatto(3);

julia> f(x) = 1 + x^2;

julia> c = [f(x) for x in nodes(l)];        # the coefficients are just the values

julia> sum(c[j] * l[0.3, j] for j in eachindex(l)) ≈ f(0.3)
true
```

The last line is exact up to rounding because ``f`` is a polynomial of degree 2 and the basis
spans degree 2.

## Cardinality and partition of unity

```jldoctest
julia> l = Lagrange([0.0, 0.1, 0.7, 1.0]);

julia> all(l[nodes(l)[i], j] ≈ (i == j) for i in 1:4, j in 1:4)
true

julia> sum(l[0.42, j] for j in eachindex(l)) ≈ 1
true
```

The cardinal functions do take negative values between the nodes, so unlike
[Bernstein](@ref "Bernstein basis") this partition of unity is not a convex combination.

## Which nodes to use

Equidistant nodes are a poor choice at higher degree — see [Choice of nodes](@ref) — so two
better sets are provided directly:

```jldoctest
julia> collect(nodes(LagrangeGauß(3)))      # strictly inside the interval
3-element Vector{Float64}:
 0.11270166537925831
 0.5
 0.8872983346207417

julia> collect(nodes(LagrangeLobatto(3)))   # includes both endpoints
3-element Vector{Float64}:
 0.0
 0.5
 1.0
```

The nodes are stored in an `SVector`, hence the `collect` above; `nodes` itself returns the
static vector.

Use [`LagrangeGauß`](@ref) when the expansion is going to be integrated, and
[`LagrangeLobatto`](@ref) when the boundary values have to be imposed or read off.

## Derivative

Differentiating the product gives

```math
\ell_j'(x) = \sum_{l \neq j} \frac{1}{x_j - x_l} \prod_{i \neq j, l} \frac{x - x_i}{x_j - x_i} ,
```

which is what [`LagrangeDerivative`](@ref) evaluates, reusing the table of node differences
the basis caches at construction.

```jldoctest
julia> l = LagrangeLobatto(3);

julia> d = Derivative(axes(l, 1));

julia> sum((d*l)[0.3, j] for j in eachindex(l)) |> abs < 1e-14      # rows sum to zero
true
```

## Element type and precision

All internal quantities — the node differences and the reciprocal denominators — are formed
in the element type of the nodes, so a `BigFloat` basis carries its full precision. Before
version 0.3.0 those buffers were allocated in `Float64` regardless, so a `BigFloat` basis
reported `BigFloat` while carrying only double precision, in the basis and in its derivative
alike.

```jldoctest
julia> setprecision(BigFloat, 256) do
           l = Lagrange(gauss_legendre_nodes(BigFloat, 4))
           d = Derivative(axes(l, 1))
           z = BigFloat(1) / 7
           (abs(sum(l[z, j] for j in eachindex(l)) - 1) < 1e-70,
            abs(sum((d*l)[z, j] for j in eachindex(l))) < 1e-70)
       end
(true, true)
```

## Pitfalls

- **The nodes must be distinct.** The denominators are products of node differences, so a
  repeated node is a division by zero. This is rejected rather than silently producing `Inf`:

  ```jldoctest
  julia> Lagrange([0.0, 0.5, 0.5, 1.0])
  ERROR: ArgumentError: the nodes of a Lagrange basis must be distinct, got [0.0, 0.5, 0.5, 1.0]
  [...]
  ```

  What is tested is the differences themselves, not the node list under `isequal` — the
  comparison `allunique` uses, and a different question from the one the denominators ask.
  `0.0` and `-0.0` are distinct under `isequal`, yet their difference is zero:

  ```jldoctest
  julia> allunique([0.0, -0.0])
  true

  julia> Lagrange([0.0, -0.0])
  ERROR: ArgumentError: the nodes of a Lagrange basis must be distinct, got [0.0, -0.0]
  [...]
  ```

- **The nodes must be finite.** A `NaN` or `Inf` node is `isequal` to nothing else and so
  looks perfectly distinct, while poisoning every difference it takes part in:

  ```jldoctest
  julia> Lagrange([0.0, NaN, 1.0])
  ERROR: ArgumentError: the nodes of a Lagrange basis must be finite, got [0.0, NaN, 1.0]
  [...]
  ```

- **The product of the differences must be representable.** Nodes that are distinct and
  finite can still multiply out to zero or to an infinity — packed into a narrow range the
  product underflows, spread over a wide one it overflows — and either way the denominator
  is unusable. The nodes are not at fault there, so the message says so and points at the
  remedy:

  ```jldoctest
  julia> Lagrange([0.0, 1e160, 2e160, 3e160])
  ERROR: ArgumentError: the nodes of a Lagrange basis are distinct and finite, but the product of the differences from node 1 is -Inf in Float64, so the denominator it gives is not usable; rescale the nodes or widen the element type
  [...]
  ```

  Scaling such a node set to `[0,1]` and mapping the argument along with it is the way out;
  a wider element type buys room too, but only a fixed amount of it.

- **Indexing starts at 1**, unlike the other three bases.
- **Equidistant nodes at high degree** will diverge; see [Choice of nodes](@ref).
- The nodes are stored in an `SVector`, so the number of nodes is part of the type. Building
  a basis with very many nodes is therefore a compile-time cost, not just a run-time one.

For the barycentric reformulation, which is the preferred way to evaluate a Lagrange
interpolant at many points and is more stable than the product form, see
[Berrut:2004:BarycentricLagrange](@cite) and [Higham:2004:BarycentricStability](@cite).
