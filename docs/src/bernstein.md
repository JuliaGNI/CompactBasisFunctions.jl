```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# Bernstein basis

## Definition

The Bernstein basis of degree ``p`` on ``[0,1]`` is

```math
B_{j,p}(x) = \binom{p}{j} \, x^j \, (1-x)^{p-j} , \qquad j = 0, \dots, p .
```

For its history and its role in computer-aided geometric design, see
[Farouki:2012:BernsteinPolynomialBasis](@cite).

## Domain

No shift is needed: the basis is defined on ``[0,1]`` to begin with, which is also this
package's reference interval.

```jldoctest
julia> b = Bernstein(3);

julia> axes(b, 1)
Inclusion(0 .. 1)

julia> nbasis(b), order(b), degree(b)
(3, 3, 2)
```

## No nodes

Bernstein is a **modal** basis. Its coefficients are not values of the function at points —
they are the *control values* of a Bézier curve, which the curve approaches but in general
does not pass through. So it has no nodes:

```jldoctest
julia> grid(Bernstein(4))
ERROR: Bernstein is a modal basis and has no nodes, so grid is not defined for it.
[...]
```

See [Nodal and modal bases](@ref).

## Partition of unity, positivity, convex hull

On ``[0,1]`` the basis functions are non-negative and sum to one:

```jldoctest
julia> b = Bernstein(5);

julia> all(b[x, j] ≥ 0 for x in 0:0.1:1, j in eachindex(b))
true

julia> sum(b[0.42, j] for j in eachindex(b)) ≈ 1
true
```

Together these give the **convex-hull property**: an expansion ``\sum_j c_j B_{j,p}(x)`` is a
convex combination of its coefficients, so its values lie between ``\min_j c_j`` and
``\max_j c_j`` for every ``x`` in the interval. No other basis here has this. It is what makes
Bernstein the basis of choice when a bound on the represented function has to be guaranteed
rather than merely observed, and it is the reason Bézier curves stay inside their control
polygon.

```jldoctest
julia> b = Bernstein(5);  c = [0.0, 1.0, -2.0, 0.5, 3.0];

julia> vals = [sum(c[j+1] * b[x, j] for j in eachindex(b)) for x in 0:0.01:1];

julia> minimum(c) ≤ minimum(vals) && maximum(vals) ≤ maximum(c)
true
```

## Endpoint interpolation

The first and last basis functions interpolate the endpoints, and all the others vanish at
both:

```jldoctest
julia> b = Bernstein(4);

julia> [b[0.0, j] for j in eachindex(b)]
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
 1.0
 0.0
 0.0
 0.0

julia> [b[1.0, j] for j in eachindex(b)]
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
 0.0
 0.0
 0.0
 1.0
```

So an expansion begins at its first coefficient and ends at its last — convenient when
boundary values have to be imposed, even though the basis is modal.

## Derivative

The derivative of a Bernstein polynomial is a Bernstein expansion of one degree less:

```math
B_{j,p}'(x) = p \, \big( B_{j-1,p-1}(x) - B_{j,p-1}(x) \big) ,
```

with the convention that ``B_{j,q}`` vanishes outside ``0 \le j \le q``, which is what makes
the ``j = 0`` and ``j = p`` cases come out right — and what makes the degenerate `Bernstein(1)`
case, where ``p-1 = -1``, give zero as it should.

```jldoctest
julia> b = Bernstein(3);  d = Derivative(axes(b, 1));

julia> [(d*b)[0.5, j] for j in eachindex(b)]
3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:
 -1.0
  0.0
  1.0

julia> sum((d*b)[0.42, j] for j in eachindex(b)) |> abs < 1e-14      # rows sum to zero
true
```

## Evaluation in closed form

The two-index recurrence

```math
B_{j,p} = (1-x) \, B_{j,p-1} + x \, B_{j-1,p-1}
```

cannot be collapsed into a single upward sweep for one ``j``, so this basis is evaluated from
the closed form instead, with the binomial coefficient accumulated as ``\binom{p-j+k}{k}``,
multiplying before dividing. The partial products are integers, so in exact arithmetic every
division comes out even; in `Float64` the accumulation reproduces `binomial(p, j)` exactly for
``p \le 54``, and from ``p = 55`` the intermediate product outgrows the integers `Float64`
represents exactly and the coefficient picks up rounding.

Note that `/` promotes, so an integer element type does not survive the accumulation: the
values are floating point whenever ``j \ge 1``.

!!! warning "Changed in 0.3.0"
    Earlier versions descended the recurrence, which recomputes shared subtrees and cost
    ``O(2^p)`` per value — 45 ms for a single degree-24 value. Values now differ from those in
    the last bit for some arguments, by at most ``5 \times 10^{-16}`` relative, and are on
    balance closer to exact.

## Pitfalls

- **No nodes**: `nodes`, `nnodes` and `grid` all throw.
- **The coefficients are not function values.** An expansion generally does not pass through
  its coefficients — that is what makes the convex-hull property possible. If interpolation is
  wanted, use [`Lagrange`](@ref).
- **Indexing starts at 0**, as for Chebyshev and Legendre, but not Lagrange.
- Positivity and the convex-hull property hold **on ``[0,1]`` only**. Outside it the basis
  extrapolates and the functions do take negative values; see
  [Evaluation outside the domain](@ref).
- The binomial coefficient is exact in `Float64` up to ``p = 54``, which covers any degree at
  which this basis is numerically useful; see [Evaluation in closed form](@ref) above.
