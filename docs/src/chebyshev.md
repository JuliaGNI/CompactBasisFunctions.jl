```@meta
CurrentModule = CompactBasisFunctions
DocTestSetup = quote
    using CompactBasisFunctions
    using QuadratureRules
end
```

# Chebyshev basis

## Definition

The Chebyshev polynomials of the first and second kind are

```math
T_j(\tilde{x}) = \cos(j\theta) , \qquad
U_j(\tilde{x}) = \frac{\sin\big((j+1)\theta\big)}{\sin\theta} , \qquad
\tilde{x} = \cos\theta ,
```

on ``\tilde{x} \in [-1,+1]``. Both satisfy the same three-term recurrence,

```math
\phi_j = 2\tilde{x} \, \phi_{j-1} - \phi_{j-2} ,
```

differing only in the starting values: ``T_0 = 1, T_1 = \tilde{x}`` and
``U_0 = 1, U_1 = 2\tilde{x}``. This is how they are evaluated here — iterated upwards, which
is the stable direction; see [Evaluation by recurrence](@ref).

For the classical theory see [Mason:2003:ChebyshevPolynomials](@cite), and for their role in
spectral methods [Boyd:2001:ChebyshevFourier](@cite).

## Domain and the shift

The basis functions are the polynomials evaluated at ``\tilde{x} = 2x-1``, so the basis lives
on ``[0,1]`` like everything else in this package; see [The reference interval](@ref).

```jldoctest
julia> t = ChebyshevT(3);

julia> axes(t, 1)
Inclusion(0 .. 1)

julia> t[0.0, 2], t[0.5, 2], t[1.0, 2]      # T₂(x̃) at x̃ = -1, 0, +1
(1.0, -1.0, 1.0)
```

!!! warning "Changed in 0.3.0"
    Before version 0.3.0 the nodes and basis functions lived on ``[-1,+1]`` while `axes`
    already declared `Inclusion(0..1)`, so `grid(C)` did not lie inside `axes(C,1)`. Every
    node, basis function and derivative value therefore changed in 0.3.0, and the nodes are
    now returned in ascending rather than descending order.

## Nodes

Unlike the modal bases, a Chebyshev basis carries nodes: the Chebyshev points of the matching
kind, obtained from `QuadratureRules`, shifted onto ``[0,1]`` and returned **ascending**.

Those of the **first kind** are the roots of ``T_n``,
``\tilde{x}_i = \cos\big((2i-1)\pi/2n\big)``, and lie strictly inside the interval. Those of
the **second kind** are the extrema of ``T_{n-1}``,
``\tilde{x}_i = \cos\big((i-1)\pi/(n-1)\big)``, and include both endpoints.

```jldoctest
julia> nodes(ChebyshevT(3))                 # strictly interior
3-element Vector{Float64}:
 0.06698729810778067
 0.5
 0.9330127018922193

julia> nodes(ChebyshevU(3))                 # includes 0 and 1
3-element Vector{Float64}:
 0.0
 0.5
 1.0

julia> issorted(nodes(ChebyshevT(5)))
true
```

Because the second-kind points are the extrema of ``T_{n-1}``, they need at least two points
to exist, so `ChebyshevU(1)` is an error.

These are the node sets that make high-degree interpolation well behaved: their Lebesgue
constant grows logarithmically in the degree rather than exponentially, which is what defeats
the [Runge phenomenon](@ref "Choice of nodes").

## Derivative

The first kind differentiates into the second,

```math
T_j'(\tilde{x}) = j \, U_{j-1}(\tilde{x}) ,
```

which is division-free and used directly. The second kind has the closed form

```math
U_j'(\tilde{x}) = \frac{(j+1) T_{j+1}(\tilde{x}) - \tilde{x} \, U_j(\tilde{x})}{\tilde{x}^2 - 1} ,
```

which is ``0/0`` at ``\tilde{x} = \pm 1`` — that is, at both endpoints of ``[0,1]``. It is
therefore *not* used; the derivative is obtained from the differentiated recurrence

```math
U_j' = 2 U_{j-1} + 2\tilde{x} \, U_{j-1}' - U_{j-2}' ,
```

carried alongside ``U_j`` itself, which has no singularity anywhere. Both carry the
chain-rule factor `2` from the shift.

```jldoctest
julia> u = ChebyshevU(4);  d = Derivative(axes(u, 1));

julia> [(d*u)[1.0, i] for i in eachindex(u)]    # finite at the endpoint
4-element OffsetArray(::Vector{Float64}, 0:3) with eltype Float64 with indices 0:3:
  0.0
  4.0
 16.0
 40.0

julia> [(d*u)[1.0, i] for i in eachindex(u)] == [2i*(i+1)*(i+2)//3 for i in eachindex(u)]
true
```

!!! warning "Fixed in 0.3.0"
    Earlier versions evaluated the closed form, so `(d*ChebyshevU(n))[0.0, j]` and
    `[1.0, j]` returned `NaN` for every `j`.

## Orthogonality

The Chebyshev polynomials are orthogonal with respect to a *weighted* inner product, with
weight ``(1-\tilde{x}^2)^{-1/2}`` for the first kind and ``(1-\tilde{x}^2)^{+1/2}`` for the
second. They are **not** orthogonal in the unweighted ``L^2`` inner product, so the mass
matrix of this basis is not diagonal — unlike [Legendre](@ref "Legendre basis").

## Element type

`T` is the element type of the nodes and must be able to represent them, which in general
means a floating-point type. The conversion happens at construction, so a type that cannot
represent the nodes fails there rather than on first use:

```jldoctest
julia> nodes(ChebyshevU(Integer, 2))        # exactly 0 and 1, so this works
2-element Vector{Integer}:
 0
 1

julia> try ChebyshevT(Integer, 3) catch e; typeof(e) end
InexactError
```

The node values are computed in `BigFloat` by `QuadratureRules` and rounded once to `T`, so
they are correctly rounded whatever `T` is.

`T` also governs the arithmetic of an evaluation, which runs in the wider of `T` and the type
of the point; see [Element types and precision](@ref).

## Pitfalls

- **`ChebyshevU(1)` is an error**: second-kind points need ``n \ge 2``.
- **Values changed in 0.3.0** with the move onto ``[0,1]``, including the order of the nodes.
- **Indexing starts at 0**, as for Bernstein and Legendre, but not Lagrange.
- Evaluation outside ``[0,1]`` is not rejected; it extrapolates the polynomial. See
  [Evaluation outside the domain](@ref).
