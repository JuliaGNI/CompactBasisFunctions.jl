
_not_implemented(f, b) = error("$(f) is not implemented for $(typeof(b)).")

"""
Error for the accessors that a modal basis cannot answer: its coefficients are those of an
expansion, not values at points, so it has no nodes and no grid. Saying so beats the
`MethodError` about ContinuumArrays' `grid_axis` that `grid` would otherwise produce.
"""
_no_nodes(b, f) =
    error("$(nameof(typeof(b))) is a modal basis and has no nodes, so $(f) is not defined for it.")

"""
The arithmetic an evaluation is carried out in: the wider of the basis's element type `T` and
the type `S` of the point it is evaluated at.

A basis has to evaluate its recurrence in `T` even at an argument of lower precision, or it
returns a value whose type claims a precision the value does not carry — a `Legendre{BigFloat}`
evaluated at a `Float64` point used to run Bonnet's recurrence entirely in `Float64` and widen
only the trailing normalisation factor. The promotion never narrows: an argument wider than `T`
keeps its own precision.

The conversion has to come before any arithmetic on the argument, in particular before the
shift `2x-1` onto ``[-1,1]`` that the Chebyshev and Legendre recurrences want. Converting the
shifted value instead leaves the shift itself running in the argument's type, and widening its
rounded result only records that rounding in more digits — `2x-1` is exact in `Float64` for
`x ≥ 0.25`, so the omission hides from any test that samples only the upper part of the domain.

For an abstract `T` such as the `Integer` of `ChebyshevU(Integer, 2)` this is an abstract type,
and the conversion is then a no-op that leaves the argument as it is.
"""
@inline _evaltype(::Type{T}, ::Type{S}) where {T, S} = promote_type(T, S)

# `nodes` and `nnodes` are deliberately left unimplemented for the modal bases, which have
# no nodes at all; those carry their own methods with a message saying so. A generic
# `grid(::Basis)` is not defined here on purpose: ContinuumArrays dispatches `grid` on
# `Any`, so a method on its `Basis` supertype would also capture its own spline bases.

@doc raw"""
    basis(b::Basis)

Return the collection of basis functions of `b`, indexed as `b` itself is.

Each element is a callable evaluating one basis function, so that
`basis(b)[j](x) == b[x,j]`. Prefer indexing `b` directly; this accessor exists for the
rare case that the individual functions are needed as values.

Note that ContinuumArrays exports a different `basis`, which for a basis object returns the
basis itself. Loading both packages with `using` therefore makes the name ambiguous; qualify
it, or import the one that is wanted.
"""
basis(b::Basis) = _not_implemented("basis", b)

@doc raw"""
    nodes(b::Basis)

Return the nodes of `b`, the points ``x_i \in [0,1]`` that its construction is based on.

Defined for the nodal bases [`Lagrange`](@ref) and [`Chebyshev`](@ref). The modal bases
[`Bernstein`](@ref) and [`Legendre`](@ref) have no nodes and throw an error; see
[Nodal and modal bases](@ref).

See also [`nnodes`](@ref) for their number and `grid` for the same vector under the name
ContinuumArrays uses.
"""
nodes(b::Basis) = _not_implemented("nodes", b)

@doc raw"""
    nbasis(b::Basis)

Return the number of basis functions of `b`, i.e. the length of [`basis`](@ref).

This equals [`nnodes`](@ref) for the nodal bases, where each basis function belongs to one
node, and is defined for the modal bases too, where `nnodes` is not.
"""
nbasis(b::Basis) = _not_implemented("nbasis", b)

@doc raw"""
    nnodes(b::Basis)

Return the number of nodes of `b`, i.e. the length of [`nodes`](@ref).

Defined only for the nodal bases; see [`nodes`](@ref).
"""
nnodes(b::Basis) = _not_implemented("nnodes", b)

@doc raw"""
    order(b::Basis)

Return the order of `b`, the number of basis functions, hence the number of coefficients an
expansion in `b` has.

`order(b) == nbasis(b)` and `order(b) == degree(b) + 1` for every basis here. The name is
the one the ecosystem uses for the accuracy of an approximation: a basis of order ``p``
spans the polynomials of degree ``\le p-1`` and so reproduces them exactly.
"""
order(b::Basis) = _not_implemented("order", b)

@doc raw"""
    degree(b::Basis)

Return the degree of `b`, the highest polynomial degree it spans.

`degree(b) == order(b) - 1 == nbasis(b) - 1` for every basis here.
"""
degree(b::Basis) = _not_implemented("degree", b)
