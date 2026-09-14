
@doc raw"""
    PolynomialBasis{T} <: Basis{T}

Supertype of the four bases of this package: a ContinuumArrays `Basis` of element type `T`,
spanning the polynomials of degree ``\le p`` on the reference interval ``[0,1]``.

A subtype supplies one internal method, `_eval`, evaluating one basis function at one point,
and the data that method needs. Everything a caller sees follows from it here —
[`basis`](@ref), [`nbasis`](@ref), [`order`](@ref),
[`degree`](@ref), `eachindex`, `axes`, the four indexing forms, equality and the derivative
product — so that a family's own file holds only what is peculiar to it.

The hierarchy splits by whether a basis carries nodes, into [`NodalBasis`](@ref) and
[`ModalBasis`](@ref); see [Nodal and modal bases](@ref).
"""
abstract type PolynomialBasis{T} <: Basis{T} end

@doc raw"""
    NodalBasis{T} <: PolynomialBasis{T}

A [`PolynomialBasis`](@ref) built from a set of nodes, which it stores in a field `x`:
[`Lagrange`](@ref) and [`Chebyshev`](@ref). It has one basis function per node.

[`nodes`](@ref), [`nnodes`](@ref) and `grid` are answered from that field. See
[Nodal and modal bases](@ref).
"""
abstract type NodalBasis{T} <: PolynomialBasis{T} end

@doc raw"""
    ModalBasis{T} <: PolynomialBasis{T}

A [`PolynomialBasis`](@ref) without nodes, built from a number of basis functions, which it
stores in a field `n`: [`Bernstein`](@ref) and [`Legendre`](@ref).

Its coefficients are those of an expansion rather than values at points, so [`nodes`](@ref),
[`nnodes`](@ref) and `grid` have no answer and throw. See [Nodal and modal bases](@ref).
"""
abstract type ModalBasis{T} <: PolynomialBasis{T} end

"""
Error for the accessors that a modal basis cannot answer: its coefficients are those of an
expansion, not values at points, so it has no nodes and no grid. Saying so beats the
`MethodError` about ContinuumArrays' `grid_axis` that `grid` would otherwise produce.
"""
function _no_nodes(b::ModalBasis, f)
    error("$(nameof(typeof(b))) is a modal basis and has no nodes, so $(f) is not defined for it.")
end

"""
The arithmetic an evaluation is carried out in: the wider of the basis's element type `T` and
the type `S` of the point it is evaluated at.

A basis has to evaluate its recurrence in `T` even at an argument of lower precision, or it
returns a value whose type claims a precision the value does not carry. The promotion never
narrows: an argument wider than `T` keeps its own precision.

Everything an evaluation forms belongs in this type, including a constant factor such as the
``\\sqrt{2j+1}`` of a Legendre basis function: a factor formed in `T` alone holds the whole
result to the precision of `T`, whatever the argument carries.

The conversion has to come before any arithmetic on the argument, in particular before the
shift `2x-1` onto ``[-1,1]`` that the Chebyshev and Legendre recurrences want. Converting the
shifted value instead leaves the shift itself running in the argument's type, and widening its
rounded result only records that rounding in more digits. `2x-1` is exact in `Float64` for
`x ≥ 0.25`, so such an omission hides from any test that samples only the upper part of the
domain.

For an abstract `T` such as the `Integer` of `ChebyshevU(Integer, 2)` this is an abstract type,
and the conversion is then a no-op that leaves the argument as it is.
"""
@inline _evaltype(::Type{T}, ::Type{S}) where {T, S} = promote_type(T, S)

"""
Evaluate basis function `j` at the point `x`, in the wider of the basis's element type and the
argument's, cf. [`_evaltype`](@ref).

One method per family for a basis, and one per family for its derivative — that formula is
what distinguishes the families, and it is all they have to supply. The index is in range,
having been checked by the `getindex` methods below, which are shared.
"""
function _eval end

@doc raw"""
    basis(b::PolynomialBasis)

Return the basis functions of `b` as callables, indexed as `b` itself is, so that
`basis(b)[j](x) == b[x,j]`.

The functions are built on each call. Prefer indexing `b` directly; this accessor exists for
the rare case that the individual functions are needed as values.

Note that ContinuumArrays exports a different `basis`, which for a basis object returns the
basis itself. Loading both packages with `using` therefore makes the name ambiguous; qualify
it, or import the one that is wanted.
"""
basis(b::PolynomialBasis) = [x -> b[x, j] for j in eachindex(b)]

@doc raw"""
    nbasis(b::PolynomialBasis)

Return the number of basis functions of `b`.

This equals [`nnodes`](@ref) for the nodal bases, where each basis function belongs to one
node, and is defined for the modal bases too, where `nnodes` is not.
"""
nbasis(b::ModalBasis) = b.n
nbasis(b::NodalBasis) = nnodes(b)

@doc raw"""
    order(b::PolynomialBasis)

Return the order of `b`, the number of basis functions, hence the number of coefficients an
expansion in `b` has.

`order(b) == nbasis(b)` and `order(b) == degree(b) + 1` for every basis here. The name is
the one the ecosystem uses for the accuracy of an approximation: a basis of order ``p``
spans the polynomials of degree ``\le p-1`` and so reproduces them exactly.
"""
order(b::PolynomialBasis) = nbasis(b)

@doc raw"""
    degree(b::PolynomialBasis)

Return the degree of `b`, the highest polynomial degree it spans.

`degree(b) == order(b) - 1 == nbasis(b) - 1` for every basis here.
"""
degree(b::PolynomialBasis) = nbasis(b) - 1

@doc raw"""
    nodes(b::NodalBasis)

Return the nodes of `b`, the points ``x_i \in [0,1]`` that its construction is based on.

Defined for the nodal bases [`Lagrange`](@ref) and [`Chebyshev`](@ref). The modal bases
[`Bernstein`](@ref) and [`Legendre`](@ref) have no nodes and throw an error; see
[Nodal and modal bases](@ref).

See also [`nnodes`](@ref) for their number and `grid` for the same vector under the name
ContinuumArrays uses.
"""
nodes(b::NodalBasis) = b.x

@doc raw"""
    nnodes(b::NodalBasis)

Return the number of nodes of `b`, i.e. the length of [`nodes`](@ref).

Defined only for the nodal bases; see [`nodes`](@ref).
"""
nnodes(b::NodalBasis) = length(nodes(b))

ContinuumArrays.grid(b::NodalBasis) = nodes(b)

nodes(b::ModalBasis) = _no_nodes(b, "nodes")
nnodes(b::ModalBasis) = _no_nodes(b, "nnodes")
ContinuumArrays.grid(b::ModalBasis) = _no_nodes(b, "grid")

"""
What a basis is made of: the family it belongs to, and the data that family is built from —
the number of functions for a [`ModalBasis`](@ref), the nodes for a [`NodalBasis`](@ref).
Two bases are equal when their keys are.

The family is part of the key so that bases of different families never compare equal, and
the element type is not, so that `Bernstein(Float64, 3) == Bernstein(Integer, 3)`; `isequal`
is what tells those two apart.
"""
function _key end

Base.hash(b::PolynomialBasis, h::UInt) = hash(_key(b), h)
Base.:(==)(b1::PolynomialBasis, b2::PolynomialBasis) = _key(b1) == _key(b2)
function Base.isequal(b1::PolynomialBasis, b2::PolynomialBasis)
    eltype(b1) == eltype(b2) && b1 == b2
end

function Base.isapprox(b1::PolynomialBasis, b2::PolynomialBasis; kwargs...)
    family1, data1 = _key(b1)
    family2, data2 = _key(b2)
    family1 == family2 && _isapprox(data1, data2; kwargs...)
end

# a tolerance applies to the nodes of a nodal basis, and to nothing a modal basis is made of:
# two Bernstein bases of different degree are not approximately the same basis, however loose
# the tolerance, so the count is compared exactly and the keywords are dropped
_isapprox(x1, x2; kwargs...) = isapprox(x1, x2; kwargs...)
_isapprox(n1::Integer, n2::Integer; kwargs...) = n1 == n2

# the basis functions are numbered from 0, as the degrees they carry are. `Lagrange` is the
# exception, numbering its functions from 1 along with the nodes they belong to.
Base.eachindex(b::PolynomialBasis) = IdOffsetRange(Base.OneTo(nbasis(b)), -1)
Base.axes(b::PolynomialBasis) = (Inclusion(0..1), eachindex(b))

(b::PolynomialBasis)(x::Number, j::Integer) = b[x, j]

function Base.getindex(b::PolynomialBasis, x::Number, j::Integer)
    @boundscheck j ∈ eachindex(b) || throw(BoundsError(b, j))
    _eval(b, x, j)
end

function Base.getindex(b::PolynomialBasis, x::Number, ::Colon)
    [_eval(b, x, j) for j in eachindex(b)]
end

function Base.getindex(b::PolynomialBasis, X::AbstractVector, j::Integer)
    @boundscheck j ∈ eachindex(b) || throw(BoundsError(b, j))
    [_eval(b, x, j) for x in X]
end

function Base.getindex(b::PolynomialBasis, X::AbstractVector, ::Colon)
    [_eval(b, x, j) for x in X, j in eachindex(b)]
end

## Derivative

@simplify *(D::Derivative, b::PolynomialBasis) = Mul(D, b)

"""
    PolynomialBasisDerivative

The type of `Derivative(axes(b,1)) * b` for any [`PolynomialBasis`](@ref) `b`, equivalently
of `b'`, and the supertype of the four per-family derivative types.

A lazy product: it stores the basis and evaluates the derivative on indexing, in the same
four forms as the basis itself. See [Derivatives](@ref).
"""
const PolynomialBasisDerivative = QMul2{<:Derivative, <:PolynomialBasis}

Base.adjoint(b::PolynomialBasis) = Derivative(axes(b, 1)) * b

function Base.getindex(D::PolynomialBasisDerivative, x::Number, j::Integer)
    @boundscheck j ∈ eachindex(D.B) || throw(BoundsError(D.B, j))
    _eval(D, x, j)
end

function Base.getindex(D::PolynomialBasisDerivative, x::Number, ::Colon)
    [_eval(D, x, j) for j in eachindex(D.B)]
end

function Base.getindex(D::PolynomialBasisDerivative, X::AbstractVector, j::Integer)
    @boundscheck j ∈ eachindex(D.B) || throw(BoundsError(D.B, j))
    [_eval(D, x, j) for x in X]
end

function Base.getindex(D::PolynomialBasisDerivative, X::AbstractVector, ::Colon)
    [_eval(D, x, j) for x in X, j in eachindex(D.B)]
end
