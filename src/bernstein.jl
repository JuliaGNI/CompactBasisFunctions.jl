
using OffsetArrays

"""
Bernstein polynomial ``B_{j,p}(x) = \\binom{p}{j} x^j (1-x)^{p-j}`` on the interval [0..1].

Evaluated in the closed form rather than through the recurrence
``B_{j,p} = (1-x) B_{j,p-1} + x B_{j-1,p-1}``: that recurrence has two indices, so
descending it recomputes shared subtrees and costs `O(2^p)` evaluations, which is 45 ms per
value at `p = 24`.

The binomial coefficient is accumulated as ``\\binom{p-j+k}{k}``, multiplying before
dividing: the partial products are integers, so in exact arithmetic every division comes
out even. In `Float64` this reproduces `binomial(p, j)` exactly for `p ≤ 54`; from `p = 55`
the intermediate product outgrows the integers `Float64` represents exactly and the
coefficient picks up rounding. That covers any degree at which a Bernstein basis is
numerically useful.

Note that `/` promotes, so the accumulator leaves an integer `T`: the return type is
floating point whenever `j ≥ 1`.

Returns zero outside `0 ≤ j ≤ p`, which is what makes the `p = n-2` calls from the
derivative work for `n = 1`.
"""
@inline function _bernstein(j::Int, p::Int, x::T) where {T}
    (j < 0 || j > p) && return zero(T)

    local c = one(T)

    for k in 1:j
        c = c * (p - j + k) / k
    end

    return c * x^j * (1-x)^(p-j)
end

@doc raw"""
    Bernstein(n)
    Bernstein(T, n)

The Bernstein basis of `n` functions, i.e. of degree ``p = n-1``, on the interval ``[0,1]``,

```math
B_{j,p}(x) = \binom{p}{j} \, x^j \, (1-x)^{p-j} , \qquad j = 0, \dots, p ,
```

indexed from `0`. `T` is the element type and defaults to `Float64`.

The basis is *modal*: its coefficients are not values at points, so it has no
[`nodes`](@ref) and no `grid`, and those accessors throw. The coefficients are the control
values of a Bézier curve.

On ``[0,1]`` the basis functions are non-negative and form a partition of unity,
``\sum_j B_{j,p}(x) = 1``, so an expansion is a convex combination of its coefficients and
therefore lies in their convex hull. The first and last function interpolate the endpoints,
``B_{0,p}(0) = B_{p,p}(1) = 1``, while every other function vanishes at both.

```jldoctest
julia> b = Bernstein(3);

julia> b[0.0, 0], b[0.5, 0], b[1.0, 0]     # B₀ interpolates the left endpoint
(1.0, 0.25, 0.0)

julia> b[0.0, 2], b[0.5, 2], b[1.0, 2]     # B₂ interpolates the right endpoint
(0.0, 0.25, 1.0)

julia> sum(b[0.3, j] for j in eachindex(b)) ≈ 1
true
```

The derivative is again a Bernstein expansion, of one degree less,

```math
B_{j,p}'(x) = p \, \big( B_{j-1,p-1}(x) - B_{j,p-1}(x) \big) ,
```

and is obtained as `Derivative(axes(b,1)) * b`; see [Derivatives](@ref).

See also [`Legendre`](@ref) for the other modal basis, and [Bernstein basis](@ref) for the
full discussion.
"""
struct Bernstein{T, BT} <: Basis{T}
    b::BT
    n::Int

    function Bernstein{T}(n::Integer) where {T}
        p = n-1
        # evaluated in the wider of T and the argument type, cf. `_evaltype`
        b = OffsetArray([y -> _bernstein(i, p, _evaltype(T, typeof(y))(y)) for i in 0:p], 0:p)
        new{T, typeof(b)}(b, n)
    end
end

Bernstein(::Type{T}, n::Integer) where {T} = Bernstein{T}(n)
Bernstein(n::Integer) = Bernstein(Float64, n)

(B::Bernstein)(x::Number, j::Integer) = B.b[j](x)

basis(B::Bernstein) = B.b
nbasis(B::Bernstein) = B.n
order(B::Bernstein) = nbasis(B)
degree(B::Bernstein) = nbasis(B) - 1

nodes(B::Bernstein) = _no_nodes(B, "nodes")
nnodes(B::Bernstein) = _no_nodes(B, "nnodes")
ContinuumArrays.grid(B::Bernstein) = _no_nodes(B, "grid")

Base.eltype(::Bernstein{T}) where {T} = T
Base.eachindex(B::Bernstein) = eachindex(B.b)
Base.axes(B::Bernstein) = (Inclusion(0..1), eachindex(B))

Base.hash(B::Bernstein, h::UInt) = hash(B.n, h)
Base.:(==)(B1::Bernstein, B2::Bernstein) = (B1.n == B2.n)
Base.isequal(B1::Bernstein{T1}, B2::Bernstein{T2}) where {T1, T2} = (T1 == T2 && B1 == B2)
Base.isapprox(B1::Bernstein, B2::Bernstein; kwargs...) = (B1.n == B2.n)

Base.getindex(B::Bernstein, x::Number, j::Integer) = B(x, j)
Base.getindex(B::Bernstein, x::Number, ::Colon) = [b(x) for b in B.b]
Base.getindex(B::Bernstein, X::AbstractVector, j::Integer) = B.(X, j)
Base.getindex(B::Bernstein, X::AbstractVector, ::Colon) = [b(x) for x in X, b in B.b]

## Derivative

function _eval_derivative(b::Bernstein{T}, x::DT, i::Int) where {T, DT}
    @boundscheck i ≥ 0 && i < b.n || throw(BoundsError(b, i))
    local x̃ = _evaltype(T, DT)(x)
    (b.n-1) * (_bernstein(i-1, b.n-2, x̃) - _bernstein(i, b.n-2, x̃))
end

@simplify *(D::Derivative, B::Bernstein) = Mul(D, B)

"""
    BernsteinDerivative

The type of `Derivative(axes(b,1)) * b` for a [`Bernstein`](@ref) basis `b`, equivalently of
`b'`.

A lazy product: it stores the basis and evaluates the derivative on indexing, so
`(d*b)[x,j]` is ``B_{j,p}'(x)``. See [Derivatives](@ref).
"""
const BernsteinDerivative = QMul2{<:Derivative, <:Bernstein}

Base.getindex(D::BernsteinDerivative, x::Number, j::Integer) = _eval_derivative(D.B, x, j)
function Base.getindex(D::BernsteinDerivative, x::Number, ::Colon)
    [_eval_derivative(D.B, x, j) for j in eachindex(D.B)]
end
function Base.getindex(D::BernsteinDerivative, X::AbstractVector, j::Integer)
    [_eval_derivative(D.B, x, j) for x in X]
end
function Base.getindex(D::BernsteinDerivative, X::AbstractVector, ::Colon)
    [_eval_derivative(D.B, x, j) for x in X, j in eachindex(D.B)]
end

Base.adjoint(B::Bernstein) = Derivative(axes(B, 1)) * B
