
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

Returns zero outside `0 ≤ j ≤ p`, which is what makes the calls from the derivative, one
degree lower, work for a basis of a single function.
"""
@inline function _bernstein(j::Integer, p::Integer, x::T) where {T}
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
struct Bernstein{T} <: ModalBasis{T}
    n::Int
end

Bernstein(::Type{T}, n::Integer) where {T} = Bernstein{T}(n)
Bernstein(n::Integer) = Bernstein(Float64, n)

_key(B::Bernstein) = (Bernstein, B.n)

function _eval(B::Bernstein, x, j::Integer)
    local T = _evaltype(eltype(B), typeof(x))
    _bernstein(j, degree(B), T(x))
end

## Derivative

"""
    BernsteinDerivative

The type of `Derivative(axes(b,1)) * b` for a [`Bernstein`](@ref) basis `b`, equivalently of
`b'`. See [`PolynomialBasisDerivative`](@ref) and [Derivatives](@ref).
"""
const BernsteinDerivative = QMul2{<:Derivative, <:Bernstein}

function _eval(D::BernsteinDerivative, x, j::Integer)
    local p = degree(D.B)
    local T = _evaltype(eltype(D.B), typeof(x))
    local x̃ = T(x)
    p * (_bernstein(j-1, p-1, x̃) - _bernstein(j, p-1, x̃))
end
