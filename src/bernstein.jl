
using OffsetArrays

"""
Bernstein polynomial ``B_{j,p}(x) = \\binom{p}{j} x^j (1-x)^{p-j}`` on the interval [0..1].

Evaluated in the closed form rather than through the recurrence
``B_{j,p} = (1-x) B_{j,p-1} + x B_{j-1,p-1}``: that recurrence has two indices, so
descending it recomputes shared subtrees and costs `O(2^p)` evaluations, which is 45 ms per
value at `p = 24`.

The binomial coefficient is accumulated as ``\\binom{p-j+k}{k}``, multiplying before
dividing so that every intermediate is an integer and the division is exact. Values are
exact while they stay below `2^53` for a floating-point `T`, which covers any degree at
which a Bernstein basis is numerically useful.

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


"""
Bernstein basis on the interval [0..1].
"""
struct Bernstein{T, BT} <: Basis{T}
    b::BT
    n::Int

    function Bernstein{T}(n::Integer) where {T}
        p = n-1
        b = OffsetArray([y -> _bernstein(i, p, y) for i in 0:p], 0:p)
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

Base.eltype(::Bernstein{T}) where {T} = T
Base.eachindex(B::Bernstein) = eachindex(B.b)
Base.axes(B::Bernstein) = (Inclusion(0..1), eachindex(B))

Base.hash(B::Bernstein, h::UInt) = hash(B.n, h)
Base.:(==)(B1::Bernstein, B2::Bernstein) = (B1.n == B2.n)
Base.isequal(B1::Bernstein{T1}, B2::Bernstein{T2}) where {T1,T2} = (T1 == T2 && B1 == B2)

Base.getindex(B::Bernstein, x::Number, j::Integer) = B(x,j)
Base.getindex(B::Bernstein, x::Number,  ::Colon) = [b(x) for b in B.b]
Base.getindex(B::Bernstein, X::AbstractVector, j::Integer) = B.(X,j)
Base.getindex(B::Bernstein, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in B.b]


## Derivative

function _eval_derivative(b::Bernstein{BT}, x::DT, i::Int) where {BT,DT}
    @assert i ≥ 0 && i < b.n
    (b.n-1) * ( _bernstein(i-1, b.n-2, x) - _bernstein(i, b.n-2, x) )
end

@simplify *(D::Derivative, B::Bernstein) = Mul(D,B)

const BernsteinDerivative = QMul2{<:Derivative,<:Bernstein}

Base.getindex(D::BernsteinDerivative, x::Number, j::Integer) = _eval_derivative(D.B, x, j)
Base.getindex(D::BernsteinDerivative, x::Number,  ::Colon) = [_eval_derivative(D.B, x, j) for j in eachindex(D.B)]
Base.getindex(D::BernsteinDerivative, X::AbstractVector, j::Integer) = [_eval_derivative(D.B, x, j) for x in X]
Base.getindex(D::BernsteinDerivative, X::AbstractVector,  ::Colon) = [_eval_derivative(D.B, x, j) for x in X, j in eachindex(D.B)]

Base.adjoint(B::Bernstein) = Derivative(axes(B,1)) * B
