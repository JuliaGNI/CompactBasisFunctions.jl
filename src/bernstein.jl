
using OffsetArrays

@inline function _bernstein(j::Int, p::Int, x::T) where {T}
    if j < 0 || j > p
        return zero(T)
    else
        if p == 0
            return one(T)
        else
            return _bernstein(j, p-1, x) * (1-x) + _bernstein(j-1, p-1, x) * x
        end
    end
end


"""
Bernstein basis on the interval [0..1].
"""
struct Bernstein{T, BT} <: Basis{T}
    b::BT
    n::Integer

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
eachbasis(B::Bernstein) = eachindex(B.b)
order(B::Bernstein) = nbasis(B)
degree(B::Bernstein) = nbasis(B) - 1

Base.axes(B::Bernstein) = (Inclusion(0..1), eachbasis(B))

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
Base.getindex(D::BernsteinDerivative, x::Number,  ::Colon) = [_eval_derivative(D.B, x, j) for j in eachbasis(D.B)]
Base.getindex(D::BernsteinDerivative, X::AbstractVector, j::Integer) = [_eval_derivative(D.B, x, j) for x in X]
Base.getindex(D::BernsteinDerivative, X::AbstractVector,  ::Colon) = [_eval_derivative(D.B, x, j) for x in X, j in eachbasis(D.B)]

Base.adjoint(B::Bernstein) = Derivative(axes(B,1)) * B
