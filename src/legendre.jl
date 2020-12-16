
using OffsetArrays

"""
Legendre polynomial on the interval [-1..+1].
"""
@inline function _legendre(j::Int, x::T) where {T}
    if j <= 0
        return one(T)
    elseif j == 1
        return x
    else
        return ( (2j-1) * _legendre(j-1, x) * x - (j-1) * _legendre(j-2, x) ) / j
    end
end

"""
Legendre basis on the interval [0..+1].
"""
struct Legendre{T, LT} <: Basis{T}
    b::LT
    n::Int

    function Legendre{T}(n::Integer) where {T}
        p = n-1
        b = OffsetArray([y -> _legendre(i, 2y-1) for i in 0:p], 0:p)
        new{T, typeof(b)}(b, n)
    end

end

Legendre(::Type{T}, n::Integer) where {T} = Legendre{T}(n)
Legendre(n::Integer) = Legendre(Float64, n)

function _eval(l::Legendre{LT}, x::DT, j::Int) where {LT,DT}
    @assert j ≥ 0 && j < l.n
    _legendre(j, 2x-1)
end

(L::Legendre)(x::Number, j::Integer) = L.b[j](x)

basis(L::Legendre) = L.b
nbasis(L::Legendre) = L.n
order(L::Legendre) = nbasis(L)
degree(L::Legendre) = nbasis(L) - 1

Base.eltype(::Legendre{T}) where {T} = T
Base.eachindex(L::Legendre) = eachindex(L.b)
Base.axes(L::Legendre) = (Inclusion(0..1), eachindex(L))

Base.hash(L::Legendre, h::UInt) = hash(L.n, h)
Base.:(==)(L1::Legendre, L2::Legendre) = (L1.n == L2.n)
Base.isequal(L1::Legendre{T1}, L2::Legendre{T2}) where {T1,T2} = (T1 == T2 && L1 == L2)

Base.getindex(L::Legendre, x::Number, j::Integer) = L(x,j)
Base.getindex(L::Legendre, x::Number,  ::Colon) = [b(x) for b in L.b]
Base.getindex(L::Legendre, X::AbstractVector, j::Integer) = L.(X,j)
Base.getindex(L::Legendre, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in L.b]


## Derivative

"""
Derivative of Legendre polynomial on the interval [-1..+1].
"""
@inline function _legendre_derivative(j::Int, x::T) where {T}
    if j <= 0
        return zero(T)
    elseif j == 1
        return one(T)
    else
        return ( (2j-1) * _legendre(j-1, x) + (2j-1) * _legendre_derivative(j-1, x) * x - (j-1) * _legendre_derivative(j-2, x) ) / j
    end
end

@simplify *(D::Derivative, L::Legendre) = Mul(D,L)

const LegendreDerivative = QMul2{<:Derivative,<:Legendre}

"""
Evaluate derivative of Legendre polynomial on the interval [0..+1].
"""
function _eval(D::LegendreDerivative, x::DT, j::Int) where {DT}
    @assert j ≥ 0 && j < nbasis(D.B)
    _legendre_derivative(j, promote_type(eltype(D.B), DT)(2x-1)) * 2
end

Base.getindex(D::LegendreDerivative, x::Number, j::Integer) = _eval(D, x, j)
Base.getindex(D::LegendreDerivative, x::Number,  ::Colon) = [_eval(D, x, j) for j in eachindex(D.B)]
Base.getindex(D::LegendreDerivative, X::AbstractVector, j::Integer) = [_eval(D, x, j) for x in X]
Base.getindex(D::LegendreDerivative, X::AbstractVector,  ::Colon) = [_eval(D, x, j) for x in X, j in eachindex(D.B)]

Base.adjoint(L::Legendre) = Derivative(axes(L,1)) * L
