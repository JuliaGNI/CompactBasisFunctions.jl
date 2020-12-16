
import FastTransforms: chebyshevpoints#, transform, itransform
import OffsetArrays: OffsetArray


@inline function _chebyshev(::Val{1}, j::Int, x::T) where {T}
    if j < 0
        return zero(T)
    elseif j == 0
        return one(T)
    elseif j == 1
        return x
    else
        return _chebyshev(Val(1), j-1, x) * 2x - _chebyshev(Val(1), j-2, x)
    end
end

@inline function _chebyshev(::Val{2}, j::Int, x::T) where {T}
    if j < 0
        return zero(T)
    elseif j == 0
        return one(T)
    elseif j == 1
        return 2x
    else
        return _chebyshev(Val(2), j-1, x) * 2x - _chebyshev(Val(2), j-2, x)
    end
end


"""
Chebyshev basis on the interval [-1..+1].
"""
struct Chebyshev{kind, T, BT, XT <: AbstractVector{T}} <: Basis{T}
    b::BT
    x::XT

    function Chebyshev{kind, T}(n::Int) where {kind, T}
        p = n-1
        x = chebyshevpoints(T, n, Val(kind))
        b = OffsetArray([y -> _chebyshev(Val(kind), i, y) for i in 0:p], 0:p)
        new{kind, T, typeof(b), typeof(x)}(b, x)
    end

    Chebyshev{kind}(::Type{T}, n::Int) where {kind, T} = Chebyshev{kind, T}(n)
    Chebyshev{kind}(n::Int) where {kind} = Chebyshev{kind, Float64}(n)
end

Chebyshev(::Type{T}, n::Int, ::Val{kind}) where {kind, T} = Chebyshev{kind, T}(n)
Chebyshev(n::Int, ::Val{kind}) where {kind} = Chebyshev(Float64, n, Val(kind))

const ChebyshevT = Chebyshev{1}
const ChebyshevU = Chebyshev{2}

(C::Chebyshev)(x::Number, j::Int) = C.b[j](x)

basis(C::Chebyshev) = C.b
nodes(C::Chebyshev) = C.x
nbasis(C::Chebyshev) = length(basis(C))
nnodes(C::Chebyshev) = length(nodes(C))
eachbasis(C::Chebyshev) = eachindex(basis(C))
eachnode(C::Chebyshev)  = eachindex(nodes(C))
order(C::Chebyshev)  = nnodes(C)
degree(C::Chebyshev) = nnodes(C) - 1

Base.eltype(::Chebyshev{kind,T}) where {kind,T} = T
Base.axes(C::Chebyshev) = (Inclusion(0..1), eachbasis(C))
ContinuumArrays.grid(C::Chebyshev) = nodes(C)

Base.hash(C::Chebyshev{kind}, h::UInt) where {kind} = hash(C.x, hash(kind, h))
Base.:(==)(C1::Chebyshev{kind1}, C2::Chebyshev{kind2}) where {kind1,kind2} = (C1.x == C2.x && kind1 == kind2)
Base.isequal(C1::Chebyshev{kind1,T1}, C2::Chebyshev{kind2,T2}) where {kind1,T1,kind2,T2} = (T1 == T2 && C1 == C2)
Base.isapprox(C1::Chebyshev, C2::Chebyshev; kwargs...) = isapprox(C1.x, C2.x; kwargs...)

Base.getindex(C::Chebyshev, x::Number, j::Int) = C(x,j)
Base.getindex(C::Chebyshev, x::Number,  ::Colon) = [b(x) for b in C.b]
Base.getindex(C::Chebyshev, X::AbstractVector, j::Int) = C.(X,j)
Base.getindex(C::Chebyshev, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in C.b]


## Derivative

@simplify *(D::Derivative, C::Chebyshev) = Mul(D,C)

const ChebyshevDerivative  = QMul2{<:Derivative,<:Chebyshev}
const ChebyshevTDerivative = QMul2{<:Derivative,<:ChebyshevT}
const ChebyshevUDerivative = QMul2{<:Derivative,<:ChebyshevU}

function _eval(D::ChebyshevTDerivative, x::DT, i::Int) where {DT}
    local C = D.B
    local x̃ = promote_type(eltype(C), DT)(x)
    @assert i ≥ 0 && i < nbasis(C)
    _chebyshev(Val(2), i-1, x̃) * i
end

function _eval(D::ChebyshevUDerivative, x::DT, i::Int) where {DT}
    local C = D.B
    local x̃ = promote_type(eltype(C), DT)(x)
    @assert i ≥ 0 && i < nbasis(C)
    ( _chebyshev(Val(1), i+1, x̃) * (i+1) - _chebyshev(Val(2), i, x̃) * x̃ ) / (x̃^2 - 1)
end

Base.getindex(D::ChebyshevDerivative, x::Number, j::Integer) = _eval(D, x, j)
Base.getindex(D::ChebyshevDerivative, x::Number,  ::Colon) = [_eval(D, x, j) for j in eachbasis(D.B)]
Base.getindex(D::ChebyshevDerivative, X::AbstractVector, j::Integer) = [_eval(D, x, j) for x in X]
Base.getindex(D::ChebyshevDerivative, X::AbstractVector,  ::Colon) = [_eval(D, x, j) for x in X, j in eachbasis(D.B)]

Base.adjoint(C::Chebyshev) = Derivative(axes(C,1)) * C
