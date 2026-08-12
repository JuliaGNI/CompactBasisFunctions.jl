
using OffsetArrays

"""
Legendre polynomial on the interval [-1..+1].

Evaluated by iterating Bonnet's recurrence ``j P_j = (2j-1) x P_{j-1} - (j-1) P_{j-2}``
upwards. The recursive formulation of the same recurrence descends into two subproblems per
step and recomputes shared subtrees, which costs `O(φʲ)` evaluations rather than `O(j)`.
"""
@inline function _legendre(j::Int, x::T) where {T}
    j < 0  && return zero(T)
    j == 0 && return one(T)

    local p₂ = one(T)
    local p₁ = x

    for k in 2:j
        p₂, p₁ = p₁, ( (2k-1) * p₁ * x - (k-1) * p₂ ) / k
    end

    return p₁
end

"""
Legendre basis on the interval [0..+1].
"""
struct Legendre{T, LT} <: Basis{T}
    b::LT
    n::Int

    function Legendre{T}(n::Integer) where {T}
        p = n-1
        b = OffsetArray([y -> _legendre(i, 2y-1) * sqrt(T(2i+1)) for i in 0:p], 0:p)
        new{T, typeof(b)}(b, n)
    end

end

Legendre(::Type{T}, n::Integer) where {T} = Legendre{T}(n)
Legendre(n::Integer) = Legendre(Float64, n)

(L::Legendre)(x::Number, j::Integer) = L.b[j](x)

basis(L::Legendre) = L.b
nbasis(L::Legendre) = L.n
order(L::Legendre) = nbasis(L)
degree(L::Legendre) = nbasis(L) - 1

nodes(L::Legendre) = _no_nodes(L, "nodes")
nnodes(L::Legendre) = _no_nodes(L, "nnodes")
ContinuumArrays.grid(L::Legendre) = _no_nodes(L, "grid")

Base.eltype(::Legendre{T}) where {T} = T
Base.eachindex(L::Legendre) = eachindex(L.b)
Base.axes(L::Legendre) = (Inclusion(0..1), eachindex(L))

Base.hash(L::Legendre, h::UInt) = hash(L.n, h)
Base.:(==)(L1::Legendre, L2::Legendre) = (L1.n == L2.n)
Base.isequal(L1::Legendre{T1}, L2::Legendre{T2}) where {T1,T2} = (T1 == T2 && L1 == L2)
Base.isapprox(L1::Legendre, L2::Legendre; kwargs...) = (L1.n == L2.n)

Base.getindex(L::Legendre, x::Number, j::Integer) = L(x,j)
Base.getindex(L::Legendre, x::Number,  ::Colon) = [b(x) for b in L.b]
Base.getindex(L::Legendre, X::AbstractVector, j::Integer) = L.(X,j)
Base.getindex(L::Legendre, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in L.b]


## Derivative

"""
Derivative of Legendre polynomial on the interval [-1..+1].

Obtained by differentiating Bonnet's recurrence term by term,

```math
j P_j' = (2j-1) P_{j-1} + (2j-1) x P_{j-1}' - (j-1) P_{j-2}' ,
```

and iterated upwards alongside ``P_j`` itself, which the right-hand side needs.
"""
@inline function _legendre_derivative(j::Int, x::T) where {T}
    j <= 0 && return zero(T)
    j == 1 && return one(T)

    local p₂ = one(T)    # P_{k-2}
    local p₁ = x         # P_{k-1}
    local d₂ = zero(T)   # P_{k-2}'
    local d₁ = one(T)    # P_{k-1}'

    for k in 2:j
        # evaluated before p₁ advances, since the derivative needs P_{k-1}
        local d = ( (2k-1) * p₁ + (2k-1) * d₁ * x - (k-1) * d₂ ) / k
        p₂, p₁ = p₁, ( (2k-1) * p₁ * x - (k-1) * p₂ ) / k
        d₂, d₁ = d₁, d
    end

    return d₁
end

@simplify *(D::Derivative, L::Legendre) = Mul(D,L)

const LegendreDerivative = QMul2{<:Derivative,<:Legendre}

"""
Evaluate derivative of Legendre polynomial on the interval [0..+1].
"""
function _eval(D::LegendreDerivative, x::DT, j::Int) where {DT}
    @boundscheck j ≥ 0 && j < nbasis(D.B) || throw(BoundsError(D.B, j))
    _legendre_derivative(j, promote_type(eltype(D.B), DT)(2x-1)) * 2
end

Base.getindex(D::LegendreDerivative, x::Number, j::Integer) = _eval(D, x, j)
Base.getindex(D::LegendreDerivative, x::Number,  ::Colon) = [_eval(D, x, j) for j in eachindex(D.B)]
Base.getindex(D::LegendreDerivative, X::AbstractVector, j::Integer) = [_eval(D, x, j) for x in X]
Base.getindex(D::LegendreDerivative, X::AbstractVector,  ::Colon) = [_eval(D, x, j) for x in X, j in eachindex(D.B)]

Base.adjoint(L::Legendre) = Derivative(axes(L,1)) * L
