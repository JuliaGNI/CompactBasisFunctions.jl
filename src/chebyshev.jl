
import QuadratureRules: chebyshev_nodes
import OffsetArrays: OffsetArray


"""
Chebyshev polynomial of the first kind on the interval [-1..+1].

Evaluated by iterating the three-term recurrence ``T_j = 2x T_{j-1} - T_{j-2}`` upwards.
The recursive formulation of the same recurrence descends into two subproblems per step and
recomputes shared subtrees, which costs `O(φʲ)` evaluations rather than `O(j)`.
"""
@inline function _chebyshev(::Val{1}, j::Int, x::T) where {T}
    j < 0  && return zero(T)
    j == 0 && return one(T)

    local t₂ = one(T)
    local t₁ = x

    for _ in 2:j
        t₂, t₁ = t₁, t₁ * 2x - t₂
    end

    return t₁
end

"""
Chebyshev polynomial of the second kind on the interval [-1..+1].

Evaluated by iterating ``U_j = 2x U_{j-1} - U_{j-2}``, cf. [`_chebyshev`](@ref).
"""
@inline function _chebyshev(::Val{2}, j::Int, x::T) where {T}
    j < 0  && return zero(T)
    j == 0 && return one(T)

    local u₂ = one(T)
    local u₁ = 2x

    for _ in 2:j
        u₂, u₁ = u₁, u₁ * 2x - u₂
    end

    return u₁
end


@doc raw"""
    Chebyshev{kind}(n)
    Chebyshev{kind}(T, n)
    ChebyshevT(n)     # kind = 1
    ChebyshevU(n)     # kind = 2

Chebyshev basis of the first (`kind = 1`) or second (`kind = 2`) kind on the interval
``[0,1]``, of `n` functions, i.e. of degree ``p = n-1``, indexed from `0`.

The basis functions are the Chebyshev polynomials evaluated at `2x-1`, and the nodes are the
Chebyshev points shifted onto [0..1] and returned in ascending order.

```math
T_j(\tilde{x}) = \cos(j \arccos \tilde{x}) , \qquad
U_j(\tilde{x}) = \frac{\sin\big((j+1) \arccos \tilde{x}\big)}{\sin(\arccos \tilde{x})} ,
\qquad \tilde{x} = 2x-1 ,
```

both satisfying the recurrence ``\phi_j = 2\tilde{x} \phi_{j-1} - \phi_{j-2}``, which is how
they are evaluated; see [The reference interval](@ref) for the shift and its chain-rule
factor `2` on derivatives.

Unlike the modal bases, a Chebyshev basis carries [`nodes`](@ref): the Chebyshev points of
the corresponding kind, from `QuadratureRules`, shifted onto ``[0,1]`` and **ascending**.
Those of the first kind are the roots of ``T_n`` and lie strictly inside the interval; those
of the second kind are the extrema of ``T_{n-1}`` and include both endpoints, which is why
`kind = 2` requires `n ≥ 2`.

```jldoctest
julia> t = ChebyshevT(3);

julia> t[0.0, 2], t[0.5, 2], t[1.0, 2]      # T₂(x̃) at x̃ = -1, 0, +1
(1.0, -1.0, 1.0)

julia> nodes(ChebyshevU(3))                 # includes both endpoints
3-element Vector{Float64}:
 0.0
 0.5
 1.0

julia> all(x -> x ∈ axes(t, 1), nodes(t))   # the nodes lie in the domain
true
```

`T` is the element type of the nodes, and must be able to represent them. This means a
floating-point type in general; an integer-like `T` only works in the special cases where the
nodes are exactly representable, such as `ChebyshevU(Integer, 2)`, whose nodes are 0 and 1.
Otherwise the constructor throws an `InexactError`.

See also [Chebyshev basis](@ref) for the full discussion, and [`Lagrange`](@ref) for the
other nodal basis.
"""
struct Chebyshev{kind, T, BT, XT <: AbstractVector{T}} <: Basis{T}
    b::BT
    x::XT

    function Chebyshev{kind, T}(n::Int) where {kind, T}
        p = n-1
        # chebyshev_nodes returns the points on [0,1]; shift_nodes widens integer
        # element types, so convert back to keep XT <: AbstractVector{T}
        x = convert(Vector{T}, chebyshev_nodes(T, n, Val(kind)))
        # evaluated in the wider of T and the argument type, as the derivatives below already are
        b = OffsetArray([y -> _chebyshev(Val(kind), i, _evaltype(T, typeof(y))(2y-1)) for i in 0:p], 0:p)
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
order(C::Chebyshev)  = nnodes(C)
degree(C::Chebyshev) = nnodes(C) - 1

Base.eltype(::Chebyshev{kind,T}) where {kind,T} = T
Base.eachindex(C::Chebyshev) = eachindex(C.b)
Base.axes(C::Chebyshev) = (Inclusion(0..1), eachindex(C))
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

"""
    ChebyshevDerivative

The type of `Derivative(axes(b,1)) * b` for a [`Chebyshev`](@ref) basis `b` of either kind,
equivalently of `b'`.

A lazy product: it stores the basis and evaluates the derivative on indexing. The two kinds
need different formulas, so the evaluation dispatches on the narrower
[`ChebyshevTDerivative`](@ref) and [`ChebyshevUDerivative`](@ref). See
[Derivatives](@ref).
"""
const ChebyshevDerivative  = QMul2{<:Derivative,<:Chebyshev}

"""
    ChebyshevTDerivative

[`ChebyshevDerivative`](@ref) narrowed to the first kind, evaluated through
``T_j' = j \\, U_{j-1}``.
"""
const ChebyshevTDerivative = QMul2{<:Derivative,<:ChebyshevT}

"""
    ChebyshevUDerivative

[`ChebyshevDerivative`](@ref) narrowed to the second kind, evaluated through the
differentiated recurrence rather than the closed form, which is singular at the endpoints.
"""
const ChebyshevUDerivative = QMul2{<:Derivative,<:ChebyshevU}

"""
Derivative of the Chebyshev polynomial of the first kind on the interval [-1..+1],
via ``T_j' = j \\, U_{j-1}``.
"""
@inline function _chebyshev_derivative(::Val{1}, j::Int, x::T) where {T}
    j ≤ 0 && return zero(T)
    return _chebyshev(Val(2), j-1, x) * j
end

"""
Derivative of the Chebyshev polynomial of the second kind on the interval [-1..+1].

Obtained by differentiating the recurrence ``U_j = 2x U_{j-1} - U_{j-2}`` term by term,

```math
U_j' = 2 U_{j-1} + 2x U_{j-1}' - U_{j-2}' ,
```

carried alongside ``U_j`` itself. The closed form
``U_j' = ((j+1) T_{j+1} - x U_j) / (x^2 - 1)`` is not used because it is ``0/0`` at
``x = \\pm 1``, i.e. at both endpoints of the interval, where the derivative is perfectly
finite.
"""
@inline function _chebyshev_derivative(::Val{2}, j::Int, x::T) where {T}
    j ≤ 0 && return zero(T)

    local u₂ = one(T)      # U_0
    local u₁ = 2x          # U_1
    local d₂ = zero(T)     # U_0'
    local d₁ = T(2)        # U_1'

    for _ in 2:j
        u₂, u₁ = u₁, u₁ * 2x - u₂
        d₂, d₁ = d₁, 2u₂ + d₁ * 2x - d₂
    end

    return d₁
end


"""
Evaluate derivative of Chebyshev polynomial of the first kind on the interval [0..1].
"""
function _eval(D::ChebyshevTDerivative, x::DT, i::Int) where {DT}
    local C = D.B
    local x̃ = promote_type(eltype(C), DT)(2x-1)
    @boundscheck i ≥ 0 && i < nbasis(C) || throw(BoundsError(C, i))
    _chebyshev_derivative(Val(1), i, x̃) * 2
end

"""
Evaluate derivative of Chebyshev polynomial of the second kind on the interval [0..1].
"""
function _eval(D::ChebyshevUDerivative, x::DT, i::Int) where {DT}
    local C = D.B
    local x̃ = promote_type(eltype(C), DT)(2x-1)
    @boundscheck i ≥ 0 && i < nbasis(C) || throw(BoundsError(C, i))
    _chebyshev_derivative(Val(2), i, x̃) * 2
end

Base.getindex(D::ChebyshevDerivative, x::Number, j::Integer) = _eval(D, x, j)
Base.getindex(D::ChebyshevDerivative, x::Number,  ::Colon) = [_eval(D, x, j) for j in eachindex(D.B)]
Base.getindex(D::ChebyshevDerivative, X::AbstractVector, j::Integer) = [_eval(D, x, j) for x in X]
Base.getindex(D::ChebyshevDerivative, X::AbstractVector,  ::Colon) = [_eval(D, x, j) for x in X, j in eachindex(D.B)]

Base.adjoint(C::Chebyshev) = Derivative(axes(C,1)) * C
