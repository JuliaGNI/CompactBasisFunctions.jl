
"""
Chebyshev polynomial of the first (`kind = 1`) or second (`kind = 2`) kind on the interval
[-1..+1].

Both kinds obey the same three-term recurrence ``\\phi_j = 2x \\phi_{j-1} - \\phi_{j-2}`` and
differ only in where it starts, ``T_1 = x`` against ``U_1 = 2x`` — the `kind * x` below. It is
iterated upwards; the recursive formulation descends into two subproblems per step and
recomputes shared subtrees, which costs `O(φʲ)` evaluations rather than `O(j)`.
"""
@inline function _chebyshev(::Val{kind}, j::Integer, x::T) where {kind, T}
    j < 0 && return zero(T)
    j == 0 && return one(T)

    local p₂ = one(T)
    local p₁ = kind * x

    for _ in 2:j
        p₂, p₁ = p₁, p₁ * 2x - p₂
    end

    return p₁
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
struct Chebyshev{kind, T, BT} <: NodalBasis{T}
    b::BT
    x::Vector{T}

    function Chebyshev{kind, T}(n::Integer) where {kind, T}
        p = n-1
        # chebyshev_nodes returns the points on [0,1]; shift_nodes widens integer
        # element types, so convert back to keep the nodes in T
        x = convert(Vector{T}, chebyshev_nodes(T, n, Val(kind)))
        # evaluated in the wider of T and the argument type, as the derivatives below are.
        # The conversion precedes the shift onto [-1,1]: 2y-1 rounds in the argument's type,
        # and widening the rounded result would freeze that error in.
        b = OffsetArray(
            [y -> _chebyshev(Val(kind), i, 2 * _evaltype(T, typeof(y))(y) - 1)
             for i in 0:p], 0:p)
        new{kind, T, typeof(b)}(b, x)
    end

    Chebyshev{kind}(::Type{T}, n::Integer) where {kind, T} = Chebyshev{kind, T}(n)
    Chebyshev{kind}(n::Integer) where {kind} = Chebyshev{kind, Float64}(n)
end

Chebyshev(::Type{T}, n::Integer, ::Val{kind}) where {kind, T} = Chebyshev{kind, T}(n)
Chebyshev(n::Integer, ::Val{kind}) where {kind} = Chebyshev(Float64, n, Val(kind))

const ChebyshevT = Chebyshev{1}
const ChebyshevU = Chebyshev{2}

"""
The kind of a Chebyshev basis, as the `Val` that the evaluators dispatch on.
"""
_kind(::Chebyshev{kind}) where {kind} = Val(kind)

_key(C::Chebyshev{kind}) where {kind} = (Chebyshev{kind}, C.x)

## Derivative

"""
    ChebyshevDerivative

The type of `Derivative(axes(b,1)) * b` for a [`Chebyshev`](@ref) basis `b` of either kind,
equivalently of `b'`. The two kinds need different formulas, and
[`ChebyshevTDerivative`](@ref) and [`ChebyshevUDerivative`](@ref) name them apart. See
[`PolynomialBasisDerivative`](@ref) and [Derivatives](@ref).
"""
const ChebyshevDerivative = QMul2{<:Derivative, <:Chebyshev}

"""
    ChebyshevTDerivative

[`ChebyshevDerivative`](@ref) narrowed to the first kind, evaluated through
``T_j' = j \\, U_{j-1}``.
"""
const ChebyshevTDerivative = QMul2{<:Derivative, <:ChebyshevT}

"""
    ChebyshevUDerivative

[`ChebyshevDerivative`](@ref) narrowed to the second kind, evaluated through the
differentiated recurrence rather than the closed form, which is singular at the endpoints.
"""
const ChebyshevUDerivative = QMul2{<:Derivative, <:ChebyshevU}

"""
Derivative of the Chebyshev polynomial of the first kind on the interval [-1..+1],
via ``T_j' = j \\, U_{j-1}``.
"""
@inline function _chebyshev_derivative(::Val{1}, j::Integer, x::T) where {T}
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
@inline function _chebyshev_derivative(::Val{2}, j::Integer, x::T) where {T}
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

function _eval(D::ChebyshevDerivative, x, j::Integer)
    # converted before the shift onto [-1,1], so that 2x-1 is not rounded in the argument's
    # type and the rounding then widened along with it, cf. `_evaltype`
    local x̃ = 2 * _evaltype(eltype(D.B), typeof(x))(x) - 1
    _chebyshev_derivative(_kind(D.B), j, x̃) * 2
end
