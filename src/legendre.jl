
"""
Legendre polynomial on the interval [-1..+1].

Evaluated by iterating Bonnet's recurrence ``j P_j = (2j-1) x P_{j-1} - (j-1) P_{j-2}``
upwards. The recursive formulation of the same recurrence descends into two subproblems per
step and recomputes shared subtrees, which costs `O(φʲ)` evaluations rather than `O(j)`.
"""
@inline function _legendre(j::Integer, x::T) where {T}
    j < 0 && return zero(T)
    j == 0 && return one(T)

    local p₂ = one(T)
    local p₁ = x

    for k in 2:j
        p₂, p₁ = p₁, ((2k-1) * p₁ * x - (k-1) * p₂) / k
    end

    return p₁
end

@doc raw"""
    Legendre(n)
    Legendre(T, n)

The Legendre basis of `n` functions, i.e. of degree ``p = n-1``, on the interval ``[0,1]``,

```math
L_j(x) = \sqrt{2j+1} \, P_j(2x-1) , \qquad j = 0, \dots, p ,
```

indexed from `0`. `T` is the element type and defaults to `Float64`.

The Legendre polynomials ``P_j`` are defined on ``[-1,+1]``, so the argument is shifted by
``\tilde{x} = 2x-1``; see [The reference interval](@ref).

The basis is *modal*: its coefficients are those of an expansion rather than values at
points, so it has no [`nodes`](@ref) and no `grid`, and those accessors throw.

The factor ``\sqrt{2j+1}`` makes the basis **orthonormal** on ``[0,1]``: since
``\int_0^1 P_i(2x-1) P_j(2x-1) \, dx = \delta_{ij} / (2j+1)``,

```math
\int_0^1 L_i(x) \, L_j(x) \, dx = \delta_{ij} ,
```

so the mass matrix is the identity and the coefficients of a projection are just the inner
products against the basis functions.

```jldoctest
julia> l = Legendre(3);

julia> l[0.5, 0], l[0.5, 1], l[0.5, 2]      # the midpoint is x̃ = 0
(1.0, 0.0, -1.118033988749895)

julia> l[1.0, 1] ≈ sqrt(3)                  # L₁(1) = √3 P₁(1) = √3
true

julia> using QuadratureRules

julia> quad = GaussLegendreQuadrature(8);

julia> sum(weights(quad)[k] * l[nodes(quad)[k], 1]^2 for k in eachindex(nodes(quad))) ≈ 1
true
```

The derivative follows from differentiating Bonnet's recurrence, with the chain-rule factor
`2` from the shift, and is obtained as `Derivative(axes(l,1)) * l`; see
[Derivatives](@ref).

See also [`Bernstein`](@ref) for the other modal basis, and [Legendre basis](@ref) for the
full discussion.
"""
struct Legendre{T} <: ModalBasis{T}
    n::Int
end

Legendre(::Type{T}, n::Integer) where {T} = Legendre{T}(n)
Legendre(n::Integer) = Legendre(Float64, n)

_key(L::Legendre) = (Legendre, L.n)

function _eval(L::Legendre, x, j::Integer)
    local T = _evaltype(eltype(L), typeof(x))
    # the conversion precedes the shift onto [-1,1], and the normalisation factor is formed
    # in the same arithmetic as the recurrence, cf. `_evaltype`
    _legendre(j, 2 * T(x) - 1) * sqrt(T(2j+1))
end

## Derivative

"""
Derivative of Legendre polynomial on the interval [-1..+1].

Obtained by differentiating Bonnet's recurrence term by term,

```math
j P_j' = (2j-1) P_{j-1} + (2j-1) x P_{j-1}' - (j-1) P_{j-2}' ,
```

and iterated upwards alongside ``P_j`` itself, which the right-hand side needs.
"""
@inline function _legendre_derivative(j::Integer, x::T) where {T}
    j <= 0 && return zero(T)
    j == 1 && return one(T)

    local p₂ = one(T)    # P_{k-2}
    local p₁ = x         # P_{k-1}
    local d₂ = zero(T)   # P_{k-2}'
    local d₁ = one(T)    # P_{k-1}'

    for k in 2:j
        # evaluated before p₁ advances, since the derivative needs P_{k-1}
        local d = ((2k-1) * p₁ + (2k-1) * d₁ * x - (k-1) * d₂) / k
        p₂, p₁ = p₁, ((2k-1) * p₁ * x - (k-1) * p₂) / k
        d₂, d₁ = d₁, d
    end

    return d₁
end

"""
    LegendreDerivative

The type of `Derivative(axes(l,1)) * l` for a [`Legendre`](@ref) basis `l`, equivalently of
`l'`. See [`PolynomialBasisDerivative`](@ref) and [Derivatives](@ref).
"""
const LegendreDerivative = QMul2{<:Derivative, <:Legendre}

function _eval(D::LegendreDerivative, x, j::Integer)
    local T = _evaltype(eltype(D.B), typeof(x))
    # the sqrt(2j+1) is part of the basis function, so it belongs to its derivative too;
    # the 2 is the chain rule of the shift onto [0,1]
    _legendre_derivative(j, 2 * T(x) - 1) * 2 * sqrt(T(2j+1))
end
