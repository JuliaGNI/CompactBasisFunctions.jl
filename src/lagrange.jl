
import QuadratureRules: gauss_legendre_nodes, lobatto_legendre_nodes

@doc raw"""
    Lagrange(x)
    Lagrange{T}(x)

The Lagrange basis on the nodes `x`, generally taken in the interval ``[0,1]``,

```math
\ell_j(x) = \prod_{i \neq j} \frac{x - x_i}{x_j - x_i} ,
```

indexed from `1`, unlike the other three bases, which are indexed from `0`.

This is the *cardinal* basis of the nodes: ``\ell_j(x_i) = \delta_{ij}``, so the
coefficients of an expansion are the values of the function at the nodes, and no linear
system is needed to interpolate. The basis also forms a partition of unity,
``\sum_j \ell_j(x) = 1``.

```jldoctest
julia> l = Lagrange([0.0, 0.5, 1.0]);

julia> all(l[nodes(l)[i], j] ≈ (i == j) for i in 1:3, j in 1:3)   # ℓⱼ(xᵢ) = δᵢⱼ
true

julia> sum(l[0.3, j] for j in eachindex(l)) ≈ 1                   # partition of unity
true
```

The nodes must be **distinct** and **finite**, since the denominators are products of their
differences; anything that leaves such a product zero or non-finite throws an `ArgumentError`.
That covers a repeated node, the pair `0.0` and `-0.0`, whose difference is zero although the
two are not `isequal`, and a `NaN` or `Inf` node. They need not be sorted, nor confined to
``[0,1]``, although the declared domain is ``[0,1]``.

Which nodes to use matters: equidistant nodes make high-degree interpolation diverge near the
endpoints (the Runge phenomenon), so the Gauß-Legendre and Lobatto-Legendre node sets are
provided as [`LagrangeGauß`](@ref) and [`LagrangeLobatto`](@ref); see
[Choice of nodes](@ref).

The element type `T` is taken from the nodes, and all internal quantities are formed in it,
so an arbitrary-precision basis carries its full precision:

```jldoctest
julia> setprecision(BigFloat, 256) do
           l = Lagrange(BigFloat[0, 1//4, 1])
           abs(sum(l[BigFloat(1)/3, j] for j in eachindex(l)) - 1) < 1e-70
       end
true
```

See also [`Chebyshev`](@ref) for the other nodal basis, and [Lagrange basis](@ref) for the
full discussion.
"""
struct Lagrange{T, BT, XT <: AbstractVector{T}} <: Basis{T}
    b::BT
    x::XT

    denom::XT
    diffs::Matrix{T}

    function Lagrange{T}(x::XT) where {T, XT <: SVector}
        n = length(x)

        denom = zeros(T, n)
        diffs = zeros(T, n, n)

        for i in eachindex(x)
            local p = one(T)
            for j in eachindex(x)
                diffs[i,j] = x[i] - x[j]
                if i ≠ j
                    p *= diffs[i,j]
                end
            end

            # the product of the node differences is what has to be invertible, so it is
            # tested rather than the node list: `allunique` compares with `isequal`, which
            # holds 0.0 and -0.0 to be distinct although their difference is zero, and lets
            # a lone NaN through to poison every difference silently
            iszero(p) && throw(ArgumentError(
                "the nodes of a Lagrange basis must be distinct, got $(x)"))
            isfinite(p) || throw(ArgumentError(
                "the nodes of a Lagrange basis must be finite, got $(x)"))

            denom[i] = 1/p
        end

        sdenom = SVector{n}(denom)

        b = collect(y -> sdenom[j] * mapreduce(i -> i ≠ j ? (y - x[i]) : one(T), *, eachindex(x)) for j in eachindex(sdenom))

        new{T, typeof(b), typeof(x)}(b, x, sdenom, diffs)
    end

    Lagrange{T}(x::Vector) where {T} = Lagrange{T}(SVector{length(x),T}(x))
    Lagrange{T}(x::AbstractVector) where {T} = Lagrange{T}(collect(x))
end

Lagrange(x::AbstractVector{T}) where {T} = Lagrange{T}(x)

"""
    LagrangeGauß(n)

The [`Lagrange`](@ref) basis on the `n` Gauß-Legendre nodes of ``[0,1]``.

These lie strictly inside the interval, which suits a basis whose expansion is integrated
rather than matched at the boundary. Compare [`LagrangeLobatto`](@ref), whose nodes include
the endpoints.

```jldoctest
julia> nodes(LagrangeGauß(2)) ≈ [(1 - 1/sqrt(3)) / 2, (1 + 1/sqrt(3)) / 2]
true
```
"""
LagrangeGauß(n) = Lagrange(gauss_legendre_nodes(n))

"""
    LagrangeLobatto(n)

The [`Lagrange`](@ref) basis on the `n` Lobatto-Legendre nodes of ``[0,1]``.

These include both endpoints, so an expansion has coefficients that are the boundary values
themselves — what a method needs when it has to impose or read off conditions there.
Compare [`LagrangeGauß`](@ref), whose nodes lie strictly inside.

```jldoctest
julia> nodes(LagrangeLobatto(3)) == [0.0, 0.5, 1.0]
true
```
"""
LagrangeLobatto(n) = Lagrange(lobatto_legendre_nodes(n))

(L::Lagrange)(x::Number, j::Integer) = L.b[j](x)

basis(L::Lagrange) = L.b
nodes(L::Lagrange) = L.x
nbasis(L::Lagrange) = length(basis(L))
nnodes(L::Lagrange) = length(nodes(L))
order(L::Lagrange)  = nnodes(L)
degree(L::Lagrange) = nnodes(L) - 1

Base.eltype(::Lagrange{T}) where {T} = T
Base.eachindex(L::Lagrange) = eachindex(L.b)
Base.axes(L::Lagrange) = (Inclusion(0..1), eachindex(L))
ContinuumArrays.grid(L::Lagrange) = nodes(L)

Base.hash(L::Lagrange, h::UInt) = hash(L.x, h)
Base.:(==)(L1::Lagrange, L2::Lagrange) = (L1.x == L2.x)
Base.isequal(L1::Lagrange{T1}, L2::Lagrange{T2}) where {T1,T2} = (T1 == T2 && L1 == L2)
Base.isapprox(L1::Lagrange, L2::Lagrange; kwargs...) = isapprox(L1.x, L2.x; kwargs...)

Base.getindex(L::Lagrange, x::Number, j::Integer) = L(x,j)
Base.getindex(L::Lagrange, x::Number,  ::Colon) = [b(x) for b in L.b]
Base.getindex(L::Lagrange, X::AbstractVector, j::Integer) = L.(X,j)
Base.getindex(L::Lagrange, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in L.b]


## Derivative

@simplify *(D::Derivative, L::Lagrange) = Mul(D,L)

"""
    LagrangeDerivative

The type of `Derivative(axes(l,1)) * l` for a [`Lagrange`](@ref) basis `l`, equivalently of
`l'`.

A lazy product: it stores the basis and evaluates the derivative on indexing, from the node
differences the basis caches. See [Derivatives](@ref).
"""
const LagrangeDerivative = QMul2{<:Derivative,<:Lagrange}

function _eval(D::LagrangeDerivative, x::DT, j::Int) where {DT}
    local L = D.B
    local T = promote_type(eltype(L), DT)
    local d::T = 0

    for l in eachindex(L)
        if l ≠ j
            z = 1 / L.diffs[j,l]
            for i in eachindex(L)
                i ≠ j && i ≠ l ? z *= (x - L.x[i]) / L.diffs[j,i] : nothing
            end
            d += z
        end
    end
    return d
end

Base.getindex(D::LagrangeDerivative, x::Number, j::Integer) = _eval(D, x, j)
Base.getindex(D::LagrangeDerivative, x::Number,  ::Colon) = [_eval(D, x, j) for j in eachindex(D.B)]
Base.getindex(D::LagrangeDerivative, X::AbstractVector, j::Integer) = [_eval(D, x, j) for x in X]
Base.getindex(D::LagrangeDerivative, X::AbstractVector,  ::Colon) = [_eval(D, x, j) for x in X, j in eachindex(D.B)]

Base.adjoint(L::Lagrange) = Derivative(axes(L,1)) * L
