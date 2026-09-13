
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
struct Lagrange{T, BT, XT <: AbstractVector{T}} <: NodalBasis{T}
    b::BT
    x::XT

    denom::XT
    diffs::Matrix{T}

    function Lagrange{T}(x::XT) where {T, XT <: SVector}
        n = length(x)

        # a non-finite node poisons every difference it enters, and is `isequal` to nothing
        # at all, so it is caught before the differences are formed — and before their
        # product, so that the message names the fault that is actually present. A degenerate
        # product does not by itself say which: it is equally what an unrepresentable product
        # of perfectly good nodes gives, and reporting that one as a repeated node sends the
        # reader after the wrong thing.
        all(isfinite, x) || throw(ArgumentError(
            "the nodes of a Lagrange basis must be finite, got $(x)"))

        denom = zeros(T, n)
        diffs = zeros(T, n, n)

        for i in 1:n
            for j in 1:n
                diffs[i, j] = x[i] - x[j]
            end

            # the product of the node differences is what has to be invertible, so it is
            # tested rather than the node list: `allunique` compares with `isequal`, which
            # holds 0.0 and -0.0 to be distinct although their difference is zero.
            any(j -> j ≠ i && iszero(diffs[i, j]), 1:n) && throw(ArgumentError(
                "the nodes of a Lagrange basis must be distinct, got $(x)"))

            local p = prod(diffs[i, j] for j in 1:n if j ≠ i; init = one(T))

            (iszero(p) || !isfinite(p)) && throw(ArgumentError(
                "the nodes of a Lagrange basis are distinct and finite, but the product of " *
                "the differences from node $(i) is $(p) in $(T), so the denominator it " *
                "gives is not usable; rescale the nodes or widen the element type"))

            denom[i] = 1/p
        end

        sdenom = SVector{n}(denom)

        b = collect(y -> sdenom[j] *
                         mapreduce(i -> i ≠ j ? (y - x[i]) : one(T), *, eachindex(x))
        for j in eachindex(sdenom))

        new{T, typeof(b), typeof(x)}(b, x, sdenom, diffs)
    end

    Lagrange{T}(x::AbstractVector) where {T} = Lagrange{T}(SVector{length(x), T}(collect(x)))
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

_key(L::Lagrange) = (Lagrange, L.x)

## Derivative

"""
    LagrangeDerivative

The type of `Derivative(axes(l,1)) * l` for a [`Lagrange`](@ref) basis `l`, equivalently of
`l'`, evaluated from the node differences the basis caches. See
[`PolynomialBasisDerivative`](@ref) and [Derivatives](@ref).
"""
const LagrangeDerivative = QMul2{<:Derivative, <:Lagrange}

function _eval(D::LagrangeDerivative, x, j::Integer)
    local L = D.B
    local T = _evaltype(eltype(L), typeof(x))
    local d::T = 0

    # the product rule applied to ∏_{i≠j} (x - xᵢ) / (xⱼ - xᵢ): one summand per factor
    # differentiated, with the remaining factors left as they are
    for l in eachindex(L)
        if l ≠ j
            z = 1 / L.diffs[j, l]
            for i in eachindex(L)
                if i ≠ j && i ≠ l
                    z *= (x - L.x[i]) / L.diffs[j, i]
                end
            end
            d += z
        end
    end

    return d
end
