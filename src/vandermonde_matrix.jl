
using LinearAlgebra

@doc raw"""
    vandermonde_matrix(x)

The Vandermonde matrix of the nodes `x`, ``V_{ij} = x_i^{j-1}``.

`V` is the matrix of the monomial basis evaluated at the nodes, so `V * c` are the values at
the nodes of the polynomial with monomial coefficients `c`, and `V \ y` are the monomial
coefficients of the polynomial interpolating the values `y`.

```jldoctest
julia> vandermonde_matrix([0.0, 0.5, 1.0])
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 1.0  0.5  0.25
 1.0  1.0  1.0
```

See also [`vandermonde_matrix_inverse`](@ref) and [The Vandermonde matrix](@ref).
"""
vandermonde_matrix(x::Vector) = [x[i]^(j-1) for i in eachindex(x), j in eachindex(x)]
vandermonde_matrix(x::AbstractVector{T}) where {T} = vandermonde_matrix(collect(x))

@doc raw"""
    vandermonde_matrix_inverse(x)

The inverse of [`vandermonde_matrix`](@ref), computed in closed form.

`V` is factored as ``V^{-1} = U L`` with `L` lower and `U` upper triangular, whose entries
are known explicitly from the nodes: `L` holds the reciprocals of the products of node
differences that also appear in the Lagrange cardinal functions, and `U` is built by a
two-term recurrence. This avoids a general factorisation of a matrix that is notoriously
ill-conditioned, though the result is still limited by that conditioning.

`V⁻¹ * y` are the monomial coefficients of the polynomial interpolating the values `y` at
the nodes `x`.

```jldoctest
julia> vandermonde_matrix_inverse([0.0, 0.5, 1.0])
3×3 Matrix{Float64}:
  1.0   0.0   0.0
 -3.0   4.0  -1.0
  2.0  -4.0   2.0
```

The nodes must be distinct, as for [`Lagrange`](@ref); repeated nodes divide by zero.
"""
function vandermonde_matrix_inverse(x::Vector{T}) where {T}
    local n = length(x)
    local L = zeros(T, n, n)
    local U = Matrix{T}(I, n, n)

    L[1, 1] = one(T)
    for i in 2:n
        for j in 1:i
            p = one(T)
            for k in 1:i
                if k ≠ j
                    p *= (x[j] - x[k])
                end
            end
            L[i, j] = 1/p
        end
    end

    i = 1
    for j in (i + 1):n
        U[i, j] = - U[i, j - 1] * x[j - 1]
    end

    for i in 2:n
        for j in (i + 1):n
            U[i, j] = U[i - 1, j - 1] - U[i, j - 1] * x[j - 1]
        end
    end

    return U * L
end

function vandermonde_matrix_inverse(x::AbstractVector{T}) where {T}
    vandermonde_matrix_inverse(collect(x))
end
