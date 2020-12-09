
using LinearAlgebra

vandermonde_matrix(x::Vector) = hcat([x.^(i-1) for i in 1:length(x)]...)
vandermonde_matrix(x::AbstractVector{T}) where {T} = vandermonde_matrix(collect(x))


function vandermonde_matrix_inverse(x::Vector{T}) where {T}
    local n = length(x)
    local L::Matrix{T} = zeros(n,n)
    local U::Matrix{T} = Matrix{T}(I, n, n)

    L[1,1] = 1
    for i in 2:n
        for j in 1:i
            p = 1
            for k in 1:i
                if k ≠ j
                    p *= (x[j] - x[k])
                end
            end
            L[i,j] = 1/p
        end
    end

    i = 1
    for j in i+1:n
        U[i,j] = - U[i,j-1] * x[j-1]
    end

    for i in 2:n
        for j in i+1:n
            U[i,j] = U[i-1,j-1] - U[i,j-1] * x[j-1]
        end
    end

    return *(U,L)
end

vandermonde_matrix_inverse(x::AbstractVector{T}) where {T} = vandermonde_matrix_inverse(collect(x))
