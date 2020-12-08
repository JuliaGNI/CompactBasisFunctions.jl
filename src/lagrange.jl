"""
Lagrange basis on the interval [0..1].
"""
struct Lagrange{T, BT, XT <: AbstractVector{T}} <: Basis{T}
    b::BT
    x::XT

    denom::XT
    diffs::Matrix{T}
    vdminv::Matrix{T}

    function Lagrange{T}(x̂) where {T}
        n = length(x̂)
        x = SVector{n}(x̂)

        denom = zeros(n)
        diffs = zeros(n,n)

        for i in eachindex(x)
            local p = 1
            for j in eachindex(x)
                diffs[i,j] = x[i] - x[j]
                if i ≠ j
                    p *= diffs[i,j]
                end
            end
            denom[i] = 1/p
        end

        sdenom = SVector{n}(denom)

        b = collect(y -> sdenom[j] * mapreduce(i -> i ≠ j ? (y - x[i]) : one(T), *, eachindex(x)) for j in eachindex(sdenom))

        new{T, typeof(b), typeof(x)}(b, x, sdenom, diffs, vandermonde_matrix_inverse(x))
    end
end

Lagrange(x::AbstractVector{T}) where {T} = Lagrange{T}(x)

# LagrangeGauß(n) = Lagrange()
# LagrangeLobatto(n) = Lagrange()

(L::Lagrange)(x::Number, j::Integer) = L.b[j](x)

nodes(L::Lagrange) = L.x
nnodes(L::Lagrange) = length(nodes(L))
nbasis(L::Lagrange) = nnodes(L)
degree(L::Lagrange) = nnodes(L) - 1
order(L::Lagrange) = nnodes(L)
eachnode(L::Lagrange) = Base.OneTo(nnodes(L))

Base.axes(L::Lagrange) = (Inclusion(0..1), eachnode(L))
ContinuumArrays.grid(L::Lagrange) = nodes(L)

Base.hash(L::Lagrange, h::UInt) = hash(L.x, h)
Base.:(==)(L1::Lagrange, L2::Lagrange) = (L1.x == L2.x)
Base.isequal(L1::Lagrange{T1}, L2::Lagrange{T2}) where {T1,T2} = (T1 == T2 && L1 == L2)

Base.getindex(L::Lagrange, x::Number, j::Integer) = L(x,j)
Base.getindex(L::Lagrange, x::Number,  ::Colon) = [b(x) for b in L.b]
Base.getindex(L::Lagrange, X::AbstractVector, j::Integer) = L.(X,j)
Base.getindex(L::Lagrange, X::AbstractVector,  ::Colon) = [b(x) for x in X, b in L.b]


## Derivative

function _eval_derivative(L::Lagrange{LT}, x::DT, j::Int) where {LT,DT}
    local T = promote_type(LT,DT)
    local d::T = 0

    for l in eachnode(L)
        if l ≠ j
            z = 1 / L.diffs[j,l]
            for i in eachnode(L)
                i ≠ j && i ≠ l ? z *= (x - L.x[i]) / L.diffs[j,i] : nothing
            end
            d += z
        end
    end
    return d
end

@simplify *(D::Derivative, L::Lagrange) = Mul(D,L)

const LagrangeDerivative = QMul2{<:Derivative,<:Lagrange}

Base.getindex(D::LagrangeDerivative, x::Number, j::Integer) = _eval_derivative(D.B, x, j)
Base.getindex(D::LagrangeDerivative, x::Number,  ::Colon) = [_eval_derivative(D.B, x, j) for j in eachnode(D.B)]
Base.getindex(D::LagrangeDerivative, X::AbstractVector, j::Integer) = [_eval_derivative(D.B, x, j) for x in X]
Base.getindex(D::LagrangeDerivative, X::AbstractVector,  ::Colon) = [_eval_derivative(D.B, x, j) for x in X, j in eachnode(D.B)]
