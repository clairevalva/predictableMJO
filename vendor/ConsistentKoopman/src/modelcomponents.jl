export normW, computeDiffusionEig

"""Normalize the kernel and form its symmetric diffusion matrix."""
function normW(X::Union{Matrix{Float64}, SparseMatrixCSC{Float64, Int64}})

    D = sum(X, dims = 2)[:]
    Dinv = Diagonal(D.^(-1))
    S = sum(X ./ (D'), dims = 2)[:]
    Sneghalf = Diagonal(S.^(-1/2))

    K̂ = Dinv * X * Sneghalf
    K̃ = K̂ * (K̂')

    return K̃
end

"""Leading diffusion eigenpairs, scaled by the first eigenvector entry."""
function computeDiffusionEig(K, L::Integer = 0)

    if K isa SparseMatrixCSC{Float64, Int64}
        K = Matrix(K)
    end
    η, φ = eigen(K, sortby = x -> -real(x))
    if L > 0
        η = η[1:L]
        φ = φ[:,1:L]
    end
    w = φ[:,1].^2
    normφ = φ ./ φ[1,1]
    return η, normφ, w
end
