export centeredDiff, makeVf, makeΔ, diffV, doKoopman_diff

"""Fourth-order centered time derivative, with zero boundary values."""
function centeredDiff(f::AbstractVector, dt::Real)

    Vf = 1/12*(circshift(f, 2) - circshift(f, -2)) +  2/3*(circshift(f, -1) - circshift(f, 1))
    Vf[1:2] .= 0
    Vf[(end - 1):end] .= 0

    return Vf/dt
end

"""Antisymmetric Galerkin matrix of the time derivative."""
function makeVf(fs::AbstractMatrix, dt::Real)
    Vf_unsym = mapslices(x -> centeredDiff(x, dt), fs, dims = 1)
    N = size(fs, 1)
    Ṽ = (fs' * Vf_unsym) / N
    Ṽ_sym = (Ṽ - Ṽ')/2

    return Ṽ_sym
end

"""Diffusion regularization normalized by the second eigenvalue."""
function makeΔ(κs::AbstractVector)
    E = 1 ./κs .- 1
    E ./= E[2]
    return Diagonal(E)
end

"""Regularized finite-difference generator."""
function diffV(V::AbstractMatrix, Δ::AbstractMatrix, ϵ::Real)
    return V - ϵ*Δ
end

"""Generator eigenpairs and their time-series representation."""
function doKoopman_diff(fs::AbstractMatrix, κs::AbstractVector, dt::Real, ϵ::Real)
    V = makeVf(fs, dt)
    Δ = makeΔ(κs)
    W = diffV(V, Δ, ϵ)

    λ, u = eigen(W, sortby = x -> -1*real(x))
    nu = λ / (2*π)
    gs = 1im * zeros(size(fs))
    L = size(fs, 2)

    Threads.@threads for j = 1:L
        for i = 1:L
            gs[:, j] .+= u[i, j] * fs[:, i]
        end
    end

    return λ, u, nu, gs
end
