export delayembed, distNN, rmNN, xcorr, sortautocorr

"""Delay windows are columns, with the newest coordinate first."""
function delayembed(X, numdelays::Int)
        numdelays -= 1
        if numdelays == 0
            return X
        end

        n = size(X,1);
        m = size(X,2);
        Xemb = zeros(typeof(X[1]), n*(1 + numdelays), m - numdelays);
        for del in 0:numdelays
            Xemb[del*n .+ (1:n), :] = X[:, (numdelays - del + 1):(m -  del)];
        end

        return Xemb
end

"""Delay-window distances and the nearest-neighbor mask."""
function distNN(X::Matrix, NN::Integer, emb::Integer; usenorm = euclidean)
    _, nT = size(X)
    D = pairwise(usenorm, X, dims=2)

    nT_emb = nT - emb + 1
    D_emb = zeros(Float64, nT_emb, nT_emb)

    for k = 0:(emb - 1)
        D_emb += D[(emb - k):(end - k), (emb - k):(end - k)].^2
    end
    D_emb = sqrt.(D_emb)

    if NN > 0
        N = rmNN(D_emb, NN)
        return D_emb, N
    else
        N = ones(Bool, nT_emb, nT_emb)
        return D_emb, N
    end
end

"""Select the nearest neighbors of each sample."""
function rmNN(D::Matrix, NN::Integer = 0)
    nT = size(D, 1)
    N = ones(Bool, nT, nT)
    if NN == 0
        return D
    else
        sortsD = sortperm(D, dims = 2)
        N[sortsD[:, (NN + 1):end]] .= false

    end

    return N
end

"""Complex lag correlation with the original normalization."""
function xcorr(x::AbstractVector, y::AbstractVector, lag::Int; coeff = "innerprod")
    rawcorr = sum(x[(lag + 1):end] .* conj.(y[1:(end - lag)]))

    if coeff == "innerprod"
        xnorm = sum(x[(lag + 1):end] .* conj.(x[(lag + 1):end]))
        ynorm = sum(y[1:(end - lag)] .* conj.(y[1:(end - lag)]))

        usec = sqrt(xnorm * ynorm)

    elseif coeff == "normalized"
        xnorm = sum(x .* conj.(x))
        ynorm = sum(y .* conj.(y))
        usec = usec = sqrt(xnorm * ynorm)
    else
        usec = 1
    end

    return rawcorr / usec
end

"""Rank modes by their maximal autocorrelation error."""
function sortautocorr(ζ::Matrix{ComplexF64}, frequencies::Vector{Float64}, nLags::Int64, dt::Float64; Flim::Float64 = 30., returnall = false)
    L = size(ζ, 2)
    acs = 1im*zeros(nLags + 1, L)
    goalacs = 1im*zeros(nLags + 1, L)

    for j in 1:L
        goalacs[:, j] = exp.(-1im*frequencies[j]*(0:nLags)*dt)
        for lag in 0:nLags
            acs[lag + 1, j] = xcorr(ζ[:,j], ζ[:,j], lag, coeff = "innerprod")
        end
    end

    eps_t = real(1 .- abs.(goalacs .* acs))
    max_eps = maximum(eps_t, dims = 1)[1,:] + (abs.(frequencies) .> Flim)
    sortinds = sortperm(vec(max_eps))
    if returnall
        return max_eps, sortinds, acs, goalacs
    else
        return max_eps, sortinds
    end

end
