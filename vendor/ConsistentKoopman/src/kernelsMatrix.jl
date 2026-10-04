export makeSϵ, tuneBandwidth, coneBandwidths, makeW

"""Log kernel sum used to select the bandwidth."""
function makeSϵ(D::AbstractMatrix, σ)
    function estSϵ(γ)
        γ = exp(γ)
        kentries = exp.(-1*D.^2 ./ (γ^2 * σ.^2))
        S = sum(kentries)
        return log(S)
    end
end

"""Select the bandwidth by the derivative of the log kernel sum."""
function tuneBandwidth(D::AbstractMatrix, σ, testγ::AbstractVector)

    S = makeSϵ(D, σ)
    S_grad(γ) = ForwardDiff.derivative(S, γ)
    Sgs = S_grad.(testγ)

    if !isfinite(sum(Sgs))
        println("seems like an overflow error... σ too big?")
        N = findfirst(x -> isnan(x), Sgs)

        Sgs = Sgs[1:N-1]
    end

    bestj = argmax(Sgs)
    bestϵ = exp.(testγ[bestj])
    m̂ = Sgs[bestj]

    return bestϵ, m̂
end

"""Cone bandwidth matrix after dropping the first transition sample."""
function coneBandwidths(X::AbstractMatrix, D::AbstractMatrix; ζ = 0.995)
    _, nT = size(X)
    cθ1 = zeros(nT - 1, nT - 1)
    cθ2 = zeros(nT - 1, nT - 1)

    vx = X[:, 2:nT] .- X[:, 1:(nT - 1)]
    vx = vx ./ norm.(eachcol(vx), 2)'

    Threads.@threads for i = 1:(nT - 1)
        if (i % 100) == 0
            println("iter ", i, ", ", i / nT, "%")
            flush(stdout)
        end
        for j = 1:(i-1)
            diffxy = (X[:, i + 1] - X[:, j + 1]) / D[i + 1, j + 1]
            cθ1[i, j] = -1*vx[:,i]' * diffxy
            cθ2[i, j] = vx[:,j]' * diffxy
        end
    end

    cθ1 .+= cθ1'
    cθ2 .+= cθ2'

    σinv = sqrt.((1 .- ζ*cθ1) .* (1 .- ζ*cθ2)).^-1
    return σinv
end

"""Evaluate the kernel using the selected global and pairwise bandwidths."""
function makeW(D::AbstractMatrix, γ::Real, σ)
    K = exp.(-1*D.^2 ./ (γ^2 * σ.^2))
    return K
end
