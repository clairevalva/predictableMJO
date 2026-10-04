export projectdata, reconstructfromproj

"""Project coordinate-by-time data; delay coefficients run oldest to newest."""
function projectdata(x::AbstractMatrix, zetaR, mu, nEmb, nT_eig = nothing, nD = nothing)

    if nT_eig === nothing
        nT_eig = length(zetaR)
    end

    if nD === nothing
        nD = size(x, 1)
    end

    As = 1im*zeros(nD, nEmb)
    for q = 1:nEmb
        As[:, q] = sum(x[:,q:(nT_eig + q - 1)] .* conj.(zetaR)'.* mu, dims = 2)
    end

    return conj.(As)
end

"""Reconstruct coordinate-by-time values by averaging available delay coordinates."""
function reconstructfromproj(As::AbstractMatrix, zetaL, nEmb, nT_rec = nothing, nD = nothing)

    if nT_rec === nothing
        nT_rec = length(zetaL) + nEmb - 1
    end

    nT_eig = length(zetaL)

    if nD === nothing
        nD = size(As, 1)
    end

    ytest = 1im*zeros(nD, nT_rec)
    for j = 1:nT_rec
        if j < nEmb
            for k = 0:(j - 1)
                ytest[:, j] .+= As[:,1 + k] * zetaL[j - k]
            end
            ytest[:, j] ./= j
        elseif j > nT_eig
            for k = (j - nT_eig):(nEmb - 1)
                ytest[:, j] .+= As[:,1 + k] * zetaL[j - k]
            end

            ytest[:, j] ./= (nEmb - (j - nT_eig))
        else
            for k = 0:(nEmb - 1)
                ytest[:, j] .+= As[:,1 + k] * zetaL[j - k]
            end
            ytest[:, j] ./= nEmb
        end
    end

    return ytest
end
