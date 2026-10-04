export conebw

"""Cone bandwidth for two consecutive embedded states."""
function conebw(x::Vector{Float64}, y::Vector{Float64},
     xpre::Vector{Float64}, ypre::Vector{Float64}; ζ = 0.995)

    if (x == y) & (xpre == ypre)
        return 1
    end

    vx = x .- xpre
    vy = y .- ypre

    vx = vx / norm(vx, 2)
    vy = vy / norm(vy, 2)

    diffxy = x .- y
    diffxy = diffxy / norm(diffxy, 2)

    cθ1 = -vx' * diffxy
    cθ2 = vy' * diffxy

    return sqrt((1 - ζ*cθ1) * (1 - ζ*cθ2))
end
