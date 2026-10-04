include(joinpath(@__DIR__, "utils", "utils.jl"))
include(joinpath(@__DIR__, "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))

"""Build the cone-kernel diffusion basis from the prepared normalized features."""
function compute_nlsa(; output = paper_cache("init_NLSA.h5"), overwrite = false, settings = PAPER, feature_file = FEATURE_FILE)
    require_files([feature_file])
    isfile(output) && !overwrite && error("Output exists: $output; pass --force")
    features = h5read(feature_file, "olr_u") # variable × longitude × day
    X = reshape(features, :, size(features, 3))[:, 1:settings.train_last]
    Xemb = ConsistentKoopman.delayembed(X, settings.nemb)
    D, DN = ConsistentKoopman.distNN(X, settings.neighbors, settings.nemb)
    DN = ConsistentKoopman.symDist(DN)
    sigma = ConsistentKoopman.coneBandwidths(Xemb, D)
    trimmed_D, trimmed_DN = D[2:end, 2:end], DN[2:end, 2:end]
    bw, _ = ConsistentKoopman.tuneBandwidth(trimmed_D, sigma, 2 .^ range(-40, 40; length = 100))
    W = ConsistentKoopman.makeW(trimmed_D, bw, sigma)
    W[.!trimmed_DN] .= 0
    P = ConsistentKoopman.normW(W)
    kappa, phi, _ = ConsistentKoopman.computeDiffusionEig(P, settings.diffusion_modes)
    return write_cache(output; overwrite) do f
        for (key, value) in (("W", W), ("bw", bw), ("phi", phi), ("kappa", kappa))
            f[key] = value
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_nlsa(; overwrite = "--force" in ARGS)
end
