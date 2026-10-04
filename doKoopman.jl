include(joinpath(@__DIR__, "utils", "utils.jl"))
include(joinpath(@__DIR__, "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))

function compute_koopman(; input = paper_cache("init_NLSA.h5"),
    output = paper_cache("init_Keigen.h5"), overwrite = false, settings = PAPER)
    require_files([input])
    phis = h5read(input, "phi")[:, 1:settings.nphi]
    lambdas = h5read(input, "kappa")[1:settings.nphi]
    omegas, c, _, _ = ConsistentKoopman.doKoopman_diff(phis, lambdas, settings.dt, settings.regularization)
    cR = (c * c') \ c
    zetaL, zetaR = phis * c, phis * cR
    return write_cache(output; overwrite) do f
        f["omegas"] = omegas
        f["c"] = c
        f["zetaL"] = zetaL
        f["zetaR"] = zetaR
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_koopman(; overwrite = "--force" in ARGS)
end
