# Out-of-sample evaluation, using the archived paper setting of 2,000 neighbors.
using HDF5
using Dates
using Statistics
using LinearAlgebra

include(joinpath(@__DIR__, "utils", "utils.jl"))
include(joinpath(@__DIR__, "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))

function compute_ose(; output = OSE_FILE, overwrite = false, settings = PAPER,
    feature_file = FEATURE_FILE, nlsa_file = paper_cache("init_NLSA.h5"),
    koopman_file = paper_cache("init_Keigen.h5"))
    isfile(output) && !overwrite && error("Output exists: $output; pass --force")
    tStart = Dates.value(settings.ose_start - settings.start_date) + 1
    tEnd = Dates.value(settings.ose_end - settings.start_date) + 1
    nEmb, p_NN = settings.nemb, settings.ose_neighbors
    sn_NLSA, sn_K = nlsa_file, koopman_file
    # get correct data for prediction, adjusting for embedding length
    X_load = h5read(feature_file, "olr_u")
    X = reshape(X_load, :, size(X_load, 3))

    X_new_e = X[:, (tStart - nEmb):tEnd]
    Xemb_new_e = ConsistentKoopman.delayembed(X_new_e, nEmb)
    # this should have first day as tStart - nEmb, and last day as tEnd, so that the physical time mapped to j=1 is tStart and has no future info leak.
    # appears to be true as X_new_e has size (720, 1891) and Xemb_new_e has size (720*nEmb, 1891 - nEmb + 1) = (720*64, 1891 - 64 + 1) = (46080, 1828).

    dates_new_e = settings.start_date .+ Day.(((tStart - nEmb):tEnd) .- 1)
    inds_new_e = reshape(1:length(dates_new_e), 1, :)
    inds_emb = ConsistentKoopman.delayembed(inds_new_e, nEmb)
    j_check = min(147, length(tStart:tEnd))
    ypre_dates_actual = dates_new_e[Int.(inds_emb[:, j_check])]
    y_dates_actual = dates_new_e[Int.(inds_emb[:, j_check + 1])]
    current_date_j = settings.start_date + Day(tStart + j_check - 2)

    @assert maximum(ypre_dates_actual) == current_date_j - Day(1) "Unexpected ypre date alignment in OSE transition"
    @assert maximum(y_dates_actual) == current_date_j "Unexpected y date alignment in OSE transition"
    println(
        "Delayembed transition date check passed. Target date: ",
        current_date_j,
        " | Max ypre date: ",
        maximum(ypre_dates_actual),
        " | Max y date: ",
        maximum(y_dates_actual),
    )

    # make and evaluate normalized kernel function
    tStart_o = 1
    tEnd_o = settings.train_last

    X_o = X[:, tStart_o:tEnd_o]
    if nEmb > 1
        Xemb_o = ConsistentKoopman.delayembed(X_o, nEmb)
        nT_o = size(Xemb_o, 2)
    end
    nT_o -= 1

    W = h5read(sn_NLSA, "W")
    bw = h5read(sn_NLSA, "bw")

    p = ConsistentKoopman.makeNormKernel_cone(W, Xemb_o, nT_o, bw, p_NN)

    nT_p = length(tStart:tEnd)

    p_evals_all = zeros(nT_p, nT_o)
    for j = 1:nT_p
        j % 100 == 1 && println("Computing OSE row ", j, " / ", nT_p)
        flush(stdout)
        y = Xemb_new_e[:, j + 1]
        ypre = Xemb_new_e[:, j]

        p_evals_all[j, :] = p(y, ypre, verbose2 = false)
    end
    # OSE row j is the cone-kernel transition ending at physical time tStart + j - 1.

    # Evaluate the diffusion basis and transform it to Koopman coordinates.
    phis, lambdas = h5read(sn_NLSA, "phi"), h5read(sn_NLSA, "kappa")

    # compute corresponding phi, zeta
    phi_calc_all =  zeros(settings.nphi, nT_p)

    for k = 1:settings.nphi
        for j = 1:nT_p
            phi_calc_all[k, j] =  ConsistentKoopman.evalPhi(p_evals_all[j, :], phis[:, k], lambdas[k])
        end
    end


    c = h5read(sn_K, "c")
    c = c[1:settings.nphi, 1:settings.nphi]
    cR = (c * c') \ c;


    # have zetaL = phis * c and zetaR = phis * cR
    nZ = min(51, settings.nphi)
    zetaL_vals = 1im*zeros(nZ, nT_p)
    zetaR_vals = 1im*zeros(nZ, nT_p)

    for j = 1:nT_p
        phi_vals = phi_calc_all[:, j]
        for k = 1:nZ
            zetaL_vals[k, j] = sum(phi_vals .* c[:, k])
            zetaR_vals[k, j] = sum(phi_vals .* cR[:, k])
        end
    end


    return write_cache(output; overwrite) do f
        f["p_evals"] = p_evals_all
        f["phi_vals"] = phi_calc_all
        f["zeta_vals"] = zetaL_vals
        f["zetaR_vals"] = zetaR_vals
        f["dates"] = string.(settings.ose_start:Day(1):settings.ose_end)
        f["date_convention"] = "cone_transition_end"
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_ose(; overwrite = "--force" in ARGS)
end
