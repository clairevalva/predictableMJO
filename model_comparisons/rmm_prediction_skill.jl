# Bivariate RMM forecast-skill curves for the paper figure.
#
# Merges the four curve families that used to be split across all_prediction_skill.jl and
# all_prediction_skill_direct_RMM.jl (S2S models, the zeta forecast, the KAF of RMM, and RMM
# reconstructed from forecast zeta), and caches them so plotting/rmmSkillSummary.jl can
# render without recomputing. Run this once; rerun only when the OSE or Koopman data change.

using Dates
using HDF5
using Statistics

include(joinpath(@__DIR__, "..", "utils", "utils.jl"))
include(joinpath(@__DIR__, "..", "analogForecasting.jl"))
if !isdefined(@__MODULE__, :ConsistentKoopman)
    include(joinpath(@__DIR__, "..", "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))
end

# Editable settings.
const MODE_SETS = ["j7" => [7]]
const ZETA_EIG_INDEX = 7
const MODELS_USE = ["ecmwf", "hmcr", "kma", "ncep", "ukmo"]
const ANALYSIS_START = PAPER.analysis_start
const ANALYSIS_END = PAPER.analysis_end
const MAXLAG = PAPER.maxlag
const MIN_SAMPLES = 2
const SEASONS = ["all", "DJFMA", "JJASO"]

# Row 1 is DJFMA, row 2 JJASO; "all" is the implicit first slot inside seasonal_bivar.
const SEASONS_JS = let js = zeros(2, 5)
    js[1, :] = [12, 1, 2, 3, 4]
    js[2, :] = [6, 7, 8, 9, 10]
    js
end

const EXPNAME = PAPER.experiment
const OSE_START, OSE_END = PAPER.ose_start, PAPER.ose_end
const SN_P_EVALS = OSE_FILE
const SN_NLSA = paper_cache( "$(EXPNAME)_NLSA.h5")
const SN_K = paper_cache( "$(EXPNAME)_Keigen.h5")
const FN_REC = paper_cache( "init_recon_RMM.h5")
const MODEL_DIR = paper_input("model_RMM_pred", "merged_RMM")

skill_cache_path() = paper_cache(
    "rmm_pred_skill_$(ANALYSIS_START)_$(ANALYSIS_END).h5")

# The KOSE file already stores the OSE phi evaluations; recomputing them (101 modes x ~1800
# times through evalPhi) was the dominant cost of the old scripts. Spot-check and reuse.
function load_phis_ose(orig_phis, lambdas)
    phis_ose = permutedims(h5read(SN_P_EVALS, "phi_vals"))
    ose_p_evals = h5read(SN_P_EVALS, "p_evals")
    maxdiff = 0.0
    for j in (1, 500, size(phis_ose, 1)), k in (1, 7, 50, 101)
        v = ConsistentKoopman.evalPhi(ose_p_evals[j, :], orig_phis[:, k], lambdas[k])
        maxdiff = max(maxdiff, abs(v - phis_ose[j, k]))
    end
    println("cached phi_vals spot-check max abs diff: ", maxdiff)
    @assert maxdiff < 1e-8
    return phis_ose
end

# Bivariate correlation of a (nD x nT x nlag) forecast cube against (nD x nT) truth, per lag
# and per season, restricted to the analysis window. Returns an (nlag+1) x nseasons matrix.
function skill_by_lag(forecast, truth, dates)
    nT = minimum((size(forecast, 2), size(truth, 2), length(dates)))
    forecast = forecast[:, 1:nT, :]
    truth = truth[:, 1:nT]
    dates = dates[1:nT]

    corr = fill(NaN, MAXLAG + 1, length(SEASONS))
    for lag in 0:MAXLAG
        pred = forecast[:, 1:(end - lag), lag + 1]
        obs = truth[:, (1 + lag):end]
        valid_dates = dates[(1 + lag):end]

        window = (valid_dates .>= ANALYSIS_START) .& (valid_dates .<= ANALYSIS_END)
        any(window) || continue

        seasonal, _ = seasonal_bivar(pred[:, window], obs[:, window], valid_dates[window],
            SEASONS_JS, MIN_SAMPLES)
        corr[lag + 1, :] = seasonal
    end
    return corr
end

# S2S models come with their own lead times (each model's `step`), not the Koopman 0:MAXLAG grid.
function compute_model_skill(model_path, rmm1_obs, rmm2_obs, obs_date_to_idx)
    rmm1 = read(model_path, "RMM1")
    rmm2 = read(model_path, "RMM2")
    step = read(model_path, "step")
    valid_time = convert_valid_time(read(model_path, "valid_time"),
        read(HDF5.attributes(model_path["valid_time"])["units"]))

    lags = Float64[]
    rows = Vector{Vector{Float64}}()
    for s in 1:size(valid_time, 1)
        dates = Date.(valid_time[s, :])
        pred1 = rmm1[s, :]
        pred2 = rmm2[s, :]

        window = (dates .>= ANALYSIS_START) .& (dates .<= ANALYSIS_END)
        any(window) || continue
        dates, pred1, pred2 = dates[window], pred1[window], pred2[window]

        obs_idx = [get(obs_date_to_idx, d, 0) for d in dates]
        keep = obs_idx .> 0
        any(keep) || continue
        obs_idx, dates, pred1, pred2 = obs_idx[keep], dates[keep], pred1[keep], pred2[keep]

        pred = transpose(hcat(pred1, pred2))
        obs = transpose(hcat(rmm1_obs[obs_idx], rmm2_obs[obs_idx]))
        corr, _ = seasonal_bivar(pred, obs, dates, SEASONS_JS, MIN_SAMPLES)

        push!(lags, float(step[s]))
        push!(rows, collect(corr))
    end

    order = sortperm(lags)
    corr = Matrix{Float64}(undef, length(order), length(SEASONS))
    for (i, o) in enumerate(order)
        corr[i, :] = rows[o]
    end
    return lags[order], corr
end

# Forecast of the Koopman index itself. NOTE: as in the original script, the truth is the lag-0
# KAF estimate of zeta, not the directly evaluated OSE zeta (zetaL_vals) -- a self-prediction score.
function compute_zeta_skill(orig_phis, phis_ose, zeta_L, ose_dates)
    forecast, _ = make_predictions(orig_phis[:, 1:101], phis_ose, zeta_L, MAXLAG, true)
    ZETA_EIG_INDEX <= size(forecast, 1) ||
        error("ZETA_EIG_INDEX=$(ZETA_EIG_INDEX) exceeds available modes $(size(forecast, 1))")

    corr = fill(NaN, MAXLAG + 1, length(SEASONS))
    for lag in 0:MAXLAG
        dates = ose_dates[(1 + lag):end]
        window = (dates .>= ANALYSIS_START) .& (dates .<= ANALYSIS_END)
        any(window) || continue

        pred = stack_util(forecast[ZETA_EIG_INDEX, 1:(end - lag), lag + 1][window])
        obs = stack_util(forecast[ZETA_EIG_INDEX, (1 + lag):end, 1][window])
        seasonal, _ = seasonal_bivar(pred, obs, dates[window], SEASONS_JS, MIN_SAMPLES)
        corr[lag + 1, :] = seasonal
    end
    return collect(0:MAXLAG), corr
end

# RMM reconstructed from forecast zeta, for a given set of Koopman eigenfunctions.
function compute_koopman_rmm_skill(use_k, orig_phis, phis_ose, zeta_L, ose_dates, nEmb)
    maxlag_zeta = required_zeta_maxlag(MAXLAG, nEmb)
    _, indices_rec, Acoeff = get_reconstructed_data(FN_REC; ret_coeffs = true,
        load_small_j = use_k)

    forecasted_zeta, _ = make_predictions(orig_phis[:, 1:101], phis_ose, zeta_L[:, use_k],
        maxlag_zeta, true)
    forecast = reconstruct_forecast_from_delay_coeffs(Acoeff, indices_rec, use_k,
        forecasted_zeta, nEmb, MAXLAG, 2; zeta_indices = use_k)

    rmm1, rmm2, _, _ = getRMM(ose_dates[1], ose_dates[end])
    return collect(0:MAXLAG), skill_by_lag(forecast, stack([rmm1, rmm2], dims = 1), ose_dates)
end

# Kernel analog forecast of RMM: regress raw RMM onto the NLSA eigenfunctions directly.
function compute_kaf_rmm_skill(orig_phis, phis_ose, ose_dates, nEmb)
    tStart, tEnd, _, _, data_start = getExperiment(EXPNAME)
    train_dates = data_start + Day.((tStart:tEnd) .- 1)
    rmm1_train, rmm2_train, _, _ = getRMM(train_dates[1], train_dates[end])

    forecast = real.(make_direct_issue_date_forecast(orig_phis[:, 1:101], phis_ose,
        hcat(rmm1_train, rmm2_train), nEmb, MAXLAG))

    rmm1, rmm2, _, _ = getRMM(ose_dates[1], ose_dates[end])
    return collect(0:MAXLAG), skill_by_lag(forecast, stack([rmm1, rmm2], dims = 1), ose_dates)
end

function compute_rmm_skill_cache(; output = skill_cache_path(), overwrite = false)
    isfile(output) && !overwrite && error("Output exists: $output; pass --force")
    _, _, _, nEmb, _ = getExperiment(EXPNAME)
    orig_phis = h5read(SN_NLSA, "phi")
    lambdas = h5read(SN_NLSA, "kappa")
    zeta_L = h5read(SN_K, "zetaL")
    ose_dates = get_OSE_dates(SN_P_EVALS; fallback_start = OSE_START, fallback_end = OSE_END)
    phis_ose = load_phis_ose(orig_phis, lambdas)
    @assert size(phis_ose, 1) == length(ose_dates)

    curves = Pair{String, Tuple{Vector{Float64}, Matrix{Float64}}}[]

    println("S2S models...")
    rmm1_obs, rmm2_obs, _, _ = getRMM(ANALYSIS_START, ANALYSIS_END)
    obs_dates = load_rmm_dates(ANALYSIS_START, ANALYSIS_END)
    obs_date_to_idx = Dict{Date, Int}(d => i for (i, d) in enumerate(obs_dates))
    for model in MODELS_USE
        path = joinpath(MODEL_DIR, "RMMS_$(model)_realtime.nc")
        isfile(path) || error("Missing S2S file for model $(model): $(path)")
        lags, corr = h5open(f -> compute_model_skill(f, rmm1_obs, rmm2_obs, obs_date_to_idx),
            path, "r")
        push!(curves, "s2s_$(model)" => (lags, corr))
    end

    println("Koopman zeta forecast (j=$(ZETA_EIG_INDEX))...")
    lags, corr = compute_zeta_skill(orig_phis, phis_ose, zeta_L, ose_dates)
    push!(curves, "koopman_zeta_j$(ZETA_EIG_INDEX)" => (float.(lags), corr))

    println("KAF of RMM...")
    lags, corr = compute_kaf_rmm_skill(orig_phis, phis_ose, ose_dates, nEmb)
    push!(curves, "kaf_rmm" => (float.(lags), corr))

    for (tag, use_k) in MODE_SETS
        println("RMM from Koopman index, modes $(use_k)...")
        lags, corr = compute_koopman_rmm_skill(use_k, orig_phis, phis_ose, zeta_L, ose_dates,
            nEmb)
        push!(curves, "koopman_rmm_$(tag)" => (float.(lags), corr))
    end

    path = output
    write_cache(path; overwrite) do f
        f["meta/seasons"] = SEASONS
        f["meta/curves"] = [first(c) for c in curves]
        f["meta/analysis_start"] = string(ANALYSIS_START)
        f["meta/analysis_end"] = string(ANALYSIS_END)
        for (key, (lags, corr)) in curves
            f["skill/$(key)/lags"] = lags
            f["skill/$(key)/corr"] = corr
        end
    end
    println("wrote ", path)
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_rmm_skill_cache(; overwrite = "--force" in ARGS)
end
