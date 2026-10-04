# checked functionality, (lack of) info leak on April 16, 2026.
# functions to implement analog forecasting.
# write functions, assuming kernel prediction is okay.

"""
    compute_expansion_coeffs_j(phi_j, target_variable, lag, nT = nothing, nD = nothing)

    Computes the expansion coefficients of a target variable for a given lag (prediction) for a single NLSA eigenvalue.
    Based on psuedocode in: 10.1016/j.physd.2020.132520

    # Arguments
    - `phi_j`: the j-th NLSA eigenfunction, of shape (nT, )
    - `target_variable`: the variable we want to forecast, of shape (nT, nD)
    - `lag`: the lag (in time steps) for which we want to compute the expansion coefficients. For example, if lag = 1, we are computing the expansion coefficients for the target variable at time t+1, using the NLSA eigenfunction values at time t.
    - `nT`: the number of time steps in the target variable. If not provided, it will be inferred from the shape of `target_variable`.
    - `nD`: the number of dimensions in the target variable. If not provided, it will be inferred from the shape of `target_variable`.

    # Returns
    - `coeffs`: the expansion coefficients for the target variable at the given lag, of shape (nD, )

"""

function compute_expansion_coeffs_j(phi_j, target_variable, lag, nT = nothing, nD = nothing)
    if (nT == nothing) || (nD == nothing)
        nT, nD = size(target_variable)
    end

    coeffs = zeros(nD) * 1im
    for d = 1:nD
        coeffs[d] = sum(phi_j[1:(end - lag)] .* target_variable[(1 + lag):end, d]) / (nT - lag)
    end
    
    return coeffs
end

"""
    compute_expansion_coeffs(phis, target_variable, maxlag)

    Computes the expansion coefficients of a target variable for a given set of lags (prediction) for all NLSA eigenfunctions.
    Based on psuedocode in: 10.1016/j.physd.2020.132520

    # Arguments
    - `phis`: the NLSA eigenfunctions, of shape (nT, nPhi)
    - `target_variable`: the variable we want to forecast, of shape (nT, nD)
    - `maxlag`: the maximum lag (in time steps) for which we want to compute the expansion coefficients. For example, if maxlag = 10, we are computing the expansion coefficients for the target variable at time t+1, t+2, ..., t+10, using the NLSA eigenfunction values at time t.

    # Returns
    - `expansion_coeffs`: the expansion coefficients for the target variable at the given lags, of shape (nPhi, nD, maxlag + 1). The first dimension corresponds to the NLSA eigenvalues, the second dimension corresponds to the dimensions of the target variable, and the third dimension corresponds to the lags (from 0 to maxlag).
"""
function compute_expansion_coeffs(phis, target_variable, maxlag)
    nPhi = size(phis, 2)
    nT, nD = size(target_variable)

    expansion_coeffs = zeros(nPhi, nD, maxlag + 1) * 1im
    for k = 1:nPhi
        phi_k = phis[:, k]
        for lag = 0:maxlag
            expansion_coeffs[k, :, lag + 1] = compute_expansion_coeffs_j(phi_k, target_variable, lag, nT, nD)
        end
    end

    return expansion_coeffs

end

"""
    make_predictions_from_coeffs(phis_ose, expansion_coeffs)

    From expansion coeffs, make predictions for the target variable at the out of sample times. 
    Based on psuedocode in: 10.1016/j.physd.2020.132520

    # Arguments
    - 'phis_ose': the NLSA eigenfunctions evaluated at the out of sample times, of shape (nT_ose, nPhi)
    - 'expansion_coeffs': the expansion coefficients for the target variable at the given lags, of shape (nPhi, nD, maxlag + 1)
    
    # Returns
    - 'forecasted_variable': the forecasted variable at the out of sample times for the given lags, of shape (nD, nT_ose, maxlag + 1). The first dimension corresponds to the dimensions of the target variable, the second dimension corresponds to the out of sample time steps, and the third dimension corresponds to the lags (from 0 to maxlag).
"""
function make_predictions_from_coeffs(phis_ose, expansion_coeffs)
    nPhi, nD, maxlag_plus1 = size(expansion_coeffs)
    nT_ose = size(phis_ose, 1)

    forecasted_variable = zeros(nD, nT_ose, maxlag_plus1) * 1im
    for lag = 0:(maxlag_plus1 - 1)
        for d = 1:nD
            for t = 1:nT_ose
                phi_t = phis_ose[t, :]
                forecasted_variable[d, t, lag + 1] = sum(expansion_coeffs[:, d, lag + 1] .* phi_t)
            end
        end
    end

    return forecasted_variable
end

"""
    make_predictions(phis_orig, phis_ose, target_variable, max_lag::Int, ret_exp::Bool = false)

    A wrapper function that computes the expansion coefficients and makes predictions for the target variable at the out of sample times for the given lags. Assumes that NLSA kernel evaluations have already been computed for the out of sample times.

    # Arguments
    - `phis_orig`: the NLSA eigenfunctions evaluated at the original (in sample) times, of shape (nT, nPhi)
    - `phis_ose`: the NLSA eigenfunctions evaluated at the out of sample times, of shape (nT_ose, nPhi)
    - `target_variable`: the variable we want to forecast, of shape (nT, nD)
    - `max_lag`: the maximum lag (in time steps) for which we want to make predictions. For example, if max_lag = 10, we are making predictions for the target variable at time t+1, t+2, ..., t+10, using the NLSA eigenfunction values at time t.
    - `ret_exp`: a boolean flag indicating whether to return the expansion coefficients along with the forecasted variable. If true, the function will return a tuple of (forecasted_variable, expansion_coeffs). If false, the function will return only the forecasted_variable.

    # Returns
    - `forecasted_variable`: the forecasted variable at the out of sample times for the given lags, of shape (nD, nT_ose, max_lag + 1). The first dimension corresponds to the dimensions of the target variable, the second dimension corresponds to the out of sample time steps, and the third dimension corresponds to the lags (from 0 to max_lag).
    - `expansion_coeffs` (optional): the expansion coefficients for the target variable at the given lags, of shape (nPhi, nD, max_lag + 1). The first dimension corresponds to the NLSA eigenvalues, the second dimension corresponds to the dimensions of the target variable, and the third dimension corresponds to the lags (from 0 to max_lag).

"""
function make_predictions(phis_orig, phis_ose, target_variable, max_lag, ret_exp = false)
    expansion_coeffs = compute_expansion_coeffs(phis_orig, target_variable, max_lag)
    forecasted_variable = make_predictions_from_coeffs(phis_ose, expansion_coeffs)

    if ret_exp
        return forecasted_variable, expansion_coeffs
    else
        return forecasted_variable
    end
end

"""
    make_direct_issue_date_forecast(phis_orig, phis_ose, target_by_time, nEmb, maxlag)

Directly forecast a raw physical target on the issue-date convention.

`target_by_time` is the uncropped raw target with shape `(nT_raw, nD)`. Row `t`
of `phis_orig` is paired with the newest/current coordinate of the corresponding
delay window:

```text
target_by_time[t + nEmb - 1, :]
```

The output has shape `(nD, nT_ose, maxlag + 1)`, and
`forecasted_variable[:, i, lag + 1]` targets the raw variable at OSE issue row
`i` plus `lag` days.
"""
function make_direct_issue_date_forecast(phis_orig, phis_ose, target_by_time, nEmb, maxlag)
    nT_emb = size(phis_orig, 1)
    target_start = nEmb
    target_end = nT_emb + nEmb - 1
    @assert size(target_by_time, 1) >= target_end

    target_variable = target_by_time[target_start:target_end, :]
    return make_predictions(phis_orig, phis_ose, target_variable, maxlag, false)
end

"""
    reconstruct_forecast_from_delay_coeffs(
        Acoeff, indices_rec, use_k, forecasted_zetas, nEmb, maxlag, nD;
        zeta_indices = nothing,
    )

Reconstruct a physical forecast from forecasted Koopman eigenfunctions and
delay-reconstruction coefficients without borrowing forecasts issued at later
times in the delay window.

The reconstruction coefficients use `q = 1` for the oldest coordinate in the
delay window and `q = nEmb` for the newest coordinate. The paper uses the full
diagonal average in `reconstructfromproj`. That makes RMM lead `lag` depend
on zeta leads `lag:(lag + nEmb - 1)` issued from the same OSE row; this is
forecast-clean when those zeta leads are predicted rather than read from later
OSE rows, but lag 0 depends on predicted future zeta leads.
"""
function reconstruct_forecast_from_delay_coeffs(
    Acoeff,
    indices_rec,
    use_k,
    forecasted_zetas,
    nEmb,
    maxlag,
    nD;
    zeta_indices = nothing,
)
    n_zetas, nT_forecast, n_zeta_lags = size(forecasted_zetas)

    if isnothing(zeta_indices)
        zeta_indices = collect(1:n_zetas)
    end
    @assert length(zeta_indices) == n_zetas

    max_zeta_lag = required_zeta_maxlag(maxlag, nEmb)
    @assert n_zeta_lags >= max_zeta_lag + 1

    coeff_locs_maybe = [findfirst(==(j_val), indices_rec) for j_val in use_k]
    zeta_locs_maybe = [findfirst(==(j_val), zeta_indices) for j_val in use_k]
    @assert all(j_loc -> !isnothing(j_loc), coeff_locs_maybe)
    @assert all(j_loc -> !isnothing(j_loc), zeta_locs_maybe)

    coeff_locs = Int.(coeff_locs_maybe)
    zeta_locs = Int.(zeta_locs_maybe)

    forecasted_variable = zeros(nD, nT_forecast, maxlag + 1)
    for lag = 0:maxlag
        for t = 1:nT_forecast
            rec_sum = zeros(nD) * 1im
            for i in 1:length(use_k)
                coeff_j = Acoeff[coeff_locs[i], :, :]
                zeta_loc = zeta_locs[i]

                for q in 1:nEmb
                    zeta_lag = lag + nEmb - q
                    rec_sum .+= coeff_j[:, q] .* forecasted_zetas[zeta_loc, t, zeta_lag + 1]
                end
            end

            forecasted_variable[:, t, lag + 1] = 2 * real.(rec_sum ./ nEmb)
        end
    end

    return forecasted_variable
end

"""Largest zeta lead needed for the paper's full-delay RMM reconstruction."""
required_zeta_maxlag(maxlag, nEmb) = maxlag + nEmb - 1
