using Statistics

function stack_util(complex_vector::AbstractVector)
    # stack real and imaginary parts of a complex vector into a 2 x N matrix
    real_part = real.(complex_vector)
    imag_part = imag.(complex_vector)
    stacked = transpose(hcat(real_part, imag_part))

    return stacked
end



# also introduct simple metrics from the review paper Kim et al. 2018 for RMM
function bivariate_correlation(ζpred::AbstractMatrix, ζtrue::AbstractMatrix; time_dim::Int=2, centered::Bool = true)
    # compute bivariate correlation 2d time series ζpred and ζtrue. Only accepts time_dim = 2 right now. Should match just taking the covariance between the two complex time series.
    @assert size(ζpred) == size(ζtrue) "Shapes must match"
    @assert time_dim == 2 "Only time_dim = 2 is implemented"
    @assert centered "Currently only centered=true is implemented"

    μ_pred = mean(ζpred, dims=time_dim)
    μ_true = mean(ζtrue, dims=time_dim)
    N = size(ζpred, time_dim)

    ap = ζpred .- μ_pred
    at = ζtrue .- μ_true

    
    cov = mean(ap .* conj.(at), dims=time_dim)
    varp  = mean(abs2.(ap), dims=time_dim)
    vart  = mean(abs2.(at), dims=time_dim)

    frac_top = N * sum(cov) 
    frac_bot = sqrt.(N * sum(varp) * N * sum(vart))
    
    return frac_top ./ frac_bot
end


using Dates
using CSV

function parse_units_datetime(units::AbstractString)
    match_units = match(r"^days\s+since\s+(\d{4}-\d{2}-\d{2})(?:[ T](\d{2}:\d{2}:\d{2})(?:\.0+)?)?(?:\s*(?:UTC|Z))?$"i, strip(units))
    isnothing(match_units) && error("Expected daily time units, got: $units")
    time = something(match_units.captures[2], "00:00:00")
    return DateTime(match_units.captures[1] * " " * time, dateformat"yyyy-mm-dd HH:MM:SS")
end

function convert_valid_time(valid_time_read, units::AbstractString)
    base_dt = parse_units_datetime(units)
    millis = round.(Int, valid_time_read .* 86_400_000)
    return base_dt .+ Millisecond.(millis)
end

function load_rmm_dates(startdate::Date, enddate::Date;
    fp = paper_input("climdata", "RMM_data_edit.csv"))
    rmmfile = CSV.File(fp)
    filedates = Date.(rmmfile["year"], rmmfile["month"], rmmfile["day"])
    datebool = (filedates .>= startdate) .& (filedates .<= enddate)
    return filedates[datebool]
end

function seasonal_bivar(pred_mat::AbstractMatrix, obs_mat::AbstractMatrix, dates::AbstractVector{Date}, seasons_js, min_samples::Int)
    n_seasons = size(seasons_js, 1) + 1
    corr = Vector{Float64}(undef, n_seasons)
    counts = Vector{Int}(undef, n_seasons)

    seasons_bool = makeseasons(dates, seasons_js)
    for seas_j in 1:n_seasons
        if seas_j == 1
            seas_mask = trues(length(dates))
        else
            seas_mask = seasons_bool[seas_j - 1, :]
        end
        n_samples = count(seas_mask)
        counts[seas_j] = n_samples
        if n_samples < min_samples
            corr[seas_j] = NaN
        else
            corr[seas_j] = bivariate_correlation(pred_mat[:, seas_mask], obs_mat[:, seas_mask])
        end
    end

    return corr, counts
end
