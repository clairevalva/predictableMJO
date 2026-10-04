using CSV


function data_details(var)
    if var == "olr"
        units_data = "Wm⁻²"
        nL = 1
        pressures = nothing

        data_fn = paper_input("data_clim_processed", "OLR_feb25.nc")
        data_load = "ttr"

        multiplier = -1 / (60*60) # reverse sign and convert from J m^-2 to W m^-2, accumulation period if 1 hour

    elseif var == "u"
        units_data = "ms⁻¹"
        nL = 2
        pressures = [200, 850]

        data_fn = paper_input("data_clim_processed", "u_feb25.nc")
        data_load = "u"
        multiplier = 1
    elseif var == "precip"
        units_data = "mm day⁻¹"
        nL = 1
        pressures = nothing
        data_fn = paper_input("data_clim_processed", "precip_sep25_02.nc")
        data_load = "tp"
        multiplier = 1000 # convert to mm / day
    else
        error("data variable not recognized")
    end

    return units_data, nL, data_fn, data_load, pressures, multiplier
end


function get_data(fn, data_load, multiplier; retall = false)
    X = h5read(fn, data_load)
    X = multiplier * X

    if retall
        lon = h5read(fn, "longitude")
        lat = h5read(fn, "latitude")
        if h5open(f -> haskey(f, "pressure_level"), fn, "r")
            lev = h5read(fn, "pressure_level")
        else
            lev = nothing
        end
        time = h5read(fn, "valid_time")
        return X, lon, lat, lev, time
    end

    return X
end



function getK(koop_fn)
    omegas, c = h5read(koop_fn, "omegas"), h5read(koop_fn, "c")
    zetaL, zetaR = h5read(koop_fn, "zetaL"), h5read(koop_fn, "zetaR")

    omegas = omegas / (2 * pi)
    periods = 1 ./ imag.(omegas)
    
    return omegas, c, zetaL, zetaR, periods
end

function getRMM(startdate, enddate, fp = paper_input("climdata", "RMM_data_edit.csv"))
    RMMfile = CSV.File(fp)

    fileyears = RMMfile["year"]
    filemonths = RMMfile["month"]
    filedays = RMMfile["day"]

    filedates = Date.(fileyears, filemonths, filedays)
    datebool = (filedates .>= startdate) .& (filedates .<= enddate) 
    
    select_dates(filedates[datebool], collect(startdate:Day(1):enddate))
    RMM1 = (RMMfile["RMM1"])[datebool]
    RMM2 = (RMMfile["RMM2"])[datebool]
    phases_RMM = (RMMfile["phase"])[datebool]

    index = RMM1 + 1im*RMM2

    return RMM1, RMM2, index, phases_RMM
end

function getBSISOkikuchi(startdate, enddate; fp = paper_input("climdata", "BSISO_kikuchi.csv"), all_vars = false)
    # Kikuchi reference snapshot; columns and provenance are documented in DATA.md.
    RMMfile = CSV.File(fp)

    filedates = Date.(RMMfile["date"], "mm/dd/yyyy")
    datebool = (filedates .>= startdate) .& (filedates .<= enddate) 
    select_dates(filedates[datebool], collect(startdate:Day(1):enddate))
    PC1 = (RMMfile["PCx"])[datebool]
    PC2 = (RMMfile["PCy"])[datebool]
    index = PC1 + 1im*PC2

    if all_vars
        phase = (RMMfile["phase"])[datebool]
        norm_Amp = (RMMfile["Amp(nrm)"])[datebool]
        nonorm_Amp = (RMMfile["Amp(non-nrm)"])[datebool]
        return PC1, PC2, index, filedates[datebool], phase, norm_Amp, nonorm_Amp
    else
        return PC1, PC2, index
    end

    return PC1, PC2, index
end


function makeTakensInd_forward(nEmb, nT)
    # map embedded index j strictly to physical time t + nEmb - 1 
    return nEmb:(nT + nEmb - 1)
end

function areaweights(lat::AbstractVector, lon::AbstractVector)

    # get area weights
    rad_lat = deg2rad.(lat)
    dlat = rad_lat[1] - rad_lat[2]
    dlon = deg2rad.(lon[2]) - deg2rad.(lon[1])
    w_area = (sin.(rad_lat .+ dlat / 2) .- sin.(rad_lat .- dlat / 2)) / dlon

    return abs.(w_area)
end

function makemeridavg(X, lat_ext, lat, w_area)
    # assumes that the array has dim 2 as latitude
    latj = (abs.(lat) .<= lat_ext)
    nlat = sum(latj)
    X_merid = selectdim(X, 2, latj)
    w_merid = w_area[latj]

    for j = 1:nlat
        selectdim(X_merid, 2, j) .*= w_merid[j] 
    end

    X_avg = mean(X_merid, dims = 2)
    w_avg = mean(w_merid)
    X_avg /= w_avg

    return selectdim(X_avg, 2, 1)
end

function get_OSE_dates(fn; fallback_start = nothing, fallback_end = nothing)
    dates_from_file = h5open(fn, "r") do h5f
        if haskey(h5f, "dates")
            OSE_dates = Date.(read(h5f, "dates"), "yyyy-mm-dd")
            if length(OSE_dates) == 2
                return collect(OSE_dates[1]:Day(1):OSE_dates[end])
            else
                return OSE_dates
            end
        end

        return nothing
    end

    if !isnothing(dates_from_file)
        return dates_from_file
    end

    if isnothing(fallback_start) || isnothing(fallback_end)
        error("OSE date dataset not found in $(fn), and fallback_start/fallback_end were not provided")
    end

    return collect(fallback_start:Day(1):fallback_end)
end

function get_reconstructed_data(fn_rec; ret_coeffs::Bool = false, load_small_j = nothing)
    h5open(fn_rec, "r") do f
        indices = vec(read(f, "indices"))
        positions = isnothing(load_small_j) ? collect(eachindex(indices)) :
            [something(findfirst(==(j), indices), 0) for j in load_small_j]
        all(>(0), positions) || error("Requested reconstruction modes are absent in $fn_rec")
        recs = stack([f["recons"][i, :, :] for i in positions]; dims = 1)
        if ret_coeffs
            A = stack([f["Acoeff"][i, :, :] for i in positions]; dims = 1)
            return recs, indices[positions], A
        end
        return recs, indices[positions]
    end
end

"""Indices of requested physical dates, rejecting missing, duplicate or unsorted input dates."""
function select_dates(available, wanted)
    issorted(available) && allunique(available) || error("Input dates must be unique and increasing")
    index = Dict(d => i for (i, d) in enumerate(available))
    all(d -> haskey(index, d), wanted) || error("Input does not cover all requested dates")
    return [index[d] for d in wanted]
end
