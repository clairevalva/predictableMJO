function makePhases(eigenf; nP::Int = 8, conventional = true)
    halfphase = ceil(Int, nP / 2)
    nT = length(eigenf)
    ang = mod.(angle.(eigenf), 2*pi) .- pi
    phase_defs = (-halfphase:halfphase)*(pi/halfphase);

    phase_TFs = zeros(Bool, nP, nT);

    for j = 1:nP
        if conventional
            refj = (j + 4) % 8
            if refj === 0
                refj = 8
            end
        else 
            refj = j 
        end
        phase_TFs[refj,:] = (ang .>= phase_defs[j]) .& (ang .< phase_defs[j + 1]);
    end

    return phase_TFs
end

function align_index(eof, ζeof)

    # write function to align eof
    corr_eof = [cor(eof, ζeof), cor(conj(eof), ζeof)]
    idxcorr = argmax(abs.(corr_eof)) 
    maxcorr = abs(corr_eof[idxcorr])
    
    if idxcorr == 2
        println("WARNING: probably shouldn't have this be the case, conjugate eigenfunction?")
        eof = conj.(eof)
    end
    rotation_coeff = corr_eof[idxcorr] / maxcorr
    rotation_ang = angle(rotation_coeff)

    rot_eig = rotation_coeff * ζeof
    newcorr = cor(eof, rot_eig)

    println("max corr (abs): ", maxcorr, ", new corr: ", newcorr)
    println("rot coeff, ", rotation_coeff, ", abs of ", abs.(rotation_coeff))

    return rotation_coeff * ζeof, rotation_ang
end

function makeseason(timearray::AbstractVector, seasmonths)
    nD = length(timearray)
    seasbool = zeros(Bool, nD)

    for j = 1:nD
        seasbool[j] = month(timearray[j]) ∈ seasmonths
    end

    return seasbool
end

function makeseasons(timearray::AbstractVector, seasons = nothing; makelabels = false)
    if seasons === nothing
        seasons = zeros(4, 3);
        seasons[1, :] = [12, 1, 2]
        seasons[2, :] = [3, 4, 5]
        seasons[3, :] = [6, 7, 8]
        seasons[4, :] = [9, 10, 11]
    end
    nS = size(seasons, 1)
    nD = length(timearray)

    seasons_bool = zeros(Bool, nS, nD)

    for j = 1:nS
        seasons_bool[j, :] = makeseason(timearray, seasons[j, :])
    end

    if makelabels
        seasonlabels = zeros(Int, nD)
        for j = 1:nD
            for s = 1:nS
                if seasons_bool[s, j]
                    seasonlabels[j]  = s
                end
            end
        end

        return seasons_bool, seasonlabels
    else
        return seasons_bool
    end

end

function makeseasons(timearray::AbstractVector, seasons::Vector{Vector{Int64}}; makelabels = false)
    nS = size(seasons, 1)
    nD = length(timearray)

    seasons_bool = zeros(Bool, nS, nD)

    for j = 1:nS
        seasons_bool[j, :] = makeseason(timearray, seasons[j])
    end

    if makelabels
        seasonlabels = zeros(Int, nD)
        for j = 1:nD
            for s = 1:nS
                if seasons_bool[s, j]
                    seasonlabels[j]  = s
                end
            end
        end

        return seasons_bool, seasonlabels
    else
        return seasons_bool
    end

end



# Reference NOAA snapshots retrieved 2026-09-08; supplied externally:
# https://psl.noaa.gov/mjo/mjoindex/omi.1x.txt — CPC blended OLR, 1991-01-01–2026-06-24
# SHA256: 52b5d504eb5933b638ace34fde813fc4be5194f1a05156a0879e536eb13d1b49
# https://psl.noaa.gov/mjo/mjoindex/fmo.1x.txt — 1979-01-01–2022-12-31
# SHA256: de91cab8035e7920a4929867491f3b565da5af31f3fcaaf952564b6f81ff9e80
# Both use (PC2, -PC1) for RMM-compatible geographic phases, preserving amplitude.
# OMI: https://psl.noaa.gov/mjo/mjoindex/
# FMO: https://doi.org/10.1175/JCLI-D-14-00688.1 (section 2a)
# Research citation: Kiladis et al. (2014), https://doi.org/10.1175/MWR-D-13-00301.1

# OMI has year/month/day/PC1/PC2/amplitude; FMO inserts hour after day.
function _read_noaa_mjo(path, index::Symbol)
    index in (:omi, :fmo) || throw(ArgumentError("Expected :omi or :fmo"))
    ncols = index == :omi ? 6 : 7
    dates = Date[]
    pc1, pc2, amplitude = Float64[], Float64[], Float64[]
    seen = Set{Date}()
    for (line_number, line) in enumerate(eachline(path))
        isempty(strip(line)) && continue
        fields = split(line)
        length(fields) == ncols || error("$(path):$(line_number): expected $(ncols) columns")
        date = Date(parse.(Int, fields[1:3])...)
        date in seen && error("$(path): duplicate date $(date)")
        isempty(dates) || date > last(dates) || error("$(path): dates must increase")
        if ncols == 7
            parse(Int, fields[4]) == 0 || error("$(path): expected daily observations at hour 0")
        end
        p1, p2, amp = parse.(Float64, fields[end-2:end])
        all(isfinite, (p1, p2, amp)) && amp >= 0 || error("$(path): invalid values on $(date)")
        # Each of the three values is rounded to five decimal places in NOAA's files.
        isapprox(hypot(p1, p2), amp; atol = 2e-5, rtol = 0) ||
            error("$(path): inconsistent amplitude or missing-value sentinel on $(date)")
        push!(seen, date)
        push!(dates, date)
        push!(pc1, p1)
        push!(pc2, p2)
        push!(amplitude, amp)
    end
    isempty(dates) && error("$(path): no observations")
    return (; dates, pc1, pc2, amplitude)
end

function _select_noaa_mjo(data, dates_plot)
    by_date = Dict(date => i for (i, date) in enumerate(data.dates))
    missing_dates = filter(d -> !haskey(by_date, d), dates_plot)
    isempty(missing_dates) || error("NOAA index is missing required dates: $(missing_dates)")
    inds = [by_date[d] for d in dates_plot]
    return complex.(data.pc2[inds], -data.pc1[inds])
end

# Accept valid supplied observations; never download or replace files implicitly.
function _load_noaa_mjo(index::Symbol, dates_plot; data_dir = paper_input("climdata"))
    index in (:omi, :fmo) || throw(ArgumentError("Expected :omi or :fmo"))
    path = joinpath(data_dir, "$(index).1x.txt")
    require_files([path])
    return _select_noaa_mjo(_read_noaa_mjo(path, index), dates_plot)
end

