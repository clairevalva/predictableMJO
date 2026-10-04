using TOML, SHA, HDF5, Dates
if !isdefined(@__MODULE__, :data_details)
    include(joinpath(@__DIR__, "utils", "utils.jl"))
end

# Snapshot identity is optional; the prepared-data contract is always checked.
function check_input_snapshot(path, record; checksums = false)
    if checksums
        filesize(path) == record["bytes"] || error("Wrong snapshot size: $path")
        bytes2hex(open(SHA.sha256, path)) == record["sha256"] || error("Wrong snapshot checksum: $path")
    end
    return true
end

function input_units(dataset)
    haskey(HDF5.attributes(dataset), "units") || error("Missing units on $(HDF5.name(dataset))")
    return String(read(HDF5.attributes(dataset)["units"]))
end

# Common CF spellings, including superscripts and slash notation, have the same meaning.
function check_input_units(dataset, allowed)
    units = input_units(dataset)
    normalized = lowercase(replace(strip(units), "**" => "", "^" => "", "⁻" => "-",
        "²" => "2", "¹" => "1", "_" => "", r"\s+" => ""))
    normalized in allowed || error("Incompatible units '$units' on $(HDF5.name(dataset))")
    return nothing
end

function validate_climate_input(path, variable)
    _, levels, _, key, pressures, _ = data_details(variable)
    h5open(path, "r") do f
        required = [key, "valid_time", "longitude", "latitude"]
        isnothing(pressures) || push!(required, "pressure_level")
        all(k -> haskey(f, k), required) || error("Missing climate variable or coordinates in $path")
        shape = levels == 1 ? (240, 121, 9132) : (240, 121, levels, 9132)
        size(f[key]) == shape || error("Unexpected $variable dimensions in $path")
        read(f, "valid_time") == collect(0:9131) || error("Expected daily fields from 2000-01-01 through 2024-12-31")
        parse_units_datetime(input_units(f["valid_time"])) == DateTime(PAPER.start_date) ||
            error("Expected climate time origin 2000-01-01 00:00:00")
        read(f, "longitude") ≈ collect(0:1.5:358.5) || error("Unexpected longitude grid")
        read(f, "latitude") ≈ collect(-90:1.5:90) || error("Unexpected latitude grid")
        check_input_units(f["longitude"], ("degreeseast", "degreeeast", "degreese", "degreee"))
        check_input_units(f["latitude"], ("degreesnorth", "degreenorth", "degreesn", "degreen"))
        if !isnothing(pressures)
            read(f, "pressure_level") == pressures || error("Expected pressure levels [200, 850] hPa")
            check_input_units(f["pressure_level"], ("hpa", "mbar", "millibar", "millibars"))
        end
        units = variable == "olr" ? ("jm-2", "j/m2") : variable == "u" ? ("ms-1", "m/s") :
            ("m", "meter", "meters", "metre", "metres")
        check_input_units(f[key], units)
    end
    return true
end

function validate_model_input(path)
    h5open(path, "r") do f
        for key in ("RMM1", "RMM2", "valid_time", "step")
            haskey(f, key) || error("Missing forecast variable $key in $path")
        end
        shape = size(f["RMM1"])
        length(shape) == 2 && shape == size(f["RMM2"]) == size(f["valid_time"]) ||
            error("Expected matching lead × issue-date forecast arrays in $path")
        steps = vec(read(f, "step"))
        length(steps) == shape[1] && all(x -> isfinite(x) && x >= 0, steps) ||
            error("Invalid forecast lead grid in $path")
        if haskey(HDF5.attributes(f["step"]), "units")
            check_input_units(f["step"], ("d", "day", "days"))
        end
        dates = Date.(convert_valid_time(read(f, "valid_time"), input_units(f["valid_time"])))
        any(d -> PAPER.analysis_start <= d <= PAPER.analysis_end, dates) ||
            error("Forecast dates do not overlap the analysis period in $path")
    end
    return true
end

"""Validate prepared inputs; optionally also require the recorded snapshot sizes and SHA-256 hashes."""
function validate_inputs(; checksums = false)
    records = TOML.parsefile(joinpath(REPO_ROOT, "input_manifest.toml"))["files"]
    require_files([paper_input(split(r["path"], '/')...) for r in records])
    for record in records
        check_input_snapshot(paper_input(split(record["path"], '/')...), record; checksums)
    end
    for variable in ("olr", "u", "precip")
        validate_climate_input(data_details(variable)[3], variable)
    end
    h5open(FEATURE_FILE, "r") do f
        haskey(f, "olr_u") && size(f["olr_u"]) == (3, 240, 9132) ||
            error("Expected variable × longitude × day features")
        for key in ("stdu", "stdOLR", "meanu", "meanOLR")
            haskey(f, key) || error("Missing feature normalization metadata: $key")
            values = read(f, key)
            all(isfinite, values) || error("Invalid normalization metadata: $key")
            if startswith(key, "std")
                all(>(0), values) || error("Normalization scale must be positive: $key")
            end
        end
    end
    for values in getRMM(PAPER.start_date, PAPER.ose_end)[1:2]
        all(x -> isfinite(x) && abs(x) < 999, values) || error("Missing RMM observations")
    end
    for values in getBSISOkikuchi(PAPER.start_date + Day(1), PAPER.start_date + Day(PAPER.train_last - 1))[1:2]
        all(isfinite, values) || error("Missing BSISO observations")
    end
    yotc_dates = collect(Date(2009, 10, 1):Day(1):Date(2010, 4, 1))
    for index in (:omi, :fmo)
        _load_noaa_mjo(index, yotc_dates)
    end
    for model in ("ecmwf", "hmcr", "kma", "ncep", "ukmo")
        validate_model_input(paper_input("model_RMM_pred", "merged_RMM", "RMMS_$(model)_realtime.nc"))
    end
    println("Validated ", length(records), " prepared input files", checksums ? " including snapshot checksums" : " (schema and dates)")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    all(==("--checksums"), ARGS) || error("Usage: validate_inputs.jl [--checksums]")
    validate_inputs(; checksums = "--checksums" in ARGS)
end
