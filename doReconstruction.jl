include(joinpath(@__DIR__, "utils", "utils.jl"))
include(joinpath(@__DIR__, "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))

"""Reconstruct only the observables/modes consumed by the public figures.

Input arrays use spatial coordinate × day. Stored reconstructions use
mode × coordinate × day; coefficients use mode × coordinate × delay.
RMMdata retains every nonnegative-frequency mode for the eigenperiod color scale.
The OLR, wind and RMM outputs retain mode 7.
"""
function compute_reconstruction(variable; output = paper_cache("init_recon_$(variable).h5"),
    overwrite = false, settings = PAPER, feature_file = FEATURE_FILE,
    koopman_file = paper_cache("init_Keigen.h5"))
    variable in ("RMMdata", "RMM", "olr", "u") || throw(ArgumentError("Unknown variable $variable"))
    isfile(output) && !overwrite && error("Output exists: $output; pass --force")
    omegas, _, zetaL, zetaR, _ = getK(koopman_file)
    modes = variable == "RMMdata" ? findall(imag.(omegas) .>= 0) : [7]
    dates = settings.start_date .+ Day.((settings.train_first:settings.train_last) .- 1)
    if variable == "RMM"
        rmm1, rmm2, _, _ = getRMM(first(dates), last(dates))
        X = permutedims(hcat(rmm1, rmm2))
    else
        if variable == "RMMdata"
            field = h5read(feature_file, "olr_u")
        else
            _, _, path, key, _, multiplier = data_details(variable)
            field = get_data(path, key, multiplier)
        end
        X = reshape(field, :, size(field, ndims(field)))[:, settings.train_first:settings.train_last]
    end
    nT, nD = size(zetaL, 1), size(X, 1)
    recons = zeros(ComplexF64, length(modes), nD, nT + settings.nemb - 1)
    coeffs = zeros(ComplexF64, length(modes), nD, settings.nemb)
    for (i, k) in enumerate(modes)
        @info "Reconstructing" variable mode=k
        A = ConsistentKoopman.projectdata(X, zetaR[:, k], 1 / nT, settings.nemb)
        recons[i, :, :] = ConsistentKoopman.reconstructfromproj(A, zetaL[:, k], settings.nemb)
        coeffs[i, :, :] = A
    end
    return write_cache(output; overwrite) do f
        f["recons"] = recons
        f["indices"] = modes
        f["Acoeff"] = coeffs
        f["datasize"] = collect(size(X))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    variables = filter(!=("--force"), ARGS)
    isempty(variables) && (variables = ["RMMdata", "RMM", "olr", "u"])
    for variable in variables
        compute_reconstruction(variable; overwrite = "--force" in ARGS)
        GC.gc()
    end
end
