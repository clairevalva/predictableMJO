# Each stage runs in a fresh Julia process, freeing its arrays before the next stage.
using HDF5
include(joinpath(@__DIR__, "paper_config.jl"))
include(joinpath(@__DIR__, "validate_inputs.jl"))

function run_paper(; include_large_composites = false, reuse_cache = false, overwrite = false,
    checksums = false)
    validate_inputs(; checksums)
    stages = [("doNLSA.jl", String[], ["init_NLSA.h5"]),
        ("doKoopman.jl", String[], ["init_Keigen.h5"])]
    for variable in ("RMMdata", "RMM", "olr", "u")
        push!(stages, ("doReconstruction.jl", [variable], ["init_recon_$(variable).h5"]))
    end
    push!(stages, ("computeOSE.jl", String[], [basename(OSE_FILE)]))
    push!(stages, (joinpath("model_comparisons", "rmm_prediction_skill.jl"), String[],
        ["rmm_pred_skill_$(PAPER.analysis_start)_$(PAPER.analysis_end).h5"]))
    for (script, args, products) in stages
        outputs = paper_cache.(products)
        if reuse_cache && all(isfile, outputs)
            @info "Reusing existing stage outputs" script
            continue
        end
        overwrite && push!(args, "--force")
        @info "Running stage" script
        run(`$(Base.julia_cmd()) --startup-file=no --project=$REPO_ROOT $(joinpath(REPO_ROOT, script)) $args`)
    end
    args = include_large_composites ? ["--large-composites"] : String[]
    @info "Rendering figures" include_large_composites
    run(`$(Base.julia_cmd()) --startup-file=no --project=$REPO_ROOT $(joinpath(REPO_ROOT, "plotting", "all_paper_plotting.jl")) $args`)
end

if abspath(PROGRAM_FILE) == @__FILE__
    all(a -> a in ("--large-composites", "--reuse-cache", "--force", "--checksums"), ARGS) ||
        error("Usage: run_paper.jl [--large-composites] [--reuse-cache] [--force] [--checksums]")
    run_paper(; include_large_composites = "--large-composites" in ARGS,
        reuse_cache = "--reuse-cache" in ARGS, overwrite = "--force" in ARGS,
        checksums = "--checksums" in ARGS)
end
