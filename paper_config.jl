using Dates

const REPO_ROOT = @__DIR__
const INPUT_ROOT = abspath(get(ENV, "PREDICTABLEMJO_INPUT_ROOT", joinpath(REPO_ROOT, "inputs")))
const CACHE_ROOT = abspath(get(ENV, "PREDICTABLEMJO_CACHE_ROOT", joinpath(REPO_ROOT, "data")))
const DEFAULT_FIGURE_OUTPUT_ROOT = abspath(get(ENV, "PREDICTABLEMJO_OUTPUT_ROOT", joinpath(REPO_ROOT, "figures")))
const PAPER = (experiment = "init", start_date = Date(2000, 1, 1),
    train_first = 2, train_last = 7305, nemb = 64, neighbors = 1500, ose_neighbors = 2000,
    nphi = 101, diffusion_modes = 501, dt = 1.0, regularization = 1e-3,
    ose_start = Date(2020, 1, 1), ose_end = Date(2024, 12, 31),
    analysis_start = Date(2020, 12, 1), analysis_end = Date(2024, 4, 30), maxlag = 50)
paper_input(parts...) = joinpath(INPUT_ROOT, parts...)
paper_cache(parts...) = joinpath(CACHE_ROOT, parts...)
const FEATURE_FILE = paper_input("data_clim_processed", "merid_avgs_predict_6_2025.h5")
const OSE_FILE = paper_cache("init_KOSE_$(PAPER.ose_start)_$(PAPER.ose_end).h5")

function require_files(files)
    missing = filter(!isfile, files)
    isempty(missing) || error("Missing inputs (see DATA.md):\n" * join(missing, "\n"))
end

function write_cache(writer, path; overwrite = false)
    isfile(path) && !overwrite && error("Output already exists: $path; pass --force to replace it")
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path)); close(io)
    try
        h5open(writer, temporary, "w")
        mv(temporary, path; force = overwrite)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

"""Paper experiment; training begins one day late for the cone-kernel transition."""
function getExperiment(expname = PAPER.experiment)
    expname == PAPER.experiment || throw(ArgumentError("Only the paper's init experiment is included"))
    return PAPER.train_first, PAPER.train_last, PAPER.neighbors, PAPER.nemb,
        PAPER.start_date, "None", PAPER.nphi
end
