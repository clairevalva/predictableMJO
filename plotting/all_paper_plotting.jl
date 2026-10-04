using CairoMakie, GeoMakie, HDF5, Dates, Statistics, StatsBase, LinearAlgebra
if !isdefined(@__MODULE__, :data_details)
    include(joinpath(@__DIR__, "..", "utils", "utils.jl"))
end
if !isdefined(@__MODULE__, :ConsistentKoopman)
    include(joinpath(REPO_ROOT, "vendor", "ConsistentKoopman", "src", "ConsistentKoopman.jl"))
end
const PAPER_FIGURE_WIDTH = 570
const YOTC_MARK_DATES = [Date(2009, 10, 25), Date(2009, 12, 16)]
const PHASE_REGION_LABELS = (indian_ocean = "Indian Ocean", west_hem_africa = "West. Hem. & Africa",
    west_pacific = "Western Pacific", maritime_continent = "Maritime Continent")
const DEFAULT_PLOT_TITLES = (
    eigenperiod = "Koopman eigenperiod sorting",
    rmm_complex = "MJO index comparison during YOTC period\n(Oct 2009-March 2010)",
    mjo_indices_labeled = "MJO index comparison during YOTC period\n(Oct 2009-March 2010)",
    hov = "Comparison of bandpassed observations and the Koopman reconstruction",
    mjo_composites = "Composites of RMM index and Koopman MJO",
    mjo_precip_dist = "MJO conditional precipitation distribution",
    mjo_composites_precip = "Comparison of RMM and Koopman MJO statistics",
    bsiso_composites = "Composites of BSISO index and Koopman BSISO",
    rmm_skill = "RMM bivariate forecast skill")
include("evPlots.jl")
include(joinpath(REPO_ROOT, "compareRMMindex.jl"))
include(joinpath(REPO_ROOT, "compareHovs.jl"))
include("compositesMJO.jl")
include("precipDistribution.jl")
include("compositesBSISO.jl")
include("rmmSkillSummary.jl")

"""Exact output inventory, relative to output_root. Larger composites are optional."""
function expected_figure_paths(; include_large_composites = false)
    paper = ["periods_comparison_AC_lags=45.png", "periods_comparison_AC_lags=45_labeled.png",
        "RMM_vs_eigenfunction_YOTC_complex_timecolored.png", "RMM_vs_eigenfunction_YOTC_complex_labeled.png",
        "YOTC_hovmoller_olr_u850_with_koopman.png", "MJO_composites_phases=2-6_alignment=end.png",
        "MJO_precip_dist_phases=2-6_season=DJFMA_pct=75_alignment=end.png",
        "MJO_composites_precip_phases=2-6_season=DJFMA_pct=75_alignment=end.png",
        "BSISO_composites_phases=2-6_alignment=end.png", "RMM_skill_modes=j7_season=DJFMA.png"]
    supplement = ["RMM_OMI_FMO_vs_eigenfunction_YOTC_complex_labeled.png"]
    if include_large_composites
        append!(supplement, ["$(mode)_composites_$(var)_phases=1-8_alignment=end.png"
            for mode in ("MJO", "BSISO") for var in ("olr", "u850", "u200")])
        append!(supplement, ["MJO_precip_dist_phases=$(phases)_season=DJFMA_pct=75_alignment=end.png"
            for phases in ("1-8", "2-4-6-8")])
    end
    return vcat(joinpath.("paper", paper), joinpath.("supplement", supplement))
end

"""Render 11 figures, or 19 with larger composites. All index dates use end alignment."""
function make_all_paper_figures(; include_large_composites::Bool = false,
    output_root::AbstractString = DEFAULT_FIGURE_OUTPUT_ROOT,
    figure_width::Integer = PAPER_FIGURE_WIDTH, titles = DEFAULT_PLOT_TITLES)
    figure_width > 0 || throw(ArgumentError("figure_width must be positive"))
    required = [data_details(v)[3] for v in ("olr", "u", "precip")]
    append!(required, [paper_input("climdata", name) for name in
        ("RMM_data_edit.csv", "BSISO_kikuchi.csv", "omi.1x.txt", "fmo.1x.txt")])
    append!(required, [paper_cache(name) for name in ("init_Keigen.h5", "init_recon_RMMdata.h5",
        "init_recon_olr.h5", "init_recon_u.h5", "rmm_pred_skill_2020-12-01_2024-04-30.h5")])
    require_files(required)
    paper_dir, supplement_dir = joinpath(output_root, "paper"), joinpath(output_root, "supplement")
    paths = Dict(:paper => String[], :supplement => String[])
    _plot_eigenperiod!(paths; paper_dir, figure_width, title = titles.eigenperiod)
    GC.gc()
    _plot_rmm_index!(paths; paper_dir, figure_width, title = titles.rmm_complex)
    _plot_rmm_index_labeled!(paths; output_dir = paper_dir, figure_width, title = titles.rmm_complex,
        width_scale = 1.0, mark_dates = YOTC_MARK_DATES, region_labels = PHASE_REGION_LABELS)
    _plot_mjo_indices_labeled!(paths; output_dir = supplement_dir, figure_width,
        title = titles.mjo_indices_labeled, mark_dates = YOTC_MARK_DATES, region_labels = PHASE_REGION_LABELS)
    _plot_hovmollers!(paths; paper_dir, figure_width, title = titles.hov)
    GC.gc()
    _plot_mjo_figures!(paths; paper_dir, supplement_dir, figure_width, include_large_composites,
        composite_title = titles.mjo_composites, precip_title = titles.mjo_precip_dist,
        combined_title = titles.mjo_composites_precip)
    GC.gc()
    _plot_bsiso_composites!(paths; paper_dir, supplement_dir, figure_width,
        title = titles.bsiso_composites, include_large_composites)
    GC.gc()
    _plot_rmm_skill!(paths; paper_dir, figure_width, title = titles.rmm_skill)
    actual = vcat(paths[:paper], paths[:supplement])
    expected = joinpath.(output_root, expected_figure_paths(; include_large_composites))
    length(actual) == length(unique(actual)) || error("Duplicate figure paths")
    Set(actual) == Set(expected) || error("Figure inventory differs from the specification")
    all(isfile, actual) || error("Missing output figure")
    return paths
end

if abspath(PROGRAM_FILE) == @__FILE__
    all(==("--large-composites"), ARGS) || error("Usage: all_paper_plotting.jl [--large-composites]")
    make_all_paper_figures(; include_large_composites = "--large-composites" in ARGS)
end
