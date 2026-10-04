# predictableMJO


A somewhat minimal Julia example of the Koopman MJO paper experiments [arxiv2609.19435](https://arxiv.org/abs/2609.19435). This code was ported from the original codebase avoid having a public repo with large data files and outdated experiments/numeircal tests that are not essential to the reproduction of the paper. 

**Note: This code was produced with the assistance of GPT 5.5 and 5.6, and the port from the old repo to the new was completed with the assistance of GPT-6 Astra.**

## Setup

Tested with Julia **1.11.5**. The project and manifest provide a known package environment, including CairoMakie and GeoMakie. The goal is to explain and support replication of the experiments; bit-for-bit agreement with archived results is not required.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
```

Place the separately supplied files under `inputs/` in the layout documented in [DATA.md](DATA.md), or set `PREDICTABLEMJO_INPUT_ROOT` to a directory with that layout. No data, figures, external symlinks, or source Git history are included. Equivalent prepared files are accepted when their variables, dimensions, units, and dates match the paper setup.

```sh
export PREDICTABLEMJO_INPUT_ROOT=/path/to/prepared-inputs
julia --project=. validate_inputs.jl
julia --project=. run_paper.jl
```

Validation checks the prepared-data format by default. To additionally verify the reference snapshots, run `julia --project=. validate_inputs.jl --checksums`, or add `--checksums` to the full runner. Files are never downloaded or replaced automatically.

The full workflow builds the diffusion basis, Koopman eigenfunctions, required reconstructions, out-of-sample evaluations, forecast-skill curves, and figures. Kernel fitting and out-of-sample evaluation are substantial computations. The original jobs requested 100–200 GB RAM; this is historical provisioning, not a measured minimum for this copy.

Generated intermediates go to `data/` and figures to `figures/`. Override these independently with `PREDICTABLEMJO_CACHE_ROOT` and `PREDICTABLEMJO_OUTPUT_ROOT`. Run stages in fresh Julia processes to release memory between them. Existing intermediate files are protected unless `--force` is passed; `--reuse-cache` explicitly skips stages whose outputs already exist and assumes their settings match this project.

```sh
julia --project=. run_paper.jl --reuse-cache --large-composites
# Render from existing intermediate results without refitting:
julia --project=. plotting/all_paper_plotting.jl
julia --project=. plotting/all_paper_plotting.jl --large-composites
```

Individual stages remain runnable: `doNLSA.jl`, `doKoopman.jl`, `doReconstruction.jl`, `computeOSE.jl`, and `model_comparisons/rmm_prediction_skill.jl`. Reconstruction accepts `RMMdata`, `RMM`, `olr`, or `u` as positional arguments; no arguments builds all four. `paper_config.jl` contains the shared paper settings.

## Figures

The default inventory is two eigenperiod figures; two RMM–Koopman complex-plane figures; the labeled RMM/OMI/FMO comparison; the OLR/u850 YOTC Hovmöller; MJO phases 2/6 composites, precipitation profiles, and their combined figure; a BSISO phases 2/6 composite; and the j7/DJFMA forecast-skill figure. `expected_figure_paths()` lists their exact relative filenames.

The optional larger set includes MJO and BSISO phases 1–8 composites of OLR, u200 and u850, plus precipitation profiles for phases 1–8 and 2/4/6/8. Results for each variable are reused across phase subsets and wind levels.

