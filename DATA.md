# Prepared inputs

This example starts from the prepared dataset described below. It does not download or regrid the original climate products. Supply the files externally, preserving these paths relative to `PREDICTABLEMJO_INPUT_ROOT` (default: `inputs/`). `validate_inputs.jl` checks the 13 required files, variables, dimensions, units, grids, and date coverage. The paper grid and study period are the supported setup. Equivalent prepared files need not have the same byte encoding or hashes. `input_manifest.toml` records reference file sizes and SHA-256 hashes; checking those additionally requires `--checksums`.

| Relative path | Content | Size, MiB |
|---|---|---:|
| `data_clim_processed/OLR_feb25.nc` | ERA5 top net thermal radiation, `ttr` | 1011.72 |
| `data_clim_processed/u_feb25.nc` | ERA5 zonal wind, `u`, at 200 and 850 hPa | 1224.34 |
| `data_clim_processed/precip_sep25_02.nc` | ERA5 total precipitation, `tp` | 1011.72 |
| `data_clim_processed/merid_avgs_predict_6_2025.h5` | Normalized meridional averages used as model input | 66.90 |
| `climdata/RMM_data_edit.csv` | Daily RMM1/RMM2 and phase | 1.32 |
| `climdata/BSISO_kikuchi.csv` | Kikuchi BSISO `PCx`, `PCy`, phase and amplitudes | 0.70 |
| `climdata/omi.1x.txt` | Recorded NOAA OMI snapshot | 0.43 |
| `climdata/fmo.1x.txt` | Recorded NOAA FMO snapshot | 0.77 |
| `model_RMM_pred/merged_RMM/RMMS_ecmwf_realtime.nc` | ECMWF RMM forecasts | 5.98 |
| `model_RMM_pred/merged_RMM/RMMS_hmcr_realtime.nc` | HMCR RMM forecasts | 1.83 |
| `model_RMM_pred/merged_RMM/RMMS_kma_realtime.nc` | KMA RMM forecasts | 11.31 |
| `model_RMM_pred/merged_RMM/RMMS_ncep_realtime.nc` | NCEP RMM forecasts | 9.31 |
| `model_RMM_pred/merged_RMM/RMMS_ukmo_realtime.nc` | UKMO RMM forecasts | 12.31 |

Total: approximately **3.28 GiB**. No input is copied into or tracked by this repository.

## ERA5 fields and normalized features

The prepared fields contain 9,132 daily samples from **2000-01-01 through 2024-12-31**, on a global 1.5-degree grid: 240 longitudes from 0 to 358.5 degrees and 121 latitudes from -90 to 90 degrees. `valid_time` gives integer days since 2000-01-01. Julia reads OLR and precipitation as longitude × latitude × time, and wind as longitude × latitude × pressure × time, with pressure ordered `[200, 850]` hPa.

The plotting/reconstruction loaders apply the source workflow's conversions: `-ttr/3600` for outgoing longwave radiation in W m⁻², unchanged wind in m s⁻¹, and `1000*tp` for the prepared daily precipitation in mm day⁻¹. The radiation conversion assumes the source's one-hour accumulation convention. Equivalent unit notation such as `J/m^2` and `m s^-1` is accepted, but the physical quantities must match these conversions. Precipitation in metres must already represent the prepared daily total; unit metadata alone cannot establish the upstream aggregation. The conversions do not apply unchanged to arbitrary ERA5 products.

`olr_u` has Julia dimensions `(3, 240, 9132)`: normalized OLR, u200 and u850. It is based on area-weighted meridional averages over ±15 degrees latitude. The preprocessing represented by `normalizeData.jl` in the original project removes temporal means at each longitude and normalizes using a scalar standard deviation for OLR and another shared by the two wind levels. The prepared HDF5 stores `stdu`, `stdOLR`, `meanu`, and `meanOLR` alongside the features. This copy consumes the existing normalized array directly and does not recompute preprocessing.

The model uses 2000–2019 for fitting, a 64-day embedding, 1,500 neighbors, 501 diffusion modes and 101 modes in the regularized generator. The cone transition drops the first embedded row. Under end alignment, 7,241 eigenfunction rows correspond to 2000-03-05 through 2019-12-31.  Out-of-sample rows cover 2020–2024; forecast skill is evaluated from 2020-12-01 through 2024-04-30.

Mode 7 represents MJO and mode 12 BSISO; the eigenperiod comparison also marks mode 15. Variable composites use DJF (December–February) for MJO and JJA (June–August) for BSISO. Precipitation uses DJFMA (December–April), selecting amplitudes at or above each index's own seasonal 75th percentile and subtracting the seasonal mean. Larger composites retain these settings.

Forecast scores use verification-date masks for all months, DJFMA, and JJASO (June–October); the figure shows DJFMA. Koopman/KAF leads are 0–50 days. Full-delay indirect RMM reconstruction therefore predicts eigenfunction leads through day 113, always from the same issue date. Operational models retain their supplied lead grids.

## Reference indices and operational forecasts

The RMM snapshot uses the Wheeler–Hendon-style two-component index. Kikuchi BSISO is supplied separately as a dated CSV; this is not the APCC BSISO file also present in the original research directory. The selected code does not contain the original download/conversion scripts for these two CSVs, so their reference hashes identify the supplied examples. Equivalent observations with the required columns and date coverage are accepted.

The OMI and FMO sources recorded in the original code are:

- OMI: <https://psl.noaa.gov/mjo/mjoindex/omi.1x.txt>
- FMO: <https://psl.noaa.gov/mjo/mjoindex/fmo.1x.txt>

The reference snapshots were retrieved on 2026-09-08; their optional SHA-256 checksums are recorded in `input_manifest.toml`.

Their parsers reject missing/duplicate dates, invalid amplitudes, and incompatible formats. Plotting accepts valid supplied snapshots and never downloads replacements. Daily observations are joined by date; OMI and FMO use `(PC2, -PC1)` for geographic phase orientation, without filtering, rescaling, or fitting to the Koopman index.

The operational forecast files are the original merged S2S RMM products. Each supplies `RMM1`, `RMM2`, `step`, and `valid_time` with its time-unit attribute. Original acquisition and merge scripts are outside this copy; the manifest records the reference inputs. Forecast times must use days since a stated epoch, and leads are in days. Date-only and midnight timestamp forms of the epoch are accepted. Model-specific lead grids and verification dates are preserved.

## Generated results and data distribution

The default pipeline produces diffusion and Koopman caches, mode-7 OLR/wind/RMM reconstructions, the nonnegative-frequency normalized-feature reconstructions needed for eigenperiod coloring, out-of-sample kernel/eigenfunction evaluations, and the RMM skill cache. Raw climate fields and generated caches are not Git assets.

New NLSA caches save only the four datasets needed downstream. Distance matrices, neighbor masks, pairwise bandwidths, and the normalized diffusion matrix are still computed but not saved. Older caches with extra datasets remain readable; reconstruction and OSE formats are unchanged.

All dimensions below use Julia's array order:

| Cache dataset | Dimensions and convention |
|---|---|
| NLSA `W`, `bw`, `phi`, `kappa` | 7,241 × 7,241 kernel, scalar bandwidth, 7,241 × 501 diffusion basis, and 501 eigenvalues |
| Koopman `zetaL`, `zetaR` | 7,241 end-dated samples × 101 modes; `c` is 101 × 101 |
| Reconstruction `recons` | selected mode × flattened spatial coordinate × 7,304 days, from 2000-01-02 through 2019-12-31 |
| Reconstruction `Acoeff` | selected mode × flattened spatial coordinate × 64 delay coordinates, ordered oldest to newest; `indices` records the selected mode numbers |
| OSE `p_evals`, `phi_vals` | 1,827 issue dates × 7,241 training samples, and 101 modes × 1,827 issue dates, respectively |
| OSE `zeta_vals`, `zetaR_vals` | 51 modes × 1,827 issue dates; `dates` contains every date and `date_convention` is `cone_transition_end` |
| Skill `skill/<curve>/corr` | lead × season, with season names and lead values stored separately; operational models retain their own lead grids |

Normalized-feature and RMM reconstructions are dimensionless; OLR and wind reconstructions retain the physical units above. Reconstructions store complex modal contributions; the retained nonzero-frequency mode is displayed as twice its real part. Projection uses the uniform sample weight `1/7241`, and reconstruction averages the available delay coordinates at each day, including shorter averages at the edges. Stored generator eigenvalues use inverse days with angular frequency in their imaginary parts; `getK` divides by `2π` to obtain cycles per day.

A small dataset containing only plotted arrays could support fast figure rendering: eigenperiod statistics, index trajectories, Hovmöller arrays, composites, precipitation means, and skill curves. The larger composite options require more maps than the original paper-only estimate. An export has deliberately not been made; measure and validate its size before deciding whether to add it to GitHub. Such a dataset would not replace the prepared inputs for recomputing the scientific results.
