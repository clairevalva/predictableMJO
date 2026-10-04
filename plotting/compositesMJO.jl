# Phase subsets appear in filenames: contiguous runs collapse to "1-8", the rest join as "2-4-6-8".
function _phase_tag(phases)
    phases = collect(phases)
    if length(phases) >= 3 && phases == collect(first(phases):last(phases))
        return "$(first(phases))-$(last(phases))"
    end
    return join(phases, "-")
end

function _reference_index(reference_kind::Symbol, dates)
    if reference_kind === :rmm
        _, _, reference, _ = getRMM(dates[1], dates[end])
        return reference, "RMM"
    elseif reference_kind === :kikuchi
        _, _, reference = getBSISOkikuchi(dates[1], dates[end]; all_vars = false)
        return reference, "BSISO"
    end
    throw(ArgumentError("Unsupported reference kind $(reference_kind)"))
end

function _compute_index_composites(; data_name::String, level_index::Int,
    season::String, eigenfunction_index::Int, reference_kind::Symbol,
    reference_loader = _reference_index)
    units, nlevels, data_fn, data_key, _, multiplier = data_details(data_name)
    X, lon, lat, levels, field_times = get_data(data_fn, data_key, multiplier; retall = true)
    tStart, tEnd, _, nEmb, datastart = getExperiment("init")
    dates = datastart + Day.((tStart:tEnd) .- 1)
    _, _, zetaL, _, periods = getK(paper_cache("init_Keigen.h5"))
    nT = size(zetaL, 1)
    time_inds = makeTakensInd_forward(nEmb, nT)
    takens_alignment = :end
    @assert length(time_inds) == nT

    embedded_dates = dates[time_inds]
    field_inds = select_dates(PAPER.start_date .+ Day.(field_times), embedded_dates)
    X = reshape(X, length(lon), length(lat), nlevels, :)
    # Match physical dates; the experiment begins on Jan 2, fields on Jan 1.
    X = X[:, :, :, field_inds]
    embedded_dates = dates[time_inds]

    day_numbers = dayofyear.(embedded_dates)
    daily_mean = zeros(length(lon), length(lat), nlevels, 366)
    for day in 1:366
        day_mask = day_numbers .== day
        if any(day_mask)
            daily_mean[:, :, :, day] = mean(X[:, :, :, day_mask]; dims = 4)
        end
    end
    for t in eachindex(day_numbers)
        X[:, :, :, t] .-= daily_mean[:, :, :, day_numbers[t]]
    end

    reference, reference_name = reference_loader(reference_kind, dates)
    reference = reference[time_inds]
    zeta, rotation = align_index(reference, zetaL[:, eigenfunction_index])
    println("$(reference_name) alignment rotation ($(takens_alignment)): ", rotation)
    println("eigenfunction period: ", periods[eigenfunction_index], " days")

    koopman_phases = makePhases(zeta)
    reference_phases = makePhases(reference)
    season_order = ["all", "DJF", "MAM", "JJA", "SON"]
    season_position = findfirst(==(season), season_order)
    isnothing(season_position) && throw(ArgumentError("Unsupported season $(season)"))
    if season_position == 1
        season_mask = trues(nT)
    else
        season_mask = makeseasons(embedded_dates)[season_position - 1, :]
    end

    output_shape = (8, length(lon), length(lat), nlevels)
    koopman_composites = zeros(output_shape)
    reference_composites = zeros(output_shape)
    for phase in 1:8
        koopman_mask = season_mask .& koopman_phases[phase, :]
        reference_mask = season_mask .& reference_phases[phase, :]
        any(koopman_mask) || error("No Koopman samples for phase $(phase), season $(season)")
        any(reference_mask) || error("No reference samples for phase $(phase), season $(season)")
        koopman_composites[phase, :, :, :] = mean(X[:, :, :, koopman_mask]; dims = 4)
        reference_composites[phase, :, :, :] = mean(X[:, :, :, reference_mask]; dims = 4)
    end

    return (; data_name, level_index, season, eigenfunction_index, reference_name,
        takens_alignment, units, levels, lon, lat, koopman_composites,
        reference_composites)
end

function _composite_extent(result, phases)
    return max(
        maximum(abs, result.koopman_composites[phases, :, :, result.level_index]),
        maximum(abs, result.reference_composites[phases, :, :, result.level_index]),
    )
end

function _composite_colorbar_label(result)
    level_label = if result.data_name == "u"
        " at $(Int(result.levels[result.level_index])) hPa"
    else
        ""
    end
    return "$(uppercase(result.data_name)) anomalies [$(result.units)]$(level_label)"
end

# Width-to-height ratio the composite maps actually render at, measured from the GeoAxis
# viewport. Pure Mercator geometry for ±30° gives 5.72, but the axis autolimit margins and
# the 0-357.5° data range stretch it; at 5.72 the maps come out 4% narrower than their cell.
const _MERCATOR_MAP_ASPECT = 5.47

# Longitude ticks are not an option here: GeoMakie's tick projection recurses until it
# overflows the stack on this map, for every tick specification including its own default.
# In the combined figure the precip panel below carries the longitude axis instead.
function _composite_map_axis(position; title)
    return GeoMakie.GeoAxis(position; title,
        dest = "+proj=merc +lon_0=180", xgridvisible = false, ygridvisible = false,
        xticklabelsvisible = false, yticklabelsvisible = false,
        xticksvisible = false, yticksvisible = false, titlesize = 12)
end

# Draws the Koopman and reference maps for one phase into the given cells, and returns a
# heatmap handle for the shared colorbar.
function _draw_composite_maps!(koopman_cell, reference_cell, result, phase, extent)
    hm = nothing
    panels = (
        (koopman_cell, result.koopman_composites, "Phase $(phase) Koopman"),
        (reference_cell, result.reference_composites,
            "Phase $(phase) $(result.reference_name) index"),
    )
    for (cell, composites, title) in panels
        ax = _composite_map_axis(cell; title)
        hm = GeoMakie.heatmap!(ax, result.lon, result.lat,
            composites[phase, :, :, result.level_index];
            colormap = Reverse(:tarn), colorrange = (-extent, extent))
        GeoMakie.lines!(ax, GeoMakie.coastlines(ax); color = :black, linewidth = 1)
        CairoMakie.ylims!(ax, -30, 30)
    end
    return hm
end

# One phase per row: Koopman in column 1, reference index in column 2.
function _draw_composite_row!(fig, row, result, phase, extent)
    return _draw_composite_maps!(fig[row, 1], fig[row, 2], result, phase, extent)
end

function _render_index_composites(result, phases, path; figure_width, title)
    extent = _composite_extent(result, phases)
    height = max(190, round(Int, figure_width * 80 * length(phases) / 600))
    fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 12)
    variable = result.data_name == "u" ? "u, $(Int(result.levels[result.level_index])) hPa" :
        uppercase(result.data_name)
    CairoMakie.Label(fig[0, 1:2], "$(title) ($(variable))"; fontsize = 14, font = :bold)
    hm = nothing
    for (row, phase) in enumerate(phases)
        hm = _draw_composite_row!(fig, row, result, phase, extent)
    end
    CairoMakie.Colorbar(fig[:, 3], hm; label = _composite_colorbar_label(result),
        height = CairoMakie.Relative(0.5), width = 10)
    CairoMakie.colgap!(fig.layout, 10)
    CairoMakie.rowgap!(fig.layout, 0)
    CairoMakie.save(path, fig; px_per_unit = 500 / 96)
    return path
end

# Reuse each computed field across phase subsets and, for wind, pressure levels.
function _save_large_variable_composites!(paths, result; mode_name, supplement_dir, figure_width, title)
    levels = result.data_name == "u" ? (1, 2) : (1,)
    for level_index in levels
        panel = merge(result, (; level_index))
        var = result.data_name == "u" ? "u$(Int(result.levels[level_index]))" : result.data_name
        path = joinpath(supplement_dir, "$(mode_name)_composites_$(var)_phases=1-8_alignment=end.png")
        push!(paths[:supplement], _render_index_composites(panel, 1:8, path; figure_width, title))
    end
end

function _plot_mjo_figures!(paths; paper_dir, supplement_dir, figure_width,
    composite_title, precip_title, combined_title, include_large_composites = false)
    mkpath(paper_dir)
    include_large_composites && mkpath(supplement_dir)
    composites = _compute_index_composites(; data_name = "olr", level_index = 1,
        season = "DJF", eigenfunction_index = 7, reference_kind = :rmm)
    precip = _compute_conditional_precip(; eigenfunction_index = 7, season = "DJFMA", percentile_want = 75)
    phases, ylims = [2, 6], (-3.0, 5.0)
    push!(paths[:paper], _render_index_composites(composites, phases,
        joinpath(paper_dir, "MJO_composites_phases=2-6_alignment=end.png"); figure_width, title = composite_title))
    precip_name(tag) = "MJO_precip_dist_phases=$(tag)_season=DJFMA_pct=75_alignment=end.png"
    push!(paths[:paper], _render_precip_dist(precip, phases, joinpath(paper_dir, precip_name("2-6"));
        figure_width, title = precip_title, ylims))
    push!(paths[:paper], _render_composites_with_precip(composites, precip, phases,
        joinpath(paper_dir, "MJO_composites_precip_phases=2-6_season=DJFMA_pct=75_alignment=end.png");
        figure_width, title = combined_title, ylims))
    if include_large_composites
        _save_large_variable_composites!(paths, composites; mode_name = "MJO", supplement_dir,
            figure_width, title = composite_title)
        for ps in (collect(1:8), [2, 4, 6, 8])
            push!(paths[:supplement], _render_precip_dist(precip, ps,
                joinpath(supplement_dir, precip_name(_phase_tag(ps))); figure_width, title = precip_title, ylims))
        end
        wind = _compute_index_composites(; data_name = "u", level_index = 1,
            season = "DJF", eigenfunction_index = 7, reference_kind = :rmm)
        _save_large_variable_composites!(paths, wind; mode_name = "MJO", supplement_dir,
            figure_width, title = composite_title)
    end
    return nothing
end
