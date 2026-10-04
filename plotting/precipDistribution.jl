# Conditional precipitation distributions: meridionally averaged precip anomaly against longitude,
# composited by phase and conditioned on index amplitude.

const _PRECIP_SEASON_MONTHS = Dict(
    "DJFMA" => [12, 1, 2, 3, 4],
    "JJASO" => [6, 7, 8, 9, 10],
)

function _precip_season_mask(season::AbstractString, dates)
    season == "all" && return trues(length(dates))
    haskey(_PRECIP_SEASON_MONTHS, season) ||
        throw(ArgumentError("Unsupported season $(season)"))
    months = _PRECIP_SEASON_MONTHS[season]
    return [month(date) in months for date in dates]
end

function _compute_conditional_precip(; eigenfunction_index::Int, season::String,
    percentile_want::Int, reference_kind::Symbol = :rmm,
    data_name::String = "precip", expname::String = "init",
    reference_loader = _reference_index)
    units, nlevels, data_fn, data_key, _, multiplier = data_details(data_name)
    X, lon, lat, _, field_times = get_data(data_fn, data_key, multiplier; retall = true)
    tStart, tEnd, _, nEmb, datastart = getExperiment(expname)
    dates = datastart + Day.((tStart:tEnd) .- 1)
    _, _, zetaL, _, periods = getK(paper_cache("$(expname)_Keigen.h5"))
    nT = size(zetaL, 1)
    time_inds = makeTakensInd_forward(nEmb, nT)
    takens_alignment = :end
    @assert length(time_inds) == nT
    # The same physical-time mapping has to reach the field, the reference index, and the seasons.
    embedded_dates = dates[time_inds]
    field_inds = select_dates(PAPER.start_date .+ Day.(field_times), embedded_dates)
    X = reshape(X, length(lon), length(lat), nlevels, :)
    # Match physical dates; the experiment begins on Jan 2, fields on Jan 1.
    X = X[:, :, 1, field_inds]
    embedded_dates = dates[time_inds]

    reference, reference_name = reference_loader(reference_kind, dates)
    reference = reference[time_inds]
    zeta, rotation = align_index(reference, zetaL[:, eigenfunction_index])
    println("$(reference_name) alignment rotation ($(takens_alignment)): ", rotation)
    println("eigenfunction period: ", periods[eigenfunction_index], " days")

    season_mask = _precip_season_mask(season, embedded_dates)
    any(season_mask) || error("No samples in season $(season)")
    koopman_phases = makePhases(zeta)
    reference_phases = makePhases(reference)

    # Amplitude thresholds are percentiles within the season, taken separately for each index.
    koopman_threshold = percentile(abs.(zeta[season_mask]), percentile_want)
    reference_threshold = percentile(abs.(reference[season_mask]), percentile_want)
    koopman_strong = abs.(zeta) .>= koopman_threshold
    reference_strong = abs.(reference) .>= reference_threshold

    season_mean = mean(X[:, :, season_mask]; dims = 3)
    koopman_fields = zeros(length(lon), length(lat), 8)
    reference_fields = zeros(length(lon), length(lat), 8)
    w_area = areaweights(lat, lon)
    koopman_counts = zeros(Int, 8)
    reference_counts = zeros(Int, 8)
    for phase in 1:8
        koopman_mask = koopman_phases[phase, :] .& season_mask .& koopman_strong
        reference_mask = reference_phases[phase, :] .& season_mask .& reference_strong
        any(koopman_mask) || error("No Koopman samples for phase $(phase), season $(season)")
        any(reference_mask) || error("No reference samples for phase $(phase), season $(season)")
        koopman_counts[phase] = count(koopman_mask)
        reference_counts[phase] = count(reference_mask)
        koopman_fields[:, :, phase] = mean(X[:, :, koopman_mask]; dims = 3) .- season_mean
        reference_fields[:, :, phase] = mean(X[:, :, reference_mask]; dims = 3) .- season_mean
    end

    koopman_lines = makemeridavg(koopman_fields, 15, lat, w_area)
    reference_lines = makemeridavg(reference_fields, 15, lat, w_area)

    return (; lon, units, season, percentile_want, takens_alignment, eigenfunction_index,
        reference_name, koopman_lines, reference_lines, koopman_counts, reference_counts)
end

function _precip_ylims(result, phases, ylims)
    isnothing(ylims) || return ylims
    extent = max(
        maximum(abs, result.koopman_lines[:, phases]),
        maximum(abs, result.reference_lines[:, phases]),
    )
    return (-extent, extent)
end

# Panels are short, so auto ticks crowd; aim for about four round ticks inside the limits.
function _precip_yticks(limits)
    low, high = limits
    step = max(1.0, round((high - low) / 4))
    return collect((ceil(low / step) * step):step:(floor(high / step) * step))
end

function _draw_precip_panel!(ax, result, phase, limits)
    reference_label = result.reference_name
    koopman_label = "Koopman"
    CairoMakie.lines!(ax, result.lon, result.reference_lines[:, phase];
        color = :blue, linewidth = 2, label = reference_label)
    CairoMakie.lines!(ax, result.lon, result.koopman_lines[:, phase];
        color = :red, linewidth = 2, label = koopman_label)
    CairoMakie.hlines!(ax, [0.0]; color = :black, linestyle = :dash)
    CairoMakie.ylims!(ax, limits...)
    return ax
end

function _precip_axis(position, result; title, show_xlabel, limits)
    return CairoMakie.Axis(position; title, titlesize = 12,
        xlabel = "Longitude (deg)", xlabelvisible = show_xlabel,
        xticklabelsvisible = show_xlabel,
        ylabel = result.units, yticks = _precip_yticks(limits))
end

function _render_precip_dist(result, phases, path; figure_width, title, ylims)
    limits = _precip_ylims(result, phases, ylims)
    height = round(Int, 110 * length(phases)) + 60
    fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 12)
    # tellwidth = false, or the single Auto column collapses to the title's text width.
    CairoMakie.Label(fig[0, 1], title; fontsize = 14, font = :bold, tellwidth = false)
    axis_last = nothing
    for (row, phase) in enumerate(phases)
        ax = _precip_axis(fig[row, 1], result; title = "Phase $(phase)",
            show_xlabel = row == length(phases), limits)
        axis_last = _draw_precip_panel!(ax, result, phase, limits)
    end
    CairoMakie.Legend(fig[length(phases) + 1, 1], axis_last;
        orientation = :horizontal, tellheight = true, tellwidth = false,
        framevisible = false)
    CairoMakie.rowgap!(fig.layout, 4)
    CairoMakie.save(path, fig; px_per_unit = 500 / 96)
    return path
end

# One column per phase, each holding the Koopman map, the reference map, and the precip
# profile. Sharing a column is what gives the three panels an identical plot width, so the
# longitude axes match without linking them.
function _render_composites_with_precip(composites, precip, phases, path;
    figure_width, title, ylims)
    limits = _precip_ylims(precip, phases, ylims)
    extent = _composite_extent(composites, phases)
    ncols = length(phases)
    map_height = figure_width / ncols / _MERCATOR_MAP_ASPECT
    height = round(Int, 2 * map_height + 250)
    fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 12,
        figure_padding = 6)
    CairoMakie.Label(fig[0, 1:ncols], title; fontsize = 14, font = :bold)
    hm = nothing
    axis_last = nothing
    for (col, phase) in enumerate(phases)
        hm = _draw_composite_maps!(fig[1, col], fig[2, col], composites, phase, extent)
        ax = _precip_axis(fig[4, col], precip; title = "Phase $(phase) precip",
            show_xlabel = true, limits)
        # Columns share y limits, so only the leftmost needs the axis.
        col == 1 || CairoMakie.hideydecorations!(ax; grid = false)
        axis_last = _draw_precip_panel!(ax, precip, phase, limits)
    end
    # Each key sits directly under the panels it describes: colorbar under the maps,
    # line legend under the precip profiles.
    CairoMakie.Colorbar(fig[3, 1:ncols], hm;
        label = _composite_colorbar_label(composites), vertical = false,
        flipaxis = false, width = CairoMakie.Relative(0.5), height = 10)
    CairoMakie.Legend(fig[5, 1:ncols], axis_last;
        orientation = :horizontal, labelsize = 9, patchsize = (14.0f0, 2.0f0),
        framevisible = false, tellheight = true, tellwidth = false,
        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0))
    # Aspect ties each map row's height to the real column width, so the maps fill their
    # cells instead of shrinking and floating.
    for row in 1:2
        CairoMakie.rowsize!(fig.layout, row, CairoMakie.Aspect(1, 1 / _MERCATOR_MAP_ASPECT))
    end
    CairoMakie.colgap!(fig.layout, 12)
    CairoMakie.rowgap!(fig.layout, 6)
    CairoMakie.save(path, fig; px_per_unit = 500 / 96)
    return path
end
