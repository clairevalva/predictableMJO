function _load_hovmoller_panel(data_name::String, level_index::Int;
    title::String, cache::Dict,
    plot_start::Date = Date(2009, 10, 1), plot_end::Date = Date(2010, 4, 1),
    mark_dates::Vector{Date} = Date[])
    repo_dir = @__DIR__
    tStart, tEnd, _, _, datastart = getExperiment("init")
    reconstruction_dates = (datastart + Day.(0:tEnd))[tStart:tEnd]

    if !haskey(cache, data_name)
        units, nlevels, data_fn, data_key, pressures, multiplier = data_details(data_name)
        X, lon, lat, _, X_time = get_data(data_fn, data_key, multiplier; retall = true)
        X = reshape(X, length(lon), length(lat), nlevels, length(X_time))
        observation_dates = Date(2000, 1, 1) .+ Day.(X_time)
        rec_file = paper_cache("init_recon_$(data_name).h5")
        recs, _ = get_reconstructed_data(rec_file; load_small_j = [7])
        recs = reshape(recs, size(recs, 1), length(lon), length(lat), nlevels, size(recs, 3))
        cache[data_name] = (; units, nlevels, pressures, X, lon, lat,
            observation_dates, recs)
    end
    common = cache[data_name]
    units = common.units
    nlevels = common.nlevels
    pressures = common.pressures
    X = common.X
    lon = common.lon
    lat = common.lat
    observation_dates = common.observation_dates
    recs = common.recs
    field = X[:, :, level_index, :]

    weights = areaweights(lat, lon)
    observed = makemeridavg(field .- mean(field; dims = 3), 15, lat, weights)
    bandpass = makeMJOFF(; days_low = 90, days_high = 30, sample_fs = 1, lanN = 49)
    observed_bandpass = permutedims(bandpass(permutedims(observed)))

    obs_inds = findall(d -> plot_start <= d <= plot_end, observation_dates)
    rec_inds = findall(d -> plot_start <= d <= plot_end, reconstruction_dates)
    isempty(obs_inds) && error("No $(data_name) observations found between $(plot_start) and $(plot_end)")
    length(obs_inds) == length(rec_inds) || error("Observation/reconstruction lengths differ")

    Z = observed_bandpass[:, obs_inds]
    rec = 2 .* real.(recs[1, :, :, level_index, :])
    reconstructed = makemeridavg(rec, 15, lat, weights)
    Zrec = reconstructed[:, rec_inds]
    fill_extent = maximum(abs, Z)
    rec_extent = maximum(abs, Zrec)
    fill_levels = range(-fill_extent, fill_extent; length = 21)
    rec_levels = collect(range(-rec_extent, rec_extent; length = 11))
    dates = observation_dates[obs_inds]
    wanted_ticks = (plot_start:Month(1):plot_end) .+ Day(14)
    yinds = collect(eachindex(dates))
    tick_values = [i for i in yinds if dates[i] in wanted_ticks]
    tick_labels = Dates.format.(dates[tick_values], "u. yyyy")
    pressure = isnothing(pressures) ? nothing : Int(pressures[level_index])
    label = data_name == "u" ? "u" : uppercase(data_name)
    cmap = data_name == "u" ? :balance : Reverse(:tarn)

    # Optional horizontal time-axis markers (e.g. 25 Oct / 20 Dec) drawn as row indices.
    mark_yinds = Int[]
    mark_labels = String[]
    for md in mark_dates
        pos = findfirst(==(md), dates)
        pos === nothing && continue
        push!(mark_yinds, pos)
        push!(mark_labels, Dates.format(md, "d u"))
    end

    return (; data_name, level_index, title, label, units, pressure, lon, yinds,
        tick_values, tick_labels, Z, Zrec, fill_levels, rec_levels, cmap,
        plot_start, plot_end, mark_yinds, mark_labels)
end

_lon_tick_label(x) = x == 0 ? "0°" : "$(Int(x))°E"

# Colour of the date marker lines on the dated Hovmöller (edit here to restyle).
const _MARK_LINE_COLOR = :gray35

function _draw_hov_panel!(layout, panel; overlay::Bool, show_ylabels::Bool = false,
    titlesize::Real = 12)
    ax = CairoMakie.Axis(layout;
        title = panel.title,
        titlesize = titlesize,
    )
    hm = CairoMakie.contourf!(ax, panel.lon, panel.yinds, panel.Z;
        colormap = panel.cmap, levels = panel.fill_levels)
    if overlay
        positive = filter(>(0), panel.rec_levels)
        negative = filter(<(0), panel.rec_levels)
        isempty(positive) || CairoMakie.contour!(ax, panel.lon, panel.yinds, panel.Zrec;
            levels = positive, color = :black, linewidth = 1, linestyle = :solid)
        isempty(negative) || CairoMakie.contour!(ax, panel.lon, panel.yinds, panel.Zrec;
            levels = negative, color = :black, linewidth = 1.5, linestyle = :dash)
    end
    ax.yticks = (panel.tick_values, show_ylabels ? panel.tick_labels : fill("", length(panel.tick_values)))
    ax.yticklabelrotation = pi / 12
    lon_ticks = [0, 100, 200, 300]
    ax.xticks = (lon_ticks, _lon_tick_label.(lon_ticks))
    ax.xticklabelrotation = pi / 4
    CairoMakie.ylims!(ax, 1, length(panel.yinds))
    CairoMakie.xlims!(ax, extrema(panel.lon)...)
    if !isempty(panel.mark_yinds)
        # Unobtrusive solid lines: grey keeps them distinct from the black reconstruction
        # contours without competing with the field. Dates are identified in the caption.
        CairoMakie.hlines!(ax, panel.mark_yinds; color = _MARK_LINE_COLOR,
            linestyle = :solid, linewidth = 1.4)
    end
    return hm
end

function _save_combined_hov(path, panels; figure_width, title, titlesize::Real = 12)
    height = round(Int, figure_width * 0.8)
    fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 10)
    CairoMakie.Label(fig[0, 1:length(panels)], title; fontsize = 14, font = :bold)
    for (i, panel) in enumerate(panels)
        hm = _draw_hov_panel!(fig[1, i], panel; overlay = true, show_ylabels = i == 1,
            titlesize)
        CairoMakie.Colorbar(fig[2, i], hm; vertical = false, flipaxis = false,
            label = "$(panel.label) [$(panel.units)]", width = CairoMakie.Relative(0.8))
    end
    CairoMakie.save(path, fig; px_per_unit = 5)
    return path
end

function _plot_hovmollers!(paths; paper_dir, figure_width, title)
    cache = Dict()
    olr = _load_hovmoller_panel("olr", 1; title = "OLR (bandpassed)", cache)
    u850 = _load_hovmoller_panel("u", 2; title = "Zonal wind at 850 hPa (bandpassed)", cache)
    mkpath(paper_dir)
    path = joinpath(paper_dir, "YOTC_hovmoller_olr_u850_with_koopman.png")
    push!(paths[:paper], _save_combined_hov(path, [olr, u850]; figure_width, title))
    return nothing
end
