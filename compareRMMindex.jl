# Shared YOTC-period data prep for the complex-plane figures. Returns the plotted RMM and aligned
# Koopman series using end dates for every embedded row.
function _yotc_complex_series(; plot_start::Date = Date(2009, 10, 1),
    plot_end::Date = Date(2010, 4, 1))
    tStart, tEnd, _, nEmb, datastart = getExperiment("init")
    dates = datastart + Day.((tStart:tEnd) .- 1)
    _, _, zetaL, _, _ = getK(paper_cache("init_Keigen.h5"))
    nT = size(zetaL, 1)
    time_inds = makeTakensInd_forward(nEmb, nT)
    embedded_dates = dates[time_inds]

    _, _, rmm_full, _ = getRMM(dates[1], dates[end])
    rmm = rmm_full[time_inds]
    zeta, _ = align_index(rmm, zetaL[:, 7])
    selected = findall(d -> plot_start <= d <= plot_end, embedded_dates)
    isempty(selected) && error("No embedded dates found in the YOTC plotting period")

    return (; dates_plot = embedded_dates[selected], rmm_plot = rmm[selected],
        zeta_plot = zeta[selected], rmm_full, zetaL, nEmb, nT, dates, plot_start, plot_end)
end

# Square, shared limits so both panels use an identical (centered) data range.
function _square_complex_limits(rmm_plot, zeta_plot)
    all_x = vcat(real.(rmm_plot), real.(zeta_plot))
    all_y = vcat(imag.(rmm_plot), imag.(zeta_plot))
    center_x = sum(extrema(all_x)) / 2
    center_y = sum(extrema(all_y)) / 2
    half_span = max(maximum(all_x) - minimum(all_x), maximum(all_y) - minimum(all_y)) / 2
    return (center_x - half_span, center_x + half_span), (center_y - half_span, center_y + half_span)
end

function _plot_rmm_index!(paths; paper_dir, figure_width, title, width_scale = 2 / 3)
    series = _yotc_complex_series()
    dates_plot = series.dates_plot
    rmm_plot = series.rmm_plot
    zeta_plot = series.zeta_plot
    time_color = collect(eachindex(dates_plot))
    colorrange = extrema(time_color)
    xlimits, ylimits = _square_complex_limits(rmm_plot, zeta_plot)

    # Narrower than the other paper figures. The panels are DataAspect, so the figure has to
    # be tall enough for a square to reach the column width; otherwise height caps the panel
    # size and leaves a gap between the two.
    paper_width = round(Int, figure_width * width_scale)
    height = round(Int, paper_width * 0.66)
    fig = CairoMakie.Figure(size = (paper_width, height), fontsize = 12)
    CairoMakie.Label(fig[0, 1:2], title; fontsize = 14, font = :bold)
    ax1 = CairoMakie.Axis(fig[1, 1]; xlabel = "RMM1", ylabel = "RMM2",
        title = "RMM index", aspect = CairoMakie.DataAspect())
    ax2 = CairoMakie.Axis(fig[1, 2]; xlabel = "Re(ζ)", ylabel = "Im(ζ)",
        title = "Koopman eigenfunction", aspect = CairoMakie.DataAspect())
    CairoMakie.lines!(ax1, real.(rmm_plot), imag.(rmm_plot);
        color = time_color, colormap = :hawaii, colorrange = colorrange)
    CairoMakie.lines!(ax2, real.(zeta_plot), imag.(zeta_plot);
        color = time_color, colormap = :hawaii, colorrange = colorrange)
    CairoMakie.scatter!(ax1, [real(first(rmm_plot))], [imag(first(rmm_plot))];
        color = :black, marker = :circle, markersize = 6)
    CairoMakie.scatter!(ax2, [real(first(zeta_plot))], [imag(first(zeta_plot))];
        color = :black, marker = :circle, markersize = 6)
    for ax in (ax1, ax2)
        CairoMakie.xlims!(ax, xlimits...)
        CairoMakie.ylims!(ax, ylimits...)
    end
    CairoMakie.colgap!(fig.layout, 6)
    mkpath(paper_dir)
    paper_path = joinpath(paper_dir, "RMM_vs_eigenfunction_YOTC_complex_timecolored.png")
    CairoMakie.save(paper_path, fig; px_per_unit = 5)
    push!(paths[:paper], paper_path)

    return nothing
end

# Draw the eight-phase octant boundaries (grey lines through the origin) plus the four MJO regional
# labels. Applied identically to both panels since RMM and the aligned ζ share the phase geometry.
function _draw_phase_regions!(ax, xlimits, ylimits, region_labels; label_fontsize = 12)
    reach = maximum(abs, (xlimits..., ylimits...))
    for (xs, ys) in (([-reach, reach], [0.0, 0.0]),        # Re axis
                     ([0.0, 0.0], [-reach, reach]),        # Im axis
                     ([-reach, reach], [-reach, reach]),   # +45° diagonal
                     ([-reach, reach], [reach, -reach]))   # -45° diagonal
        CairoMakie.lines!(ax, xs, ys; color = (:gray, 0.4), linewidth = 0.8)
    end
    xr = xlimits[2] - xlimits[1]
    yr = ylimits[2] - ylimits[1]
    xc = (xlimits[1] + xlimits[2]) / 2
    yc = (ylimits[1] + ylimits[2]) / 2
    # Region labels follow the Wheeler-Hendon (2004) convention, matching makePhases: Western
    # Pacific top, Maritime Continent right, Indian Ocean bottom, Western Hemisphere & Africa
    # left. Each label is flanked by its two phase numbers, aligned with the label as one unit
    # and counting counterclockwise (phase 1 lower-left; octant->phase map is makePhases in
    # utils/indexTools.jl). Top/bottom labels stay horizontal; left/right run vertically.
    label_color = :gray30
    number_color = :gray40
    noff = 0.32 * xr        # centre-to-number gap along the label (xr == yr for the square box)
    ytop = ylimits[2] - 0.04 * yr
    ybot = ylimits[1] + 0.04 * yr
    xright = xlimits[2] - 0.05 * xr
    xleft = xlimits[1] + 0.05 * xr
    number(x, y, s) = CairoMakie.text!(ax, x, y; text = s, align = (:center, :center),
        fontsize = label_fontsize, color = number_color)
    # Top:  7  Western Pacific  6
    CairoMakie.text!(ax, xc, ytop; text = region_labels.west_pacific,
        align = (:center, :top), fontsize = label_fontsize, color = label_color)
    number(xc - noff, ytop - 0.02 * yr, "7"); number(xc + noff, ytop - 0.02 * yr, "6")
    # Bottom:  2  Indian Ocean  3
    CairoMakie.text!(ax, xc, ybot; text = region_labels.indian_ocean,
        align = (:center, :bottom), fontsize = label_fontsize, color = label_color)
    number(xc - noff, ybot + 0.02 * yr, "2"); number(xc + noff, ybot + 0.02 * yr, "3")
    # Right (vertical):  5 (top) / Maritime Continent / 4 (bottom)
    CairoMakie.text!(ax, xright, yc; text = region_labels.maritime_continent,
        align = (:center, :center), rotation = -π / 2, fontsize = label_fontsize, color = label_color)
    number(xright, yc + noff, "5"); number(xright, yc - noff, "4")
    # Left (vertical):  8 (top) / West. Hem. & Africa / 1 (bottom)
    CairoMakie.text!(ax, xleft, yc; text = region_labels.west_hem_africa,
        align = (:center, :center), rotation = π / 2, fontsize = label_fontsize, color = label_color)
    number(xleft, yc + noff, "8"); number(xleft, yc - noff, "1")
    return nothing
end

# Shared drawing preserves the original labeled RMM panel, including nearest-date markers.
function _draw_labeled_index!(ax, series_vals, dates_plot, xlimits, ylimits;
    mark_dates, region_labels, region_label_fontsize = 10, mark_label_fontsize = 10)
    time_color = collect(eachindex(dates_plot))
    colorrange = extrema(time_color)
    _draw_phase_regions!(ax, xlimits, ylimits, region_labels;
        label_fontsize = region_label_fontsize)
    CairoMakie.lines!(ax, real.(series_vals), imag.(series_vals);
        color = time_color, colormap = :hawaii, colorrange = colorrange)
    for (mi, md) in enumerate(mark_dates)
        idx = argmin(abs.(Dates.value.(dates_plot .- md)))
        CairoMakie.scatter!(ax, [real(series_vals[idx])], [imag(series_vals[idx])];
            color = :black, marker = :star5, markersize = 12)
        # Stagger the labels above/below: on the Koopman panel the two dates sit closer
        # together than the labels are wide, so same-height labels run into each other.
        above = isodd(mi)
        CairoMakie.text!(ax, real(series_vals[idx]), imag(series_vals[idx]);
            text = Dates.format(md, "d u"),
            align = (:center, above ? :bottom : :top),
            fontsize = mark_label_fontsize, color = :black,
            offset = (0.0f0, above ? 7.0f0 : -7.0f0))
    end
    CairoMakie.xlims!(ax, xlimits...)
    CairoMakie.ylims!(ax, ylimits...)
    return nothing
end

# Additional labeled figure: the same complex-plane comparison with the standard MJO regional labels,
# grey phase-octant lines, and markers at the requested dates (default 25 Oct / 16 Dec 2009).
function _plot_rmm_index_labeled!(paths; output_dir, figure_width, title,
    output_group::Symbol = :paper,
    filename::AbstractString = "RMM_vs_eigenfunction_YOTC_complex_labeled.png",
    width_scale = 2 / 3,
    mark_dates::Vector{Date} = [Date(2009, 10, 25), Date(2009, 12, 16)],
    region_labels = (indian_ocean = "Indian Ocean", west_hem_africa = "West. Hem. & Africa",
        west_pacific = "Western Pacific", maritime_continent = "Maritime Continent"),
    region_label_fontsize = 10,   # slightly under the 12pt axis text
    mark_label_fontsize = 10)
    series = _yotc_complex_series()
    dates_plot = series.dates_plot
    rmm_plot = series.rmm_plot
    zeta_plot = series.zeta_plot
    xlimits, ylimits = _square_complex_limits(rmm_plot, zeta_plot)
    # Enlarge the plotted box ~20% about its centre so the region labels and phase numbers clear
    # the data spiral. Local to this figure (the shared _square_complex_limits is untouched) and
    # the figure size below is unchanged, so total figure width stays constant — the data just
    # renders smaller within the same frame.
    xmid = (xlimits[1] + xlimits[2]) / 2
    ymid = (ylimits[1] + ylimits[2]) / 2
    xhalf = (xlimits[2] - xlimits[1]) / 2 * 1.2
    yhalf = (ylimits[2] - ylimits[1]) / 2 * 1.2
    xlimits = (xmid - xhalf, xmid + xhalf)
    ylimits = (ymid - yhalf, ymid + yhalf)

    paper_width = round(Int, figure_width * width_scale)
    height = round(Int, paper_width * 0.66)
    fig = CairoMakie.Figure(size = (paper_width, height), fontsize = 12)
    CairoMakie.Label(fig[0, 1:2], title; fontsize = 14, font = :bold)
    ax1 = CairoMakie.Axis(fig[1, 1]; xlabel = "RMM1", ylabel = "RMM2",
        title = "RMM index", aspect = CairoMakie.DataAspect(),
        xgridvisible = false, ygridvisible = false)
    ax2 = CairoMakie.Axis(fig[1, 2]; xlabel = "Re(ζ)", ylabel = "Im(ζ)",
        title = "Koopman eigenfunction", aspect = CairoMakie.DataAspect(),
        xgridvisible = false, ygridvisible = false)

    for (ax, series_vals) in ((ax1, rmm_plot), (ax2, zeta_plot))
        _draw_labeled_index!(ax, series_vals, dates_plot, xlimits, ylimits;
            mark_dates, region_labels, region_label_fontsize, mark_label_fontsize)
    end
    CairoMakie.colgap!(fig.layout, 6)
    mkpath(output_dir)
    path = joinpath(output_dir, filename)
    CairoMakie.save(path, fig; px_per_unit = 5)
    push!(paths[output_group], path)
    return nothing
end

"""
    _plot_mjo_indices_labeled!(paths; output_dir, figure_width, title, ...)

Save the supplementary RMM/Koopman/OMI/FMO phase comparison. After including
`plotting/all_paper_plotting.jl`, generate only this figure with:

    paths = Dict(:supplement => String[])
    _plot_mjo_indices_labeled!(paths;
        output_dir = joinpath(DEFAULT_FIGURE_OUTPUT_ROOT, "supplement"),
        figure_width = PAPER_FIGURE_WIDTH, title = DEFAULT_PLOT_TITLES.mjo_indices_labeled)

Uses the end-aligned, RMM-rotated eigenfunction without fitting to OMI or FMO. NOAA
observations are joined by date, with no interpolation, filtering, or rescaling.
"""
function _plot_mjo_indices_labeled!(paths; output_dir, figure_width, title,
    mark_dates::Vector{Date} = [Date(2009, 10, 25), Date(2009, 12, 16)],
    region_labels = (indian_ocean = "Indian Ocean", west_hem_africa = "West. Hem. & Africa",
        west_pacific = "Western Pacific", maritime_continent = "Maritime Continent"),
    region_label_fontsize = 10, mark_label_fontsize = 10)
    figure_width > 0 || throw(ArgumentError("figure_width must be positive"))
    series = _yotc_complex_series()
    dates_plot = series.dates_plot
    dates_plot == collect(series.plot_start:Day(1):series.plot_end) ||
        error("Expected complete daily coverage of the YOTC period")
    all(d -> d in dates_plot, mark_dates) || error("Marker dates must be present in the plotted period")
    omi_plot = _load_noaa_mjo(:omi, dates_plot)
    fmo_plot = _load_noaa_mjo(:fmo, dates_plot)
    panels = (
        (series.rmm_plot, "RMM index", "RMM1", "RMM2"),
        (series.zeta_plot, "Koopman eigenfunction", "Re(ζ)", "Im(ζ)"),
        (omi_plot, "OMI index", "OMI PC2", "−OMI PC1"),
        (fmo_plot, "FMO index", "FMO PC2", "−FMO PC1"),
    )
    xlimits, ylimits = _square_complex_limits(
        vcat(series.rmm_plot, series.zeta_plot), vcat(omi_plot, fmo_plot))
    # Same 20% padding as the original labeled comparison, shared by all four panels.
    xmid, ymid = sum(xlimits) / 2, sum(ylimits) / 2
    xhalf = (xlimits[2] - xlimits[1]) / 2 * 1.2
    yhalf = (ylimits[2] - ylimits[1]) / 2 * 1.2
    xlimits = (xmid - xhalf, xmid + xhalf)
    ylimits = (ymid - yhalf, ymid + yhalf)

    fig = CairoMakie.Figure(size = (figure_width, round(Int, figure_width * 1.12)), fontsize = 12)
    CairoMakie.Label(fig[0, 1:2], title; fontsize = 14, font = :bold)
    for (i, (values, panel_title, xlabel, ylabel)) in enumerate(panels)
        row, col = div(i - 1, 2) + 1, mod(i - 1, 2) + 1
        ax = CairoMakie.Axis(fig[row, col]; xlabel, ylabel, title = panel_title,
            aspect = CairoMakie.DataAspect(), xgridvisible = false, ygridvisible = false)
        _draw_labeled_index!(ax, values, dates_plot, xlimits, ylimits;
            mark_dates, region_labels, region_label_fontsize, mark_label_fontsize)
    end
    CairoMakie.colgap!(fig.layout, 6)
    # The title occupies the first layout row, so the panel gap is the second gap.
    CairoMakie.rowgap!(fig.layout, 2, 6)
    mkpath(output_dir)
    path = joinpath(output_dir, "RMM_OMI_FMO_vs_eigenfunction_YOTC_complex_labeled.png")
    CairoMakie.save(path, fig; px_per_unit = 5)
    push!(paths[:supplement], path)
    return nothing
end
