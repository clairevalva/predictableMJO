function _plot_eigenperiod!(paths; paper_dir, figure_width, title)
    expname = "init"
    repo_dir = normpath(joinpath(@__DIR__, ".."))
    fn_koop = paper_cache("$(expname)_Keigen.h5")
    fn_rec = paper_cache("init_recon_RMMdata.h5")

    omegas, _, zetaL, _, periods = getK(fn_koop)
    good_js = [7, 12]
    bad_js = [15]
    marked_js = vcat(good_js, bad_js)
    pos_js = (1:101)[imag.(omegas) .>= 0]
    other_js = setdiff(pos_js, marked_js)

    maxeps, _ = ConsistentKoopman.sortautocorr(
        zetaL,
        imag.(omegas),
        45,
        1.0;
        Flim = 1.0,
    )
    dirichlet_energy = real.(omegas)

    rec = h5read(fn_rec, "recons")
    rec_js = vec(h5read(fn_rec, "indices"))
    rec = reshape(rec, length(rec_js), 3, 240, :)
    variances = log.(var(reshape(2 .* real.(rec), length(rec_js), :); dims = 2)[:])
    variance_by_j = Dict(Int(j) => variances[i] for (i, j) in enumerate(rec_js))
    color_values(js) = [variance_by_j[Int(j)] for j in js]
    clim = extrema(variances)

    function render(path; mark_bad, label_modes = false)
        height = round(Int, figure_width * 250 / 475)
        fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 12)
        CairoMakie.Label(fig[0, 1:2], title; fontsize = 14, font = :bold)
        ax1 = CairoMakie.Axis(
            fig[1, 1];
            xlabel = "maximal autocorrelation error",
            ylabel = "period (days)",
        )
        ax2 = CairoMakie.Axis(fig[1, 2]; xlabel = "Dirichlet energy")

        CairoMakie.scatter!(ax1, maxeps[other_js], periods[other_js];
            color = color_values(other_js), colormap = :amp, colorrange = clim, strokewidth = 0)
        sc = CairoMakie.scatter!(ax2, -dirichlet_energy[other_js], periods[other_js];
            color = color_values(other_js), colormap = :amp, colorrange = clim, strokewidth = 0)
        for ax_data in ((ax1, maxeps), (ax2, -dirichlet_energy))
            ax, xdata = ax_data
            CairoMakie.scatter!(ax, xdata[good_js], periods[good_js];
                color = color_values(good_js), colormap = :amp, colorrange = clim,
                marker = :star5, markersize = 12, strokewidth = 0)
            CairoMakie.scatter!(ax, xdata[bad_js], periods[bad_js];
                color = color_values(bad_js), colormap = :amp, colorrange = clim,
                marker = mark_bad ? :diamond : :circle,
                markersize = mark_bad ? 9 : 6, strokewidth = 0)
            CairoMakie.hlines!(ax, [30, 90]; linestyle = :dash, color = :black)
            if label_modes
                # Flank the stars on the left; stagger the crowded pair on the right.
                for (j, label) in ((7, "M"), (12, "B"))
                    if ax === ax1
                        label_align = j == 7 ? (:left, :center) : (:right, :center)
                        label_offset = j == 7 ? (6, 0) : (-6, 0)
                    else
                        label_align = j == 7 ? (:right, :bottom) : (:left, :top)
                        label_offset = j == 7 ? (-2, 0) : (2, -3)
                    end
                    if ax === ax2 && j == 12
                        # Draw the halo first so the dashed guide cannot obscure the letter.
                        CairoMakie.text!(ax, xdata[j], periods[j]; text = label,
                            color = :white, font = :bold, fontsize = 11,
                            strokecolor = :white, strokewidth = 1,
                            align = label_align, offset = label_offset,
                            xautolimits = false, yautolimits = false)
                    end
                    CairoMakie.text!(ax, xdata[j], periods[j]; text = label,
                        color = :black, font = :bold, fontsize = 11,
                        align = label_align, offset = label_offset,
                        xautolimits = false, yautolimits = false)
                end
            end
        end

        CairoMakie.linkyaxes!(ax1, ax2)
        valid_periods = filter(isfinite, periods[vcat(other_js, marked_js)])
        isempty(valid_periods) || CairoMakie.ylims!(ax1, 0, 1.1 * maximum(valid_periods))
        CairoMakie.hideydecorations!(ax2; grid = false)
        CairoMakie.Colorbar(fig[1, 3], sc; label = "log(variance)", width = 10)
        CairoMakie.save(path, fig; px_per_unit = 5)
        return path
    end

    mkpath(paper_dir)
    paper_path = joinpath(paper_dir, "periods_comparison_AC_lags=45.png")
    push!(paths[:paper], render(paper_path; mark_bad = true))
    labeled_path = joinpath(paper_dir, "periods_comparison_AC_lags=45_labeled.png")
    push!(paths[:paper], render(labeled_path; mark_bad = true, label_modes = true))

    return nothing
end
