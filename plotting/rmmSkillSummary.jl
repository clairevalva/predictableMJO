function _plot_rmm_skill!(paths; paper_dir,
    figure_width, title,
    paper_modes::AbstractString = "j7",
    paper_season::AbstractString = "DJFMA",
    mode_tags = ["j7"],
    seasons = ["DJFMA"],
    skill_file = paper_cache("rmm_pred_skill_2020-12-01_2024-04-30.h5"),
    curve_labels = Dict{String, String}())
    isfile(skill_file) || error("Missing RMM skill cache $(skill_file); run " *
        "model_comparisons/rmm_prediction_skill.jl to build it")

    # editable styling: S2S models solid, Koopman curves dashed, zeta forecast dotted
    palette = Makie.wong_colors()
    models = [ # (cache key, legend label, color)
        ("s2s_ecmwf", "ECMWF", palette[1]),
        ("s2s_hmcr", "HMCR", palette[2]),
        ("s2s_kma", "KMA", palette[3]),
        ("s2s_ncep", "NCEP", palette[4]),
        ("s2s_ukmo", "UKMO", palette[5]),
    ]
    mode_labels = Dict(
        "j7" => "Indirect RMM KAF",
        "j12_15_7" => "RMM from Koopman (j=12,15,7)",
    )
    color_kaf = palette[6]
    # not palette[7]: wong's 7th is yellow, which is illegible as a line on white
    color_koopman_rmm = CairoMakie.to_color(:darkviolet)
    skill_threshold = 0.5

    seasons_file = h5read(skill_file, "meta/seasons")

    for mode_tag in mode_tags, season in seasons
        seas_j = findfirst(==(season), seasons_file)
        isnothing(seas_j) && error("season $(season) not found in $(skill_file)")

        # Solid = operational S2S models, dotted = ours.
        series = vcat(
            [(key, label, color, :solid) for (key, label, color) in models],
            [("koopman_zeta_j7", "Koopman MJO KAF", CairoMakie.to_color(:black), :dot)],
            [("kaf_rmm", "Direct RMM KAF", color_kaf, :dot)],
            [("koopman_rmm_$(mode_tag)", get(mode_labels, mode_tag, mode_tag),
                color_koopman_rmm, :dot)],
        )
        series = [(key, get(curve_labels, key, label), color, linestyle)
            for (key, label, color, linestyle) in series]

        height = round(Int, figure_width * 0.7)
        fig = CairoMakie.Figure(size = (figure_width, height), fontsize = 12)
        CairoMakie.Label(fig[0, 1], title * " (" * season * ")";
            fontsize = 14, font = :bold, tellwidth = false)
        ax = CairoMakie.Axis(fig[1, 1];
            xlabel = "Lead time (days)",
            ylabel = "Bivariate correlation",
            xticks = 0:10:50,
            yticks = 0:0.2:1,
        )

        for (key, label, color, linestyle) in series
            lags = h5read(skill_file, "skill/$(key)/lags")
            corr = h5read(skill_file, "skill/$(key)/corr")[:, seas_j]
            CairoMakie.scatterlines!(ax, lags, corr;
                color = color, linestyle = linestyle, linewidth = 1.5, markersize = 5,
                label = label)
        end
        CairoMakie.hlines!(ax, [skill_threshold];
            color = :red, linestyle = :dash, linewidth = 1.5)
        # Curves that go negative at long lead are clipped rather than drawn.
        CairoMakie.xlims!(ax, 0, 50)
        CairoMakie.ylims!(ax, 0, 1)

        # Two columns via nbanks, which fills row-major -- so interleave the groups row by row
        # to keep our Koopman curves in the left column and the operational S2S models in the
        # right. The short (Koopman) column is padded with invisible entries. Translucent
        # because with eight entries no in-panel corner is free: the Koopman KAF tail (DJFMA)
        # and the long-lead S2S curves (JJASO) run under the box.
        operational = series[1:length(models)]
        ours = series[(length(models) + 1):end]
        n_rows = max(length(operational), length(ours))
        blank = (nothing, "", (:white, 0.0), :solid)
        entry(s, i) = i <= length(s) ? s[i] : blank
        elems = CairoMakie.LineElement[]
        labels = String[]
        for i in 1:n_rows, group in (ours, operational)
            _, label, color, linestyle = entry(group, i)
            push!(elems, CairoMakie.LineElement(color = color, linestyle = linestyle,
                linewidth = 2))
            push!(labels, label)
        end
        CairoMakie.axislegend(ax, elems, labels;
            position = :rt, labelsize = 12, framevisible = true, nbanks = 2,
            bgcolor = (:white, 0.6), rowgap = 0, colgap = 6,
            patchsize = (12.0f0, 2.0f0), patchlabelgap = 3,
            padding = (5.0f0, 5.0f0, 3.0f0, 3.0f0))

        outdir = paper_dir
        mkpath(outdir)
        path = joinpath(outdir, "RMM_skill_modes=$(mode_tag)_season=$(season).png")
        CairoMakie.save(path, fig; px_per_unit = 3)
        push!(paths[:paper], path)

    end
    return nothing
end
