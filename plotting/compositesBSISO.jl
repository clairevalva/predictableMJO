function _plot_bsiso_composites!(paths; paper_dir, supplement_dir, figure_width, title,
    include_large_composites = false)
    wind = _compute_index_composites(; data_name = "u", level_index = 2, season = "JJA",
        eigenfunction_index = 12, reference_kind = :kikuchi)
    mkpath(paper_dir)
    push!(paths[:paper], _render_index_composites(wind, [2, 6],
        joinpath(paper_dir, "BSISO_composites_phases=2-6_alignment=end.png"); figure_width, title))
    if include_large_composites
        mkpath(supplement_dir)
        _save_large_variable_composites!(paths, wind; mode_name = "BSISO", supplement_dir, figure_width, title)
        olr = _compute_index_composites(; data_name = "olr", level_index = 1, season = "JJA",
            eigenfunction_index = 12, reference_kind = :kikuchi)
        _save_large_variable_composites!(paths, olr; mode_name = "BSISO", supplement_dir, figure_width, title)
    end
    return nothing
end
