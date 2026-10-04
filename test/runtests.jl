using Test
include(joinpath(@__DIR__, "..", "plotting", "all_paper_plotting.jl"))
include(joinpath(@__DIR__, "..", "analogForecasting.jl"))
include(joinpath(@__DIR__, "..", "validate_inputs.jl"))
include("check_conventions.jl")
include("check_inputs.jl")
include("check_pipeline.jl")
