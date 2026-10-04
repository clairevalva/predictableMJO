using Dates, HDF5, Statistics, LinearAlgebra, StatsBase
if !isdefined(@__MODULE__, :PAPER)
    include(joinpath(@__DIR__, "..", "paper_config.jl"))
end
include("dataTools.jl")
include("indexTools.jl")
include("bandpass.jl")
include("prediction.jl")
