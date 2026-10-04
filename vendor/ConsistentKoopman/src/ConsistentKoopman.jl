module ConsistentKoopman

using LinearAlgebra
using SparseArrays
using Distances
using ForwardDiff

include("kernels.jl")
include("utils.jl")
include("modelcomponents.jl")
include("kernelsMatrix.jl")
include("fdGenerator.jl")
include("OSE.jl")
include("reconstruct.jl")

end
