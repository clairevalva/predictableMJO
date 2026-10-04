# Exercise the actual computation drivers on a small deterministic prepared input.
# No original repository or research dataset is needed for this check.
module SmallNLSA
include(joinpath(@__DIR__, "..", "doNLSA.jl"))
end
module SmallKoopman
include(joinpath(@__DIR__, "..", "doKoopman.jl"))
end
module SmallReconstruction
include(joinpath(@__DIR__, "..", "doReconstruction.jl"))
end
module SmallOSE
include(joinpath(@__DIR__, "..", "computeOSE.jl"))
end
module RunnerSmoke
include(joinpath(@__DIR__, "..", "run_paper.jl"))
end

@testset "portable computation pipeline" begin
    # The runner loads configuration before utilities; plotting loads them in the reverse order.
    @test isdefined(RunnerSmoke, :data_details)
    mktempdir() do dir
        settings = merge(PAPER, (train_last = 80, nemb = 8, neighbors = 12, ose_neighbors = 16,
            nphi = 11, diffusion_modes = 21, ose_start = Date(2000,3,21), ose_end = Date(2000,3,25)))
        feature_file = joinpath(dir, "features.h5")
        features = [sin(0.11*t + v*l) + 0.2*cos(0.031*t*t + l) for v in 1:3, l in 1:4, t in 1:100]
        h5open(feature_file, "w") do f
            f["olr_u"] = features
        end
        nlsa_file, koopman_file = joinpath(dir,"nlsa.h5"), joinpath(dir,"koopman.h5")
        SmallNLSA.compute_nlsa(; output=nlsa_file, feature_file, settings)
        @test h5open(f -> Set(keys(f)), nlsa_file, "r") == Set(["W", "bw", "phi", "kappa"])
        @test size(h5read(nlsa_file,"phi")) == (72,21)
        SmallKoopman.compute_koopman(; input=nlsa_file, output=koopman_file, settings)
        zeta = h5read(koopman_file,"zetaL")
        @test size(zeta) == (72,11)
        @test all(isfinite,zeta)
        rec_file = joinpath(dir,"reconstruction.h5")
        SmallReconstruction.compute_reconstruction("RMMdata"; output=rec_file,settings,feature_file,koopman_file)
        @test size(h5read(rec_file,"recons"),3) == 79
        @test all(isfinite,h5read(rec_file,"Acoeff"))
        ose_file = joinpath(dir,"ose.h5")
        SmallOSE.compute_ose(; output=ose_file,settings,feature_file,nlsa_file,koopman_file)
        @test get_OSE_dates(ose_file) == collect(settings.ose_start:Day(1):settings.ose_end)
        phi = h5read(ose_file,"phi_vals")
        @test size(phi) == (11,5)
        @test all(isfinite,phi)
        @test size(h5read(ose_file,"zeta_vals")) == (11,5)
        @test h5read(ose_file,"date_convention") == "cone_transition_end"
        @test_throws ErrorException SmallKoopman.compute_koopman(; input=nlsa_file,output=koopman_file,settings)
    end
end
