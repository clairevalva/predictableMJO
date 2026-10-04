@testset "public figure inventory" begin
    small = expected_figure_paths()
    large = expected_figure_paths(; include_large_composites = true)
    @test length(small) == 11
    @test length(large) == length(unique(large)) == 19
    @test issubset(Set(small), Set(large))
    @test all(p -> !occursin("alignment=center", p) && !occursin("sample_sd", p), large)
    @test "supplement/RMM_OMI_FMO_vs_eigenfunction_YOTC_complex_labeled.png" in small
    output_root = tempname()
    @test_throws ArgumentError make_all_paper_figures(; output_root, figure_width = 0)
    @test !ispath(output_root)
end

@testset "physical dates and field indexing" begin
    dates = PAPER.start_date .+ Day.((PAPER.train_first:PAPER.train_last) .- 1)
    inds = makeTakensInd_forward(64, 7241)
    @test first(dates[inds]) == Date(2000, 3, 5)
    @test last(dates[inds]) == Date(2019, 12, 31)
    field_dates = collect(PAPER.start_date:Day(1):PAPER.ose_end)
    field_inds = select_dates(field_dates, dates[inds])
    @test field_dates[field_inds] == dates[inds]
    @test field_inds == collect(inds) .+ 1
    @test_throws ErrorException select_dates([Date(2000, 1, 1)], [Date(2000, 1, 2)])
    @test_throws ErrorException select_dates(fill(Date(2000, 1, 1), 2), [Date(2000, 1, 1)])
end

@testset "time-axis filtering" begin
    time = 1:400
    signal = sin.(2pi .* time ./ 45)
    observations = repeat(reshape(signal, 1, :), 4, 1)
    filter = makeMJOFF(; days_low = 90, days_high = 30, lanN = 49)
    actual = permutedims(filter(permutedims(observations)))
    @test size(actual) == size(observations)
    @test actual[1, :] ≈ filter(signal)
    @test actual[4, :] ≈ actual[1, :]
end

@testset "HDF5 dates and reconstruction modes" begin
    mktempdir() do dir
        for dates in (["2020-01-01", "2020-01-03"], ["2020-01-01", "2020-01-02", "2020-01-03"])
            path = joinpath(dir, "dates.h5")
            h5open(path, "w") do f
                f["dates"] = dates
            end
            @test get_OSE_dates(path) == collect(Date(2020, 1, 1):Day(1):Date(2020, 1, 3))
        end
        path = joinpath(dir, "rec.h5")
        h5open(path, "w") do f
            f["indices"] = [7, 12]
            f["recons"] = reshape(ComplexF64.(1:24), 2, 3, 4)
            f["Acoeff"] = reshape(ComplexF64.(1:12), 2, 3, 2)
        end
        rec, modes, coeffs = get_reconstructed_data(path; ret_coeffs = true, load_small_j = [12])
        @test modes == [12]
        @test size(rec) == (1, 3, 4)
        @test size(coeffs) == (1, 3, 2)
        @test_throws ErrorException get_reconstructed_data(path; load_small_j = [99])
    end
end

@testset "forecasts use one issue date and all delay coordinates" begin
    nEmb, maxlag, nD, nT_forecast = 5, 3, 2, 6
    use_k, indices_rec = [2,4], [4,2]
    A = [complex(10i+d+q/10,i-d+q/7) for i in 1:2, d in 1:nD, q in 1:nEmb]
    zeta = [cis(0.17t+0.3i)+0.05*complex(t,-i) for i in 1:2, t in 1:(nT_forecast+maxlag+nEmb)]
    predictions = [zeta[i,t+lag] for i in 1:2, t in 1:nT_forecast, lag in 0:(maxlag+nEmb-1)]
    forecast = reconstruct_forecast_from_delay_coeffs(A,indices_rec,use_k,predictions,nEmb,maxlag,nD;
        zeta_indices=use_k)
    reconstructed = zeros(ComplexF64,nD,size(zeta,2)+nEmb-1)
    for (i,k) in enumerate(use_k)
        a = A[findfirst(==(k),indices_rec),:,:]
        reconstructed .+= ConsistentKoopman.reconstructfromproj(a,zeta[i,:],nEmb)
    end
    for lag in 0:maxlag, t in 1:nT_forecast
        @test forecast[:,t,lag+1] ≈ 2 .* real.(reconstructed[:,t+lag+nEmb-1])
    end
    # Changing forecasts issued later must not affect the first issue's forecast.
    changed = copy(predictions)
    changed[:,2:end,:] .= 1e6+1e6im
    changed_forecast = reconstruct_forecast_from_delay_coeffs(A,indices_rec,use_k,changed,nEmb,maxlag,nD;
        zeta_indices=use_k)
    @test changed_forecast[:,1,:] == forecast[:,1,:]
    @test_throws AssertionError reconstruct_forecast_from_delay_coeffs(A,indices_rec,use_k,
        predictions[:,:,1:maxlag+1],nEmb,maxlag,nD;zeta_indices=use_k)
end

@testset "direct forecasts target the end of the training window" begin
    nEmb, maxlag, nT = 5, 3, 12
    phis = hcat(ones(nT),cos.(0.2 .* (1:nT)))
    ose = [1.0 0.2; 1.0 0.8]
    raw = hcat((1:(nT+nEmb-1)).^2, sin.(1:(nT+nEmb-1)))
    forecast = make_direct_issue_date_forecast(phis,ose,raw,nEmb,maxlag)
    for lag in 0:maxlag
        coeff = phis[1:(end-lag),:]' * raw[(nEmb+lag):end,:] / (nT-lag)
        @test forecast[:,:,lag+1] ≈ transpose(ose*coeff)
    end
end
