# Large logical field dimensions use unallocated HDF5 storage; fixtures stay small.
function climate_fixture(path; units="J m**-2", time_units="days since 2000-01-01", first_day=0)
    h5open(path,"w") do f
        create_dataset(f,"ttr",Float32,(240,121,9132))
        HDF5.attributes(f["ttr"])["units"] = units
        f["valid_time"] = collect(first_day:(first_day+9131))
        HDF5.attributes(f["valid_time"])["units"] = time_units
        f["longitude"] = collect(0:1.5:358.5)
        HDF5.attributes(f["longitude"])["units"] = "degrees_east"
        f["latitude"] = collect(-90:1.5:90)
        HDF5.attributes(f["latitude"])["units"] = "degrees_north"
    end
end

@testset "prepared format and optional snapshot identity" begin
    mktempdir() do dir
        path = joinpath(dir,"field.h5")
        climate_fixture(path)
        @test filesize(path) < 200_000
        @test validate_climate_input(path,"olr")
        record = Dict("bytes"=>filesize(path),"sha256"=>bytes2hex(open(SHA.sha256,path)))
        @test check_input_snapshot(path,record;checksums=true)
        climate_fixture(path;units="J/m^2",time_units="days since 2000-01-01T00:00:00Z")
        @test check_input_snapshot(path,record)
        @test validate_climate_input(path,"olr")
        @test_throws ErrorException check_input_snapshot(path,record;checksums=true)
        # Isolate hash rejection from size rejection.
        same_size_wrong_hash = Dict("bytes"=>filesize(path),"sha256"=>repeat("0",64))
        @test_throws ErrorException check_input_snapshot(path,same_size_wrong_hash;checksums=true)
        climate_fixture(path;units="W m^-2")
        @test_throws ErrorException validate_climate_input(path,"olr")
        climate_fixture(path;time_units="hours since 2000-01-01")
        @test_throws ErrorException validate_climate_input(path,"olr")
        climate_fixture(path;time_units="days since 2001-01-01")
        @test_throws ErrorException validate_climate_input(path,"olr")
        climate_fixture(path;first_day=1)
        @test_throws ErrorException validate_climate_input(path,"olr")
    end
end

@testset "physical unit spellings and forecast time conversion" begin
    mktempdir() do dir
        for (units,allowed) in (("J m⁻²",("jm-2","j/m2")),("m s^-1",("ms-1","m/s")),
            ("m/s",("ms-1","m/s")),("metres",("m","metres")),("hPa",("hpa","mbar")))
            h5open(joinpath(dir,"units.h5"),"w") do f
                f["value"] = [1.0]
                HDF5.attributes(f["value"])["units"] = units
                @test isnothing(check_input_units(f["value"],allowed))
                @test_throws ErrorException check_input_units(f["value"],("incompatible",))
            end
        end
    end
    @test convert_valid_time([0,1],"days since 2020-01-01") == [DateTime(2020,1,1),DateTime(2020,1,2)]
    @test parse_units_datetime("days since 2020-01-01 00:00:00 UTC") == DateTime(2020,1,1)
    @test_throws ErrorException parse_units_datetime("seconds since 2020-01-01")
end

@testset "supplied NOAA indices keep their dates and orientation" begin
    mktempdir() do dir
        dates = [Date(2009,10,1),Date(2009,10,2)]
        omi = joinpath(dir,"omi.1x.txt")
        fmo = joinpath(dir,"fmo.1x.txt")
        write(omi,"2009 10 1 3.00000 4.00000 5.00000\n2009 10 2 0.00000 2.00000 2.00000\n")
        write(fmo,"2009 10 1 0 3.00000 4.00000 5.00000\n2009 10 2 0 0.00000 2.00000 2.00000\n")
        for index in (:omi,:fmo)
            @test _load_noaa_mjo(index,dates;data_dir=dir) == [4-3im,2+0im]
            @test_throws ErrorException _load_noaa_mjo(index,[Date(2009,10,3)];data_dir=dir)
        end
        for text in ("2009 10 1 3 4\n", "2009 10 1 3 4 1\n", "2009 10 1 NaN 4 5\n",
            "2009 10 1 3 4 5\n2009 10 1 3 4 5\n", "2009 10 2 3 4 5\n2009 10 1 3 4 5\n")
            write(omi,text)
            @test_throws ErrorException _read_noaa_mjo(omi,:omi)
        end
        write(fmo,"2009 10 1 12 3 4 5\n")
        @test_throws ErrorException _read_noaa_mjo(fmo,:fmo)
    end
end
