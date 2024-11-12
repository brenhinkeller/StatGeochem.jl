## --- Test PerpleX
if Sys.isunix()

    # Choose perpleX version
    # perplexversion = "perplex-6.8.7"
    # perplexversion= "perplex-stable"

    # Construct file path
    # perplexdir = joinpath(resourcepath, perplexversion)
    scratchdir = "./"

    # if Sys.islinux()
    #     # Download precompiled executable
    #     if !isfile(joinpath(perplexdir,"vertex"))
    #         @info "Downloading PerpleX to $perplexdir"
    #         run(`mkdir -p $perplexdir`)
    #         file = Downloads.download("https://storage.googleapis.com/statgeochem/$perplexversion-linux.tar.gz",joinpath(resourcepath,"$perplexversion-linux.tar.gz"))
    #         run(`tar -xzf $file -C $perplexdir`)
    #     end
    # else
    #     # Compile from source
    #     if !isfile(joinpath(perplexdir,"vertex"))
    #         # Check if there is a fortran compiler
    #         run(`gfortran -v`)

    #         # Download Perplex v6.8.7 -- known to work with interface used here
    #         file = Downloads.download("https://storage.googleapis.com/statgeochem/$perplexversion.zip", joinpath(resourcepath,"$perplexversion.zip"))

    #         # # For a more updated perplex version, you might also try
    #         # file = download("https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip", joinpath(resourcepath,"perplex-stable.zip"))

    #         run(`unzip -u $file -d $resourcepath`) # Extract
    #         system("cd $perplexdir; make") # Compile
    #     end
    # end

    # Kelemen (2014) primitive continental basalt excluding Mn and Ti since most melt models can"t handle them..
    elements =    [ "SIO2", "AL2O3",  "FEO",  "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    composition = [50.0956, 15.3224, 8.5103, 9.2520, 9.6912, 2.5472, 0.8588, 2.0000, 0.6000,]

    # Emphasis on phases from Holland and Powell -- all phases can be used with hp02ver.dat.
    # HP_solution_phases = "Omph(HP)\nOpx(HP)\nGlTrTsPg\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nfeldspar_B\nMica(CF)\nBio(TCC)\nChl(HP)\nCtd(HP)\nSapp(HP)\nSt(HP)\nIlHm(A)\nDo(HP)\nT\nB\nF\n"
    HP_solution_phases = "Omph(HP)\nOpx(HP)\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nfeldspar_B\nMica(CF)\nBio(TCC)\nCtd(HP)\nSt(HP)\nDo(HP)\nT\nB\nF\n"
    HP_excludes = ""

    ## --- # # # # # # # # # # # # # Isobaric example # # # # # # # # # # # # # # # #

    # Input parameters
    P = 1000 # Pressure, bar
    T_range = (0+273.15, 1500+273.15) # Temperature range, Kelvin

    # Configure (run build and vertex)
    melt_model = "melt(HP)"
    @time perplex_configure_isobar(scratchdir, composition, elements, P, T_range,
        dataset="hp02ver.dat",
        npoints=100,
        excludes=HP_excludes,
        solution_phases=melt_model*"\n"*HP_solution_phases
    )

    ## --- Query all properties at a single temperature -- results returned as text

    T = 850+273.15
    data_isobaric = perplex_query_point(scratchdir, T)
    system("ls ./out1/")
    # system("cat ./out1/werami.log")

    @test isfile("./out1/1_1.txt")
    @test isa(data_isobaric, String)

    ## --- Query the full isobar -- results returned as dict

    bulk = perplex_query_system(scratchdir, importas=:Tuple)
    @test isa(bulk, NamedTuple)
    @test haskey(bulk, :SIO2)
    if haskey(bulk, :SIO2)
        print("bulk.SIO2: ")
        println(bulk.SIO2)
        @test haskey(bulk, :SIO2) && all(isapprox.(bulk.SIO2, 50.66433039859823, atol=0.1))
        @test !isempty(bulk.SIO2)
    end

    melt = perplex_query_phase(scratchdir, melt_model, importas=:Tuple, manual_grid=true, npoints=100)
    @test isa(melt, NamedTuple)
    # println(melt)
    @test haskey(melt, :SIO2)

    if haskey(melt, :SIO2)
        print("melt.SiO2: ")
        println(melt.SIO2)
        @test !isempty(melt.SIO2) && !any(x->x<45, melt.SIO2) && !any(x->x>75, melt.SIO2)
        @test isapprox(melt.SIO2,  [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 66.23928873932091, 66.20069536595133, 66.1129689269046, 65.96817295304906, 65.77167961077932, 65.5254393152636, 65.23386085968349, 64.87278962035367, 64.45898710820259, 64.04463202231602, 63.2097683951158, 62.229362662382414, 61.2299, 60.24657590136964, 59.299682210095334, 58.39293503576102, 57.528805752880565, 56.697199999999995, 55.900244720195786, 55.151072424463784, 54.444161889086686, 53.77171075434216, 53.13486811907912, 52.53058424082473, 51.97033118219873, 51.615110323022066, 51.531989693602064, 51.49388455183463, 51.45425369117167, 51.4165640084052, 51.38137430931284, 51.34737946104821, 51.31412565706284, 51.282015384604605, 51.252599999999994, 51.2249358574551, 51.196284641114616, 51.169935818955075, 51.1451744274128, 51.11510000000001, 50.45745550320106, 50.40943024565815], nans=true, atol = 0.1)
    end

    modes = perplex_query_modes(scratchdir; importas=:Dict, manual_grid=true, npoints=100)
    @test isa(modes, Dict)
    @test haskey(modes, "Do(HP)")
    if haskey(modes, "Do(HP)")
        print("modes[\"Do(HP)\"]: ")
        println(modes["Do(HP)"])
        @test isapprox(modes["Do(HP)"],  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.381821359304602, 1.375915226299608, 1.3738211840050927, 1.374457738329655, 1.3705253254297285, 1.366661347624855, 1.3645275034904318, 1.3701797192769967, 1.375450810802064, 1.3797663596440084, 1.38170983102423, 1.3835908068064018, 1.3914280845348614, 1.4096547846120548, 1.4188273299579084, 1.417632212772262, 1.416519868920687, 1.415345984868765, 1.415345984868765, 1.4144863360848499, 1.4126771598516183, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true, atol = 0.1)
        # @test isapprox(modes["Omph(HP)"],[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1.5736, 7.33183, 13.3273, 13.874, 13.8044, 13.7504, 13.6605, 13.6055, 13.2465, 12.8556, 12.8012, 12.909, 12.8774, 12.8621, 12.8379, 12.8239, 12.8205, 12.839, 12.8654, 12.8914, 12.9423, 13.0084, 13.1195, 13.2487, 13.391, 13.5401, 13.7082, 13.9396, 14.1879, 14.4729, 14.754, 15.0912, 15.5081, 15.9689, 16.4671, 17.0194, 17.5064, 17.1991, 16.9685, 16.6926, 16.4602, 16.1634, 15.921, 15.659, 15.4497, 15.2485, 15.0301, 14.8809, 14.6926, 15.0711, 9.19562, NaN, NaN], nans=true)
    end

    # --- # # # # # # # # # # # Geothermal gradient example # # # # # # # # # # # #

    # Input parameters
    P_range = (280, 28000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_surf = 273.15 # Temperature of surface (K)
    geotherm = 0.01 # Geothermal gradient of 0.1 K/bar == about 28.4 K/km
    melt_model = ""

    # Configure (run build and vertex)
    @time perplex_configure_geotherm(scratchdir, composition, elements, P_range, T_surf, geotherm;
        dataset="hp02ver.dat",
        excludes=HP_excludes,
        solution_phases=HP_solution_phases,
        npoints=200,
        index=2
    )

    seismic = perplex_query_seismic(scratchdir, index=2, manual_grid=true, npoints=100)
    @test isa(seismic, Dict)
    @test haskey(seismic, "T(K)")
    @test isa(seismic["T(K)"], Vector{Float64})

    print("seismic[\"T(K)\"]: ")
    println(seismic["T(K)"])

    ## --- # # # # # # # # # # # P–T path example # # # # # # # # # # # #

    # Input parameters
    T_range = (550+273.15, 1050+273.15) #K
    PTdir = ""
    PTfilename = ""

    @time perplex_configure_path(scratchdir, composition, PTdir, PTfilename, elements, T_range, 
        dataset = "hp02ver.dat", index=1, solution_phases=HP_solution_phases, excludes=HP_excludes)
    
    modes = perplex_query_modes(scratchdir, index=1, manual_grid=true, npoints=100)
    # print(modes)
    @test isa(modes, Dict)
    @test haskey(modes,"node")
    @test haskey(modes, "O(HP)")
    if haskey(modes, "O(HP)")
        print("modes[\"O(HP)\"]: ")
        println(modes["O(HP)"])
    end

    # --- # # # # # # # # # # # Pseudosection example # # # # # # # # # # # # #

    P_range = (1000, 5000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_range = (400+273.15, 600+273.15) # Temperature range to explore, K
    melt_model = ""

    solution_phases = "Opx(HP)\nO(HP)\nF\n"
    excludes = ""

    # Configure (run build and vertex)
    @time perplex_configure_pseudosection(scratchdir, composition,
        elements, P_range, T_range, dataset="hp02ver.dat", excludes=excludes,
        solution_phases=melt_model*solution_phases, index=1, xnodes=50, ynodes=50)

    # Query modes on diagonal line across P-T space
    modes = perplex_query_modes(scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(modes, Dict) && !isempty(modes)
    @test haskey(modes,"T(K)") && all(extrema(modes["T(K)"]) .≈ T_range)

    phase = perplex_query_phase(scratchdir, "Opx(HP)", P_range, T_range, index=1, npoints=200)
    @test isa(phase, Dict) && !isempty(phase)
    @test haskey(phase,"T(K)") && all(extrema(phase["T(K)"]) .≈ T_range)

    sys = perplex_query_system(scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(sys, Dict) && !isempty(sys)
    @test haskey(sys,"T(K)") && all(extrema(sys["T(K)"]) .≈ T_range)

    # Query seismic properties on diagonal line across P-T space
    seismic = perplex_query_seismic(scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(seismic, Dict) && !isempty(seismic)
    @test haskey(seismic,"T(K)") && !isempty(seismic["T(K)"]) && all(extrema(seismic["T(K)"]) .≈ T_range)
    @test haskey(seismic, "rho,kg/m3") && !isempty(seismic["rho,kg/m3"]) && !any(x->x<2700, seismic["rho,kg/m3"]) && !any(x->x>3200, seismic["rho,kg/m3"])


    # Query properties on a manually-specified diagonal P-T line
    P = range(first(P_range), last(P_range), length=16)
    T = range(first(T_range), last(T_range), length=16)

    modes = perplex_query_modes(scratchdir, P, T, index=1)
    @test isa(modes, Dict) && !isempty(modes)
    @test haskey(modes,"T(K)") && all(isapprox.(modes["T(K)"], T, atol=0.1))

    phase = perplex_query_phase(scratchdir, "Opx(HP)", P, T, index=1)
    @test isa(phase, Dict) && !isempty(phase)
    @test haskey(phase,"T(K)") && all(isapprox.(phase["T(K)"], T, atol=0.1))

    sys = perplex_query_system(scratchdir, P, T, index=1)
    @test isa(sys, Dict) && !isempty(sys)
    @test haskey(sys,"T(K)") && all(isapprox.(sys["T(K)"], T, atol=0.1))
    @test haskey(sys, "rho,kg/m3") && !any(x->x<2700, sys["rho,kg/m3"]) && !any(x->x>3200, sys["rho,kg/m3"])

    seismic = perplex_query_seismic(scratchdir, P, T, index=1)
    @test isa(seismic, Dict) && !isempty(seismic)
    @test haskey(seismic,"T(K)")
    if haskey(seismic,"T(K)")
        @test isapprox(seismic["T(K)"], T, atol=0.01)
        print("T (K):\t")
        println(seismic["T(K)"])
    end
    @test haskey(seismic, "rho,kg/m3")
    if haskey(seismic, "rho,kg/m3")
        @test !any(x->x<2700, seismic["rho,kg/m3"]) && !any(x->x>3200, seismic["rho,kg/m3"])
        print("rho,kg/m3:\t")
        println(seismic["rho,kg/m3"])
    end

end

