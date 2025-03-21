## --- Test PerpleX
if Sys.isunix()

    scratchdir = "./"

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
    @info "perplex_configure_isobar:"
    @time perplex_configure_isobar(scratchdir, composition, elements, P, T_range,
        dataset="hp02ver.dat",
        npoints=100,
        excludes=HP_excludes,
        solution_phases=melt_model*"\n"*HP_solution_phases
    )

    ## --- Query all properties at a single temperature -- results returned as text

    T = 850+273.15
    data_isobaric = perplex_query_point(scratchdir, T)
    # system("ls ./out1/")

    @test isfile("./out1/1_1.txt")
    @test isa(data_isobaric, String)

    ## --- Query the full isobar -- results returned as dict

    bulk = perplex_query_system(scratchdir, importas=:Tuple)
    @test isa(bulk, NamedTuple)
    @test haskey(bulk, :SIO2)
    if haskey(bulk, :SIO2)
        # print("bulk.SIO2: ")
        # println(bulk.SIO2)
        @test haskey(bulk, :SIO2) && all(isapprox.(bulk.SIO2, 50.66433039859823, atol=0.1))
        @test !isempty(bulk.SIO2)
    end

    melt = perplex_query_phase(scratchdir, melt_model, importas=:Tuple, manual_grid=true, npoints=100)
    @test isa(melt, NamedTuple)
    # println(melt)
    @test haskey(melt, :SIO2)

    if haskey(melt, :SIO2)
        # print("melt.SiO2: ")
        # println(melt.SIO2)
        @test !isempty(melt.SIO2) && !any(x->x<45, melt.SIO2) && !any(x->x>75, melt.SIO2)
        @test isapprox(melt.SIO2,  [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 66.23928873932091, 66.20069536595133, 66.1129689269046, 65.96817295304906, 65.77167961077932, 65.5254393152636, 65.23386085968349, 64.87278962035367, 64.45898710820259, 64.04463202231602, 63.2097683951158, 62.229362662382414, 61.2299, 60.24657590136964, 59.299682210095334, 58.39293503576102, 57.528805752880565, 56.697199999999995, 55.900244720195786, 55.151072424463784, 54.444161889086686, 53.77171075434216, 53.13486811907912, 52.53058424082473, 51.97033118219873, 51.615110323022066, 51.531989693602064, 51.49388455183463, 51.45425369117167, 51.4165640084052, 51.38137430931284, 51.34737946104821, 51.31412565706284, 51.282015384604605, 51.252599999999994, 51.2249358574551, 51.196284641114616, 51.169935818955075, 51.1451744274128, 51.11510000000001, 50.45745550320106, 50.40943024565815], nans=true, atol = 0.1)
    end

    modes = perplex_query_modes(scratchdir; importas=:Dict, manual_grid=true, npoints=100)
    @test isa(modes, Dict)
    @test haskey(modes, "Do(HP)")
    if haskey(modes, "Do(HP)")
        # print("modes[\"Do(HP)\"]: ")
        # println(modes["Do(HP)"])
        @test isapprox(modes["Do(HP)"],  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.381821359304602, 1.375915226299608, 1.3738211840050927, 1.374457738329655, 1.3705253254297285, 1.366661347624855, 1.3645275034904318, 1.3701797192769967, 1.375450810802064, 1.3797663596440084, 1.38170983102423, 1.3835908068064018, 1.3914280845348614, 1.4096547846120548, 1.4188273299579084, 1.417632212772262, 1.416519868920687, 1.415345984868765, 1.415345984868765, 1.4144863360848499, 1.4126771598516183, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true, atol = 0.1)
    end
    if haskey(modes, "Omph(HP)")
        # print("modes[\"Omph(HP)\"]: ")
        # println(modes["Omph(HP)"])
        @test isapprox(modes["Omph(HP)"],[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0339723566900236, 14.960931973118665, 15.062263362661524, 15.164506420058737, 15.268082701516247, 15.372961492756188, 15.47935889951888, 16.023147006631838, 17.195113428726106, 18.546171502024475, 18.513626566383596, 18.921480372141023, 18.938381437801173, 18.95517180349422, 18.971868712229128, 18.987725829804052, 19.003917068704503, 19.019482276868484, 19.034156503334277, 19.04772018291329, 18.955385576325604, 18.66836547025949, 18.700498131123354, 18.672603711178727, 18.689706911851058, 18.70741316163515, 18.726321687261198, 18.74703357720072, 18.769573446275686, 18.793658820257388, 18.821264934638325, 18.852397423440955, 18.8867110753726, 18.978007255636545, 19.10394367818196, 19.245715839150208, 19.40418088373046, 19.579886180363193, 19.774964084364957, 19.991810067185163, 20.241688665240932, 20.52557287224623, 20.84732371711889, 21.210565263460555, 21.618170448738216, 22.09663695953447, 22.646817318444732, 23.25441852307975, 23.214082721363855, 22.81470048851166, 22.43291893689257, 22.066189222773463, 21.721311293034898, 21.399604352976905, 21.097623567233256, 20.8182863046844, 20.557653077344202, 20.318180516294838, 20.09688085563381, 19.885716974177978, 19.695776424220245, 19.521123834193872, 18.711834731308517, 0.0, 0.0], nans=true, atol = 0.1)
    end

    # --- # # # # # # # # # # # Geothermal gradient example # # # # # # # # # # # #

    # Input parameters
    P_range = (280, 28000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_surf = 273.15 # Temperature of surface (K)
    geotherm = 0.01 # Geothermal gradient of 0.1 K/bar == about 28.4 K/km
    melt_model = ""

    # Configure (run build and vertex)
    @info "perplex_configure_geotherm:"
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
    # print("seismic[\"T(K)\"]: ")
    # println(seismic["T(K)"])
    @test isapprox(seismic["T(K)"], [275.95, 278.75, 281.55, 284.35, 287.15, 289.95, 292.75, 295.55, 298.35, 301.15, 303.95, 306.75, 309.55, 312.35, 315.15, 317.95, 320.75, 323.55, 326.35, 329.15, 331.95, 334.75, 337.55, 340.35, 343.15, 345.95, 348.75, 351.55, 354.35, 357.15, 359.95, 362.75, 365.55, 368.35, 371.15, 373.95, 376.75, 379.55, 382.35, 385.15, 387.95, 390.75, 393.55, 396.35, 399.15, 401.95, 404.75, 407.55, 410.35, 413.15, 415.95, 418.75, 421.55, 424.35, 427.15, 429.95, 432.75, 435.55, 438.35, 441.15, 443.95, 446.75, 449.55, 452.35, 455.15, 457.95, 460.75, 463.55, 466.35, 469.15, 471.95, 474.75, 477.55, 480.35, 483.15, 485.95, 488.75, 491.55, 494.35, 497.15, 499.95, 502.75, 505.55, 508.35, 511.15, 513.95, 516.75, 519.55, 522.35, 525.15, 527.95, 530.75, 533.55, 536.35, 539.15, 541.95, 544.75, 547.55, 550.35, 553.15], nans=true, atol = 0.1)

    ## --- # # # # # # # # # # # P–T path example # # # # # # # # # # # #

    # Input parameters
    T_range = (550+273.15, 1050+273.15) #K
    PTdir = ""
    PTfilename = ""

    @info "perplex_configure_path:"
    @time perplex_configure_path(scratchdir, composition, PTdir, PTfilename, elements, T_range, 
        dataset = "hp02ver.dat", index=1, solution_phases=HP_solution_phases, excludes=HP_excludes)
    
    modes = perplex_query_modes(scratchdir, index=1, manual_grid=true, npoints=100)
    # print(modes)
    @test isa(modes, Dict)
    @test haskey(modes,"node")
    @test haskey(modes, "O(HP)")
    if haskey(modes, "O(HP)")
        # print("modes[\"O(HP)\"]: ")
        # println(modes["O(HP)"])
        @test isapprox(modes["O(HP)"], [0.13263234069934907])
    end

    # --- # # # # # # # # # # # Pseudosection example # # # # # # # # # # # # #

    P_range = (1000, 5000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_range = (400+273.15, 600+273.15) # Temperature range to explore, K
    melt_model = ""

    solution_phases = "Opx(HP)\nO(HP)\nF\n"
    excludes = ""

    # Configure (run build and vertex)
    @info "perplex_configure_pseudosection:"
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
        # print("T (K):\t")
        # println(seismic["T(K)"])
        @test isapprox(seismic["T(K)"], [673.15, 686.483, 699.817, 713.15, 726.483, 739.817, 753.15, 766.483, 779.817, 793.15, 806.483, 819.817, 833.15, 846.483, 859.817, 873.15], atol=0.1)
    end
    @test haskey(seismic, "rho,kg/m3")
    if haskey(seismic, "rho,kg/m3")
        @test !any(x->x<2700, seismic["rho,kg/m3"]) && !any(x->x>3200, seismic["rho,kg/m3"])
        # print("rho,kg/m3:\t")
        # println(seismic["rho,kg/m3"])
        @test isapprox(seismic["rho,kg/m3"], [2875.21, 2875.21, 2875.21, 2875.21, 2875.2, 2875.34, 2875.32, 2905.4, 2927.22, 2933.56, 2933.46, 2933.36, 2933.26, 2933.15, 2956.4, 2962.18], atol=0.1)
    end

end

