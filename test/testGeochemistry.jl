## --- General conversions

    # Europium anomalies
    @test eustar(6.5433, 5.9037) ≈ 2.0329978601003864
    @test eustar(34.7773, 6.5433, 5.9037, 0.8904) ≈ 2.0825737578695205

    # Iron oxide conversions
    @test feoconversion(3.5, NaN, NaN, NaN) == 3.5
    @test feoconversion(3.5, NaN, 7.5, NaN) == 7.5
    @test feoconversion(3.5, NaN, 7.5, 10) == 7.5
    @test feoconversion(3.5, NaN, NaN, 10) == 8.998102538090137
    @test feoconversion(3.5, 4.4, NaN, NaN) ≈ 7.45916511675966
    @test feoconversion(NaN, 4.4, NaN, NaN) ≈ 3.9591651167596607
    @test isnan(feoconversion(NaN, NaN, NaN, NaN))

    # Other oxide conversion
    D = ["Fe" "Mg" "Ca" "P" "FeOT" "MgO" "CaO" "P2O5"; 10000 10000 10000 10000 NaN NaN NaN NaN; 10000 10000 10000 10000 NaN NaN NaN NaN]
    M = elementify(D, importas=:Dict)
    O = oxideconversion(M)
    @test all(O["FeOT"] .≈ (molarmass["Fe"]+molarmass["O"])/molarmass["Fe"])
    @test all(O["MgO"] .≈ (molarmass["Mg"]+molarmass["O"])/molarmass["Mg"])
    @test all(O["CaO"] .≈ (molarmass["Ca"]+molarmass["O"])/molarmass["Ca"])
    @test all(O["P2O5"] .≈ (molarmass["P"]+2.5*molarmass["O"])/molarmass["P"])

    DT = deepcopy(TupleDataset(M))
    oxideconversion!(DT)
    @test all(DT.FeOT .≈ (molarmass["Fe"]+molarmass["O"])/molarmass["Fe"])
    @test all(DT.MgO .≈ (molarmass["Mg"]+molarmass["O"])/molarmass["Mg"])
    @test all(DT.CaO .≈ (molarmass["Ca"]+molarmass["O"])/molarmass["Ca"])
    @test all(DT.P2O5 .≈ (molarmass["P"]+2.5*molarmass["O"])/molarmass["P"])

    D = ["NiO" "CoO" "BaO" "SO3" "Ni" "Co" "Ba" "S"; 1 1 1 1 NaN NaN NaN NaN; 1 1 1 1 NaN NaN NaN NaN]
    O = elementify(D, importas=:Dict)
    M = metalconversion(O)
    @test all(M["Ni"] .≈ 10000molarmass["Ni"]/(molarmass["Ni"]+molarmass["O"]))
    @test all(M["Co"] .≈ 10000molarmass["Co"]/(molarmass["Co"]+molarmass["O"]))
    @test all(M["Ba"] .≈ 10000molarmass["Ba"]/(molarmass["Ba"]+molarmass["O"]))
    @test all(M["S"] .≈ 10000molarmass["S"]/(molarmass["S"]+3molarmass["O"]))


    # Carbonate conversions
    D = ["CaCO3" "MgCO3" "CaO" "MgO" "CO2" "TOC" "TIC" "TC" "C"; 1 1 NaN NaN NaN NaN NaN NaN 10000; 1 1 NaN NaN NaN NaN NaN NaN 10000]
    C = elementify(D, importas=:Dict)
    carbonateconversion!(C)
    @test all(C["MgO"] .≈ (molarmass["Mg"]+molarmass["O"])/(molarmass["Mg"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["CaO"] .≈ (molarmass["Ca"]+molarmass["O"])/(molarmass["Ca"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["CO2"] .≈ (molarmass["C"]+2molarmass["O"])/(molarmass["Ca"]+molarmass["C"]+3molarmass["O"]) + (molarmass["C"]+2molarmass["O"])/(molarmass["Mg"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["TC"] .≈ 1)
    @test all(C["TIC"] .≈ 0.9616817911685506*molarmass["C"]/(molarmass["C"] + 2molarmass["O"]))
    @test all(C["TOC"] .≈ 1 - 0.9616817911685506*molarmass["C"]/(molarmass["C"] + 2molarmass["O"]))

    D = ["CaCO3" "MgCO3" "CaO" "MgO" "CO2" "TOC" "TIC" "TC" "C"; 1 1 NaN NaN NaN NaN NaN NaN NaN; 1 1 NaN NaN NaN NaN NaN NaN NaN]
    C = elementify(D, importas=:Dict)
    carbonateconversion!(C)
    @test all(C["MgO"] .≈ (molarmass["Mg"]+molarmass["O"])/(molarmass["Mg"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["CaO"] .≈ (molarmass["Ca"]+molarmass["O"])/(molarmass["Ca"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["CO2"] .≈ (molarmass["C"]+2molarmass["O"])/(molarmass["Ca"]+molarmass["C"]+3molarmass["O"]) + (molarmass["C"]+2molarmass["O"])/(molarmass["Mg"]+molarmass["C"]+3molarmass["O"]))
    @test all(C["TIC"] .≈ 0.9616817911685506*molarmass["C"]/(molarmass["C"] + 2molarmass["O"]))
    @test all(C["C"] .≈ 1e4*0.9616817911685506*molarmass["C"]/(molarmass["C"] + 2molarmass["O"]))

    # Weathering indices
    @test CIA(14.8577, 4.5611, 3.29641, 2.3992) ≈ 47.66582778067264
    @test WIP(3.2964, 4.5611, 2.3992, 5.9121) ≈ 78.40320264846837


## -- Norms

    n = StatGeochem.cipw_norm(57.05, 0.44, 14.57, 8.02, 0, 0.18, 6.79, 10.55, 1.26, 0.49, 0.06)
    @test n.quartz ≈ 22.64883050985674
    @test n.orthoclase ≈ 2.895716060129942
    @test n.plagioclase ≈ 43.31401932989578
    @test n.corundum ≈ 0.0
    @test n.nepheline ≈ 0.0
    @test n.diopside ≈ 13.731414049010224
    @test n.orthopyroxene ≈ 3.94200752101842
    @test n.olivine ≈ 0.0
    @test n.magnetite ≈ 11.62824215584731
    @test n.ilmenite ≈ 0.8356557044661498
    @test n.apatite ≈ 0.1421161651208747


## --- Saturation models

    #          SiO2,  TiO2,  Al2O3,  FeOT,   MnO,   MgO,   CaO,   Na2O,  K2O,  P2O5
    majors = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 0.279]
    @test tzirc(majors..., 100) ≈ 603.4774053095614
    @test tzircZr(majors..., 800) ≈ 826.1071302971219
    @test mean(tzircM((repeat([m],2) for m in majors)...,)) ≈ 2.328787411099651

    @test StatGeochem.Ayers_tsphene_18(majors...) ≈ 637.7486728299519
    @test StatGeochem.Ayers_tspheneTiO2_18(majors..., 800) ≈ 2.3486842447760026
    @test mean(StatGeochem.Ayers_tspheneC.((repeat([m],2) for m in majors)...,)) ≈ 2.4263934899817188

    @test StatGeochem.Ayers_tsphene(majors...) ≈ 694.4428983217797
    @test StatGeochem.Ayers_tspheneTiO2(majors..., 800) ≈ 1.5286740880554586

    #                SiO2,   TiO2,  Al2O3,  FeOT,   MgO,  CaO,   Na2O,  K2O,   Li2O, H2O
    montel_elems = [58.509, 1.022, 14.858, 4.371, 4.561, 5.912, 3.296, 2.399, 0.01, 4.0]
    @test StatGeochem.Montel_tmonaziteREE(montel_elems..., 750.0) ≈ 11.884450325741755
    @test StatGeochem.Montel_tmonazite(montel_elems..., 100,100,100,0,0,0) ≈ 631.2376817530326

    @test StatGeochem.Rusiecka_tmonaziteREE(200, 750) ≈ 0.27430570654873154
    @test StatGeochem.Rusiecka_txenotimeY(200, 750) ≈ 41.9312030248943

    @test StatGeochem.Harrison_tapatiteP2O5(58.509, 14.858, 5.912, 3.296, 2.399, 750.) ≈ 0.10142278764336987
    @test StatGeochem.Harrison_tapatiteP(58.509, 14.858, 5.912, 3.296, 2.399, 750.) ≈ 442.6362451135793
    @test StatGeochem.Harrison_tapatiteP2O5(58.509, 750.) ≈ 0.10142278764336987
    @test StatGeochem.Harrison_tapatite(58.509, 0.1) ≈ 748.6127179814277

    #          SiO2,  TiO2,  Al2O3,  FeOT,   MgO,   CaO,   Na2O,  K2O,  P2O5
    majors = [58.509, 1.022, 14.858, 4.371, 4.561, 5.912, 3.296, 2.399, 0.279]
    @test StatGeochem.Tollari_tapatite(majors...) ≈ 528.5868109033141
    @test StatGeochem.Tollari_tapatiteP2O5(58.509,5.912,750.) ≈ 0.5011681927262436

    @test StatGeochem.Hayden_trutile(majors...) ≈ 822.7645622408794
    @test StatGeochem.Hayden_trutileTiO2(majors...,750.) ≈ 0.045228791859556305

## -- Test thermometers

    @test StatGeochem.Ferry_Zr_in_rutile(750,1) ≈ 982.8714076786658
    @test StatGeochem.Ferry_Zr_in_rutileT(982.8714076786658,1) ≈ 750.0
    @test StatGeochem.Ferry_Ti_in_zircon(750,1,1) ≈ 10.46178465494583
    @test StatGeochem.Ferry_Ti_in_zirconT(10.46178465494583, 1, 1) ≈ 750.0

    @test StatGeochem.Crisp_Ti_in_zircon(750,0,1,1) ≈ 14.08608953046849
    @test StatGeochem.Crisp_Ti_in_zircon(750,1000,1,1) ≈ 13.703165806686624
    @test StatGeochem.Crisp_Ti_in_zircon(750,10000,1,1) ≈ 10.570721163513461
    @test StatGeochem.Crisp_Ti_in_zircon(750,20000,1,1) ≈ 7.4184040803703954

## -- Test melts

if Sys.islinux() || Sys.isapple()
    # Which version of Melts to use
    if Sys.isapple()
        alphameltsversion = "macosx_alphamelts_1-9"
        alphameltsexec = "alphamelts_macosx64"
    else
        alphameltsversion = "linux_alphamelts_1-9"
        alphameltsexec = "alphamelts_linux64"
    end

    # Construct file path
    meltsdir = joinpath(resourcepath, alphameltsversion)
    filepath = joinpath(resourcepath, alphameltsversion*".zip")

    # Download precompiled executable
    if !isfile(filepath)
        @info "Downloading alphamelts to $meltsdir"
        run(`mkdir -p $meltsdir`)
        Downloads.download("https://storage.googleapis.com/statgeochem/$alphameltsversion.zip", filepath)
        run(`unzip -o $filepath -d $resourcepath`)
        run(`mv $meltsdir/$alphameltsexec $meltsdir/alphamelts`)
    end

    meltspath = joinpath(meltsdir, "run_alphamelts.command")
    scratchdir = "./"

    ## --- # # # # # # # # # # # pMelts equil. batch melting # # # # # # # # # # # #

    # Conditions
    P_range = [20000,20000]
    T_range = [1700,800]
    # Starting composition
    elements = ["SiO2",  "TiO2","Al2O3","Fe2O3","Cr2O3",  "FeO",  "MnO",  "MgO",   "NiO",  "CoO",  "CaO",  "Na2O", "K2O", "P2O5", "H2O",]
    composition=[44.8030, 0.1991, 4.4305, 0.9778, 0.3823, 7.1350, 0.1344, 37.6345, 0.2489, 0.0129, 3.5345, 0.3584, 0.0289, 0.0209, 0.15,] #mcdbse (McDonough Pyrolite)
    # Run simulation
    melts_configure(meltspath, scratchdir, composition, elements, T_range, P_range,
        batchstring="1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n",
        dT=-10, dP=0, index=1, version="pMELTS",mode="isobaric",fo2path="FMQ")

    # Read results
    melt_comp = melts_query_liquid(scratchdir, index=1, importas=:Tuple)
    solid_comp = melts_query_solid(scratchdir, index=1, importas=:Tuple)
    sys = melts_query_system(scratchdir, index=1, importas=:Tuple)
    modes = melts_query_modes(scratchdir, index=1, importas=:Tuple)
    clean_modes = melts_clean_modes(scratchdir, index=1)
    bulk = melts_query(scratchdir, index=1)

    @test isa(melt_comp, NamedTuple)
    @test isa(solid_comp, NamedTuple)
    @test isa(sys, NamedTuple)
    @test isa(modes, NamedTuple)
    @test isa(clean_modes, Dict)
    @test isa(bulk, Dict)

    @test all(melt_comp.Pressure .== 20000.0)
    # print("melt_comp.Temperature: ")
    # println(melt_comp.Temperature)
    @test melt_comp.Temperature ≈ 1624.06:-10.0:914.06 
    # print("melt_comp.SiO2: ")
    # println(melt_comp.SiO2)
    @test isapprox(melt_comp.SiO2, [44.7993, 45.1225, 45.4621, 45.9048, 46.3452, 46.7906, 45.6756, 45.2691, 44.9006, 44.5674, 44.2672, 43.9979, 43.7574, 43.5432, 43.3529, 43.184, 42.9984, 42.8732, 42.8037, 42.7829, 42.8033, 42.8579, 42.9405, 43.0453, 43.1669, 43.308, 43.5206, 43.7353, 43.9483, 44.1558, 44.3546, 44.5417, 44.7144, 44.8707, 45.0086, 45.1269, 45.2246, 45.3007, 45.3549, 45.3868, 45.3922, 45.1936, 44.9836, 44.76, 44.3911, 44.2064, 44.0531, 43.8942, 43.7199, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], nans=true)

    @test all(solid_comp.Pressure .== 20000.0)
    # print("solid_comp.Temperature: ")
    # println(solid_comp.Temperature)
    @test solid_comp.Temperature ≈ 1624.06:-10.0:914.06 
    # print("solid_comp.SiO2: ")
    # println(solid_comp.SiO2)
    @test isapprox(solid_comp.SiO2, [41.4817, 41.5019, 41.4897, 41.2256, 41.0971, 41.0055, 43.7186, 44.3496, 44.73, 44.9695, 45.1228, 45.2202, 45.28, 45.3137, 45.3295, 45.3326, 45.2928, 45.2409, 45.1861, 45.1333, 45.0848, 45.0416, 45.0037, 44.9709, 44.9427, 44.9181, 44.8932, 44.8727, 44.8558, 44.842, 44.8308, 44.8219, 44.8148, 44.8093, 44.8052, 44.8023, 44.8004, 44.7993, 44.799, 44.7993, 44.8002, 44.8057, 44.8108, 44.8155, 44.8207, 44.8183, 44.8174, 44.8169, 44.8166, 44.8167, 44.8169, 44.8171, 44.8173, 44.8175, 44.8177, 44.8179, 44.8181, 44.8182, 44.8184, 44.8186, 44.8187, 44.8188, 44.8189, 44.819, 44.8191, 44.8192, 44.8193, 44.8193, 44.8194, 44.8195, 44.8195, 44.8196], nans=true)

    @test all(sys.Pressure .== 20000.0)
    @test sys.Temperature ≈ 1624.06:-10.0:914.06 
    # print("system.aH2O: ")
    # println(sys.aH2O)
    @test isapprox(sys.aH2O, [0.00119355, 0.00143907, 0.00171302, 0.0020253, 0.00235407, 0.00269726, 0.00378774, 0.00472133, 0.00577799, 0.00696703, 0.00829895, 0.00978529, 0.0114387, 0.0132728, 0.0153034, 0.0175484, 0.0224613, 0.0285535, 0.0358671, 0.044436, 0.0542559, 0.0652938, 0.0774963, 0.0907963, 0.105118, 0.120441, 0.137091, 0.15441, 0.172315, 0.190729, 0.209584, 0.228825, 0.248403, 0.26828, 0.288424, 0.308812, 0.329427, 0.350254, 0.371286, 0.392517, 0.413978, 0.437322, 0.461004, 0.484995, 0.517313, 0.578275, 0.631048, 0.677779, 0.720523, 0.74231, 0.75994, 0.778542, 0.798166, 0.818869, 0.840709, 0.863755, 0.888076, 0.913751, 0.940863, 0.969506, 0.999777, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], nans=true)

    @test all(modes.Pressure .== 20000.0)
    # print("modes.Temperature: ")
    # println(modes.Temperature)
    @test modes.Temperature ≈ 1624.06:-10.0:914.06 
    # print("modes.liquid_0: ")
    # println(modes.liquid_0)
    @test isapprox(modes.liquid_0, [99.953083, 91.06567, 83.3387, 76.413373, 70.595706, 65.639854, 55.500696, 49.591648, 44.674084, 40.504401, 36.914672, 33.785529, 31.029693, 28.581541, 26.390266, 24.415352, 21.129398, 18.287815, 15.879336, 13.851438, 12.147987, 10.715663, 9.507479, 8.483825, 7.612121, 6.862318, 6.193587, 5.622014, 5.130007, 4.70375, 4.332245, 4.006636, 3.719721, 3.4656, 3.239405, 3.037098, 2.855318, 2.691253, 2.542541, 2.407194, 2.282555, 2.128077, 1.991449, 1.869639, 1.357812, 0.473618, 0.196578, 0.07415, 0.009783, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true)
    # print("modes.olivine_0: ")
    # println(modes.olivine_0)
    @test isapprox(modes.olivine_0, [0.004405, 8.88875, 16.612778, 23.40101, 29.115521, 33.983114, 36.014211, 38.366376, 40.308451, 41.935036, 43.312139, 44.487298, 45.495751, 46.364327, 47.113994, 47.761526, 48.589381, 49.290863, 49.881257, 50.378767, 50.799069, 51.155338, 51.45843, 51.717196, 51.938821, 52.150943, 52.510619, 52.826884, 53.106277, 53.354154, 53.574959, 53.772418, 53.949683, 54.109432, 54.253952, 54.385199, 54.50485, 54.614346, 54.714922, 54.807643, 54.892635, 54.939162, 54.987413, 55.036388, 54.729708, 53.998944, 53.756472, 53.649289, 53.596897, 53.618702, 53.654148, 53.687873, 53.719962, 53.750492, 53.779533, 53.807145, 53.833385, 53.858301, 53.881939, 53.904336, 53.925529, 53.98708, 54.051137, 54.11778, 54.187431, 54.26058, 54.337802, 54.419774, 54.507309, 54.601389, 54.703211, 54.814258], nans=true)

    # Test `clean_modes` vs `modes`
    @test clean_modes["Temperature"] == modes.Temperature
    @test clean_modes["liquid"] == modes.liquid_0
    @test clean_modes["olivine"] == modes.olivine_0
    @test clean_modes["apatite"] == modes.apatite
end


## --- End of File
