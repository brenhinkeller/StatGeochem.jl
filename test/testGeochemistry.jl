## --- General conversions

    # Europium anomalies
    @test eustar(6.5433, 5.9037) ≈ 2.0329978601003864
    @test eustar(34.7773, 6.5433, 5.9037, 0.8904) ≈ 2.0825737578695205

    # Iron oxide conversions
    @test feoconversion(3.5, NaN, NaN, NaN) == 3.5
    @test feoconversion(3.5, NaN, 7.5, NaN) == 7.5
    @test feoconversion(3.5, NaN, 7.5, 10) == 7.5
    @test feoconversion(3.5, 4.4, NaN, NaN) ≈ 7.45916511675966
    @test feoconversion(NaN, 4.4, NaN, NaN) ≈ 3.9591651167596607

    # Other oxide conversion
    D = elementify(["Fe" "Mg" "Ca" "P"; 10000 10000 10000 10000; 10000 10000 10000 10000], importas=:Dict)
    D = oxideconversion(D)
    @test all(D["FeOT"] .≈ (molarmass["Fe"]+molarmass["O"])/molarmass["Fe"])
    @test all(D["MgO"] .≈ (molarmass["Mg"]+molarmass["O"])/molarmass["Mg"])
    @test all(D["CaO"] .≈ (molarmass["Ca"]+molarmass["O"])/molarmass["Ca"])
    @test all(D["P2O5"] .≈ (molarmass["P"]+2.5*molarmass["O"])/molarmass["P"])

    # Weathering indices
    @test CIA(14.8577, 4.5611, 3.29641, 2.3992) ≈ 47.66582778067264
    @test WIP(3.2964, 4.5611, 2.3992, 5.9121) ≈ 78.40320264846837


## -- Perplex name abbreviations

    abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm", "oilm", "lime", "mft", "mt", "mang", "bunsn", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi", "Augite(G)", "Cpx(JH)", "Cpx(l)", "Cpx(h)", "Cpx(stx)", "Cpx(stx7)", "Omph(HP)", "Cpx(HP)", "Cpx(m)", "Cpx(stx8)", "Omph(GHP)", "cAmph(G)", "Cumm", "Gl", "Tr", "GlTrTsPg", "Amph(DHP)", "Amph(DPW)", "Ca-Amph(D)", "Na-Amph(D)", "Act(M)", "GlTrTsMr", "cAmph(DP)", "melt(G)", "melt(W)", "melt(HP)", "melt(HGP)", "pMELTS(G)", "mMELTS(G)", "LIQ(NK)", "LIQ(EF)", "Chl(W)", "Chl(HP)", "Chl(LWV)", "O(JH)", "O(SG)", "O(HP)", "O(HPK)", "O(stx)", "O(stx7)", "Ol(m)", "O(stx8)", "Sp(JH)", "GaHcSp", "Sp(JR)", "Sp(GS)", "Sp(HP)", "Sp(stx)", "CrSp", "Sp(stx7)", "Sp(WPC)", "Sp(stx8)", "Pl(JH)", "Pl(h)", "Pl(stx8)", "Kf", "San", "San(TH)", "Grt(JH)", "Gt(W)", "CrGt", "Gt(MPF)", "Gt(B)", "Gt(GCT)", "Gt(HP)", "Gt(EWHP)", "Gt(WPH)", "Gt(stx)", "Gt(stx8)", "Gt(WPPH)", "ZrGt(KP)", "Maj", "Opx(JH)", "Opx(W)", "Opx(HP)", "CrOpx(HP)", "Opx(stx)", "Opx(stx8)", "Mica(W)", "Pheng(HP)", "MaPa", "Mica(CF)", "Mica(CHA1)", "Mica(CHA)", "Mica+(CHA)", "Mica(M)", "Mica(SGH)", "Ctd(W)", "Ctd(HP)", "Ctd(SGH)", "St(W)", "St(HP)", "Bi(W)", "Bio(TCC)", "Bio(WPH)", "Bio(HP)", "Crd(W)", "hCrd", "Sa(WP)", "Sapp(HP)", "Sapp(KWP)", "Sapp(TP)", "Osm(HP)", "F", "F(salt)", "COH-Fluid", "Aq_solven0", "WADDAH", "T", "Scap", "Carp", "Carp(M)", "Carp(SGH)", "Sud(Livi)", "Sud", "Sud(M)", "Anth", "o-Amph", "oAmph(DP)", "feldspar", "feldspar_B", "Pl(I1,HP)", "Fsp(C1)", "Do(HP)", "M(HP)", "Do(AE)", "Cc(AE)", "oCcM(HP)", "Carb(M)", "oCcM(EF)", "dis(EF)", "IlHm(A)", "IlGkPy", "Ilm(WPH)", "Ilm(WPH0)", "Neph(FB)", "Chum", "Atg(PN)", "B", "Pu(M)", "Stlp(M)", "Wus",)
    common_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermakite", "diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "Mg-tschermakite", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "Fe-carpholite", "Fe-sapphirine(793)", "Mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite", "clinochlore", "daphnite", "Fe-sudoite", "Mn-chlorite", "sudoite", "antigorite", "chrysotile", "Fe-talc", "kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)", "lime", "magnesioferrite", "magnetite", "manganosite", "nickel oxide", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "chlorite", "chlorite", "chlorite", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "plagioclase", "plagioclase", "plagioclase", "k-feldspar", "k-feldspar", "k-feldspar", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "chloritoid", "chloritoid", "chloritoid", "staurolite", "staurolite", "biotite", "biotite", "biotite", "biotite", "cordierite", "cordierite", "sapphirine", "sapphirine", "sapphirine", "sapphirine", "osumilite", "fluid", "fluid", "fluid", "fluid", "fluid", "talc", "scapolite", "carpholite", "carpholite", "carpholite", "sudoite", "sudoite", "sudoite", "orthoamphibole", "orthoamphibole", "orthoamphibole", "ternary feldspar", "ternary feldspar", "ternary feldspar", "ternary feldspar", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "ilmenite", "ilmenite", "ilmenite", "ilmenite", "nepheline", "clinohumite", "serpentine", "brucite", "pumpellyite", "stilpnomelane", "wüstite",)

    @test perplex_common_name.(abbreviations) == common_names

    abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm","oilm","lime", "mft", "mt", "mang", "bunsn", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi",)
    full_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote(ordered)", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite-bredigite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermaks pyroxene", "Diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "mg-tschermak", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite(Na-free)", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "fe-carpholite", "fe-sapphirine(793)", "mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite(14Ang)", "clinochlore(ordered)", "daphnite", "Fe-sudoite", "Mn-chlorite", "Sudoite", "antigorite", "chrysotile", "Fe-talc", "Kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)","lime", "magnesioferrite", "magnetite", "manganosite", "nickel oxide", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "Forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica",)

    @test perplex_expand_name.(abbreviations) == full_names
    @test perplex_abbreviate_name.(full_names) == abbreviations
    @test perplex_phase_is_solid.(("melt(HGP)", "q", "diL", "andr", "T(K)")) == (false, true, false, true, false)

    @test findall(germ_perplex_name_matches.(germ_kd["minerals"], germ_kd["minerals"])) == [3, 12, 18]

## --- Saturation models

    #          SiO2,  TiO2,  Al2O3,  FeOT,   MnO,   MgO,   CaO,   Na2O,  K2O,  P2O5
    majors = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 0.279]
    @test tzirc(majors..., 100) ≈ 602.8489762809595
    @test tzircZr(majors..., 800) ≈ 832.9689080567883
    @test all(tzircM((repeat([m],2) for m in majors)...,) .≈ 2.335918319204001)

    @test StatGeochem.Ayers_tsphene(majors...) ≈ 637.139776663209
    @test StatGeochem.Ayers_tspheneTiO2(majors..., 800) ≈ 2.3545537746637324
    @test all(StatGeochem.Ayers_tspheneC.((repeat([m],2) for m in majors)...,) .≈ 2.4338232746497326)

    #                SiO2,   TiO2,  Al2O3,  FeOT,   MgO,  CaO,   Na2O,  K2O,   Li2O, H2O
    montel_elems = [58.509, 1.022, 14.858, 4.371, 4.561, 5.912, 3.296, 2.399, 0.01, 4.0]
    @test StatGeochem.Montel_tmonaziteREE(montel_elems..., 750.0) ≈ 12.03834338792398
    @test StatGeochem.Montel_tmonazite(montel_elems..., 100,100,100,0,0,0) ≈ 630.4499586271999

    @test StatGeochem.Rusiecka_tmonaziteREE(200, 750) ≈ 0.27430570654873154
    @test StatGeochem.Rusiecka_txenotimeY(200, 750) ≈ 41.9312030248943

    @test StatGeochem.Harrison_tapatiteP2O5(58.509, 14.858, 5.912, 3.296, 2.399, 750.) ≈ 0.10142278764336987
    @test StatGeochem.Harrison_tapatiteP(58.509, 14.858, 5.912, 3.296, 2.399, 750.) ≈ 442.6362451135793
    @test StatGeochem.Harrison_tapatiteP2O5(58.509, 750.) ≈ 0.10142278764336987
    @test StatGeochem.Harrison_tapatite(58.509, 0.1) ≈ 748.6127179814277

    #          SiO2,  TiO2,  Al2O3,  FeOT,   MgO,   CaO,   Na2O,  K2O,  P2O5
    majors = [58.509, 1.022, 14.858, 4.371, 4.561, 5.912, 3.296, 2.399, 0.279]
    @test StatGeochem.Tollari_tapatite(majors...) ≈ 521.0594433599132
    @test StatGeochem.Tollari_tapatiteP2O5(58.509,5.912,750.) ≈ 0.5011681927262436

## -- Test melts

if Sys.islinux()
    # Which version of Melts to use
    alphameltsversion = "linux_alphamelts_1-9"

    # Construct file path
    meltsdir = joinpath(resourcepath, alphameltsversion)
    filepath = joinpath(resourcepath, alphameltsversion*".zip")

    # Download precompiled executable
    if ~isfile(filepath)
        @info "Downloading alphamelts to $meltsdir"
        run(`mkdir -p $meltsdir`)
        Downloads.download("https://storage.googleapis.com/statgeochem/$alphameltsversion.zip", filepath)
        run(`unzip -o $filepath -d $resourcepath`)
        run(`mv $meltsdir/alphamelts_linux64 $meltsdir/alphamelts`)
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
    melt_comp = melts_query_liquid(scratchdir, index=1)
    solid_comp = melts_query_solid(scratchdir, index=1)
    modes = melts_query_modes(scratchdir, index=1)

    @test isa(melt_comp, NamedTuple)
    @test isa(solid_comp, NamedTuple)
    @test isa(modes, NamedTuple)

    @test all(melt_comp.Pressure .== 20000.0)
    println(melt_comp.Temperature)
    # @test melt_comp.Temperature ≈ 1624.07:-10:804.07
    println(melt_comp.SiO2)
    # @test isapprox(melt_comp.SiO2, [44.7991, 45.1222, 45.4618, 45.9044, 46.3448, 46.7903, 45.6759, 45.2694, 44.9009, 44.5677, 44.2675, 43.9981, 43.7575, 43.5433, 43.353, 43.1841, 42.9985, 42.8732, 42.8038, 42.7829, 42.8033, 42.8579, 42.9404, 43.0452, 43.1668, 43.3079, 43.5205, 43.7351, 43.9481, 44.1556, 44.3545, 44.5415, 44.7143, 44.8705, 45.0085, 45.1269, 45.2245, 45.3007, 45.3549, 45.3868, 45.3963, 45.3833, 45.3478, 45.1948, 44.9076, 44.5701, 44.3443, 44.1384, 43.9297, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], nans=true)

    @test all(solid_comp.Pressure .== 20000.0)
    println(solid_comp.Temperature)
    # @test solid_comp.Temperature ≈ 1624.07:-10:804.07
    println(solid_comp.SiO2)
    # @test isapprox(solid_comp.SiO2, [NaN, 41.5019, 41.4897, 41.2257, 41.0972, 41.0056, 43.718, 44.3492, 44.7297, 44.9694, 45.1228, 45.2202, 45.2799, 45.3137, 45.3295, 45.3326, 45.2928, 45.2409, 45.1861, 45.1333, 45.0848, 45.0416, 45.0037, 44.9709, 44.9427, 44.9182, 44.8933, 44.8727, 44.8558, 44.842, 44.8309, 44.8219, 44.8148, 44.8093, 44.8052, 44.8023, 44.8004, 44.7993, 44.799, 44.7993, 44.8001, 44.8014, 44.803, 44.8069, 44.813, 44.817, 44.8171, 44.8169, 44.8167, 44.8167, 44.8169, 44.8171, 44.8173, 44.8175, 44.8177, 44.8179, 44.8181, 44.8182, 44.8184, 44.8185, 44.8187, 44.8188, 44.8189, 44.819, 44.8191, 44.8192, 44.8193, 44.8193, 44.8194, 44.8195, 44.8195, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8195, 44.8196, 44.8197, 44.8198, 44.8199], nans=true)

    @test all(modes.Pressure .== 20000.0)
    println(modes.Temperature)
    # @test modes.Temperature ≈ 1624.07:-10:804.07
    println(modes.liquid_0)
    # @test isapprox(modes.liquid_0, [99.95749, 91.072189, 83.344261, 76.418306, 70.599891, 65.643435, 55.505785, 49.595839, 44.677608, 40.507414, 36.917284, 33.787818, 31.031719, 28.583347, 26.391889, 24.41682, 21.131799, 18.28986, 15.881062, 13.852888, 12.149205, 10.716689, 9.508346, 8.484562, 7.61275, 6.862884, 6.194069, 5.622428, 5.130364, 4.70406, 4.332516, 4.006874, 3.719932, 3.465787, 3.239572, 3.037248, 2.855453, 2.691375, 2.542652, 2.407295, 2.283625, 2.170214, 2.065851, 1.950841, 1.825706, 0.68822, 0.268833, 0.106921, 0.027371, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true)
    println(modes.olivine_0)
    # @test isapprox(modes.olivine_0, [0.0, 8.882234, 16.607219, 23.396165, 29.111411, 33.979596, 36.01218, 38.364714, 40.307067, 41.93387, 43.311147, 44.486449, 45.49502, 46.363698, 47.11345, 47.761057, 48.588784, 49.290361, 49.880834, 50.37841, 50.798768, 51.155082, 51.458212, 51.717009, 51.938661, 52.150643, 52.510356, 52.826652, 53.106072, 53.353971, 53.574796, 53.772272, 53.949552, 54.109313, 54.253844, 54.385101, 54.504761, 54.614264, 54.714847, 54.807573, 54.893359, 54.972993, 55.047157, 55.101579, 55.138126, 54.216269, 53.84507, 53.698077, 53.628726, 53.623437, 53.658848, 53.692537, 53.72459, 53.755083, 53.784087, 53.811662, 53.837863, 53.862741, 53.886339, 53.908697, 53.930719, 53.992584, 54.056671, 54.123351, 54.193047, 54.266249, 54.343534, 54.425581, 54.513207, 54.607395, 54.709348, 54.820555, 54.942876, 55.078676, 55.230995, 55.403811, 55.602415, 55.833981, 56.108425, 56.206913, 56.210501, 56.213232, 56.215135], nans=true)

end

## --- End of File
