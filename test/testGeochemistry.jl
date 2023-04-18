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
    @test tzirc(majors..., 100) ≈ 603.4774053095614
    @test tzircZr(majors..., 800) ≈ 826.1071302971219
    @test mean(tzircM((repeat([m],2) for m in majors)...,)) ≈ 2.328787411099651

    @test StatGeochem.Ayers_tsphene(majors...) ≈ 637.7486728299519
    @test StatGeochem.Ayers_tspheneTiO2(majors..., 800) ≈ 2.3486842447760026
    @test mean(StatGeochem.Ayers_tspheneC.((repeat([m],2) for m in majors)...,)) ≈ 2.4263934899817188

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

    @test StatGeochem.Ferry_Ti_in_zircon(750,1,1) ≈ 10.46178465494583


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
    system = melts_query_system(scratchdir, index=1, importas=:Tuple)
    modes = melts_query_modes(scratchdir, index=1, importas=:Tuple)
    clean_modes = melts_clean_modes(scratchdir, index=1)
    bulk = melts_query(scratchdir, index=1)

    @test isa(melt_comp, NamedTuple)
    @test isa(solid_comp, NamedTuple)
    @test isa(system, NamedTuple)
    @test isa(modes, NamedTuple)
    @test isa(clean_modes, Dict)
    @test isa(bulk, Dict)

    @test all(melt_comp.Pressure .== 20000.0)
    print("melt_comp.Temperature: ")
    println(melt_comp.Temperature)
    @test melt_comp.Temperature ≈ 1624.06:-10.0:914.06 # 624.07:-10:804.07
    print("melt_comp.SiO2: ")
    println(melt_comp.SiO2)
    @test isapprox(melt_comp.SiO2, [44.7993, 45.1225, 45.4621, 45.9048, 46.3452, 46.7906, 45.6756, 45.2691, 44.9006, 44.5674, 44.2672, 43.9979, 43.7574, 43.5432, 43.3529, 43.184, 42.9984, 42.8732, 42.8037, 42.7829, 42.8033, 42.8579, 42.9405, 43.0453, 43.1669, 43.308, 43.5206, 43.7353, 43.9483, 44.1558, 44.3546, 44.5417, 44.7144, 44.8707, 45.0086, 45.1269, 45.2246, 45.3007, 45.3549, 45.3868, 45.3922, 45.1936, 44.9836, 44.76, 44.3911, 44.2064, 44.0531, 43.8942, 43.7199, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], nans=true)
    # @test isapprox(melt_comp.SiO2, [44.7991, 45.1222, 45.4618, 45.9044, 46.3448, 46.7903, 45.6759, 45.2694, 44.9009, 44.5677, 44.2675, 43.9981, 43.7575, 43.5433, 43.353, 43.1841, 42.9985, 42.8732, 42.8038, 42.7829, 42.8033, 42.8579, 42.9404, 43.0452, 43.1668, 43.3079, 43.5205, 43.7351, 43.9481, 44.1556, 44.3545, 44.5415, 44.7143, 44.8705, 45.0085, 45.1269, 45.2245, 45.3007, 45.3549, 45.3868, 45.3963, 45.3833, 45.3478, 45.1948, 44.9076, 44.5701, 44.3443, 44.1384, 43.9297, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], nans=true)

    @test all(solid_comp.Pressure .== 20000.0)
    print("solid_comp.Temperature: ")
    println(solid_comp.Temperature)
    @test solid_comp.Temperature ≈ 1624.06:-10.0:914.06 # 1624.07:-10:804.07
    print("solid_comp.SiO2: ")
    println(solid_comp.SiO2)
    @test isapprox(solid_comp.SiO2, [41.4817, 41.5019, 41.4897, 41.2256, 41.0971, 41.0055, 43.7186, 44.3496, 44.73, 44.9695, 45.1228, 45.2202, 45.28, 45.3137, 45.3295, 45.3326, 45.2928, 45.2409, 45.1861, 45.1333, 45.0848, 45.0416, 45.0037, 44.9709, 44.9427, 44.9181, 44.8932, 44.8727, 44.8558, 44.842, 44.8308, 44.8219, 44.8148, 44.8093, 44.8052, 44.8023, 44.8004, 44.7993, 44.799, 44.7993, 44.8002, 44.8057, 44.8108, 44.8155, 44.8207, 44.8183, 44.8174, 44.8169, 44.8166, 44.8167, 44.8169, 44.8171, 44.8173, 44.8175, 44.8177, 44.8179, 44.8181, 44.8182, 44.8184, 44.8186, 44.8187, 44.8188, 44.8189, 44.819, 44.8191, 44.8192, 44.8193, 44.8193, 44.8194, 44.8195, 44.8195, 44.8196], nans=true)
    # @test isapprox(solid_comp.SiO2, [NaN, 41.5019, 41.4897, 41.2257, 41.0972, 41.0056, 43.718, 44.3492, 44.7297, 44.9694, 45.1228, 45.2202, 45.2799, 45.3137, 45.3295, 45.3326, 45.2928, 45.2409, 45.1861, 45.1333, 45.0848, 45.0416, 45.0037, 44.9709, 44.9427, 44.9182, 44.8933, 44.8727, 44.8558, 44.842, 44.8309, 44.8219, 44.8148, 44.8093, 44.8052, 44.8023, 44.8004, 44.7993, 44.799, 44.7993, 44.8001, 44.8014, 44.803, 44.8069, 44.813, 44.817, 44.8171, 44.8169, 44.8167, 44.8167, 44.8169, 44.8171, 44.8173, 44.8175, 44.8177, 44.8179, 44.8181, 44.8182, 44.8184, 44.8185, 44.8187, 44.8188, 44.8189, 44.819, 44.8191, 44.8192, 44.8193, 44.8193, 44.8194, 44.8195, 44.8195, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8196, 44.8195, 44.8196, 44.8197, 44.8198, 44.8199], nans=true)

    @test all(system.Pressure .== 20000.0)
    @test system.Temperature ≈ 1624.06:-10.0:914.06 # 1624.07:-10:804.07
    print("system.aH2O: ")
    println(system.aH2O)
    @test isapprox(system.aH2O, [0.00119355, 0.00143907, 0.00171302, 0.0020253, 0.00235407, 0.00269726, 0.00378774, 0.00472133, 0.00577799, 0.00696703, 0.00829895, 0.00978529, 0.0114387, 0.0132728, 0.0153034, 0.0175484, 0.0224613, 0.0285535, 0.0358671, 0.044436, 0.0542559, 0.0652938, 0.0774963, 0.0907963, 0.105118, 0.120441, 0.137091, 0.15441, 0.172315, 0.190729, 0.209584, 0.228825, 0.248403, 0.26828, 0.288424, 0.308812, 0.329427, 0.350254, 0.371286, 0.392517, 0.413978, 0.437322, 0.461004, 0.484995, 0.517313, 0.578275, 0.631048, 0.677779, 0.720523, 0.74231, 0.75994, 0.778542, 0.798166, 0.818869, 0.840709, 0.863755, 0.888076, 0.913751, 0.940863, 0.969506, 0.999777, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], nans=true)

    @test all(modes.Pressure .== 20000.0)
    print("modes.Temperature: ")
    println(modes.Temperature)
    @test modes.Temperature ≈ 1624.06:-10.0:914.06 # 1624.07:-10:804.07
    print("modes.liquid_0: ")
    println(modes.liquid_0)
    @test isapprox(modes.liquid_0, [99.953083, 91.06567, 83.3387, 76.413373, 70.595706, 65.639854, 55.500696, 49.591648, 44.674084, 40.504401, 36.914672, 33.785529, 31.029693, 28.581541, 26.390266, 24.415352, 21.129398, 18.287815, 15.879336, 13.851438, 12.147987, 10.715663, 9.507479, 8.483825, 7.612121, 6.862318, 6.193587, 5.622014, 5.130007, 4.70375, 4.332245, 4.006636, 3.719721, 3.4656, 3.239405, 3.037098, 2.855318, 2.691253, 2.542541, 2.407194, 2.282555, 2.128077, 1.991449, 1.869639, 1.357812, 0.473618, 0.196578, 0.07415, 0.009783, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true)
    # @test isapprox(modes.liquid_0, [99.95749, 91.072189, 83.344261, 76.418306, 70.599891, 65.643435, 55.505785, 49.595839, 44.677608, 40.507414, 36.917284, 33.787818, 31.031719, 28.583347, 26.391889, 24.41682, 21.131799, 18.28986, 15.881062, 13.852888, 12.149205, 10.716689, 9.508346, 8.484562, 7.61275, 6.862884, 6.194069, 5.622428, 5.130364, 4.70406, 4.332516, 4.006874, 3.719932, 3.465787, 3.239572, 3.037248, 2.855453, 2.691375, 2.542652, 2.407295, 2.283625, 2.170214, 2.065851, 1.950841, 1.825706, 0.68822, 0.268833, 0.106921, 0.027371, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nans=true)
    print("modes.olivine_0: ")
    println(modes.olivine_0)
    @test isapprox(modes.olivine_0, [0.004405, 8.88875, 16.612778, 23.40101, 29.115521, 33.983114, 36.014211, 38.366376, 40.308451, 41.935036, 43.312139, 44.487298, 45.495751, 46.364327, 47.113994, 47.761526, 48.589381, 49.290863, 49.881257, 50.378767, 50.799069, 51.155338, 51.45843, 51.717196, 51.938821, 52.150943, 52.510619, 52.826884, 53.106277, 53.354154, 53.574959, 53.772418, 53.949683, 54.109432, 54.253952, 54.385199, 54.50485, 54.614346, 54.714922, 54.807643, 54.892635, 54.939162, 54.987413, 55.036388, 54.729708, 53.998944, 53.756472, 53.649289, 53.596897, 53.618702, 53.654148, 53.687873, 53.719962, 53.750492, 53.779533, 53.807145, 53.833385, 53.858301, 53.881939, 53.904336, 53.925529, 53.98708, 54.051137, 54.11778, 54.187431, 54.26058, 54.337802, 54.419774, 54.507309, 54.601389, 54.703211, 54.814258], nans=true)
    # @test isapprox(modes.olivine_0, [0.0, 8.882234, 16.607219, 23.396165, 29.111411, 33.979596, 36.01218, 38.364714, 40.307067, 41.93387, 43.311147, 44.486449, 45.49502, 46.363698, 47.11345, 47.761057, 48.588784, 49.290361, 49.880834, 50.37841, 50.798768, 51.155082, 51.458212, 51.717009, 51.938661, 52.150643, 52.510356, 52.826652, 53.106072, 53.353971, 53.574796, 53.772272, 53.949552, 54.109313, 54.253844, 54.385101, 54.504761, 54.614264, 54.714847, 54.807573, 54.893359, 54.972993, 55.047157, 55.101579, 55.138126, 54.216269, 53.84507, 53.698077, 53.628726, 53.623437, 53.658848, 53.692537, 53.72459, 53.755083, 53.784087, 53.811662, 53.837863, 53.862741, 53.886339, 53.908697, 53.930719, 53.992584, 54.056671, 54.123351, 54.193047, 54.266249, 54.343534, 54.425581, 54.513207, 54.607395, 54.709348, 54.820555, 54.942876, 55.078676, 55.230995, 55.403811, 55.602415, 55.833981, 56.108425, 56.206913, 56.210501, 56.213232, 56.215135], nans=true)

    # Test `clean_modes` vs `modes`
    @test clean_modes["Temperature"] == modes.Temperature
    @test clean_modes["liquid"] == modes.liquid_0
    @test clean_modes["olivine"] == modes.olivine_0
    @test clean_modes["apatite"] == modes.apatite
end


## --- Test PerpleX

if Sys.isunix()

    # Choose perpleX version
    perplexversion = "perplex-6.8.7"

    # Construct file path
    perplexdir = joinpath(resourcepath, perplexversion)
    scratchdir = "./"

    if Sys.islinux()
        # Download precompiled executable
        if !isfile(joinpath(perplexdir,"vertex"))
            @info "Downloading PerpleX to $perplexdir"
            run(`mkdir -p $perplexdir`)
            file = Downloads.download("https://storage.googleapis.com/statgeochem/$perplexversion-linux.tar.gz",joinpath(resourcepath,"$perplexversion-linux.tar.gz"))
            run(`tar -xzf $file -C $perplexdir`)
        end
    else
        # Compile from source
        if !isfile(joinpath(perplexdir,"vertex"))
            # Check if there is a fortran compiler
            run(`gfortran -v`)

            # Download Perplex v6.8.7 -- known to work with interface used here
            file = Downloads.download("https://storage.googleapis.com/statgeochem/$perplexversion.zip", joinpath(resourcepath,"$perplexversion.zip"))

            # # For a more updated perplex version, you might also try
            # file = download("https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip", joinpath(resourcepath,"perplex-stable.zip"))

            run(`unzip -u $file -d $resourcepath`) # Extract
            system("cd $perplexdir; make") # Compile
        end
    end

    # Kelemen (2014) primitive continental basalt excluding Mn and Ti since most melt models can"t handle them..
    elements =    [ "SIO2", "AL2O3",  "FEO",  "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    composition = [50.0956, 15.3224, 8.5103, 9.2520, 9.6912, 2.5472, 0.8588, 2.0000, 0.6000,]

    # Emphasis on phases from Holland and Powell -- all phases can be used with hp02ver.dat.
    HP_solution_phases = "Omph(HP)\nOpx(HP)\nGlTrTsPg\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nfeldspar_B\nMica(CF)\nBio(TCC)\nChl(HP)\nCtd(HP)\nSapp(HP)\nSt(HP)\nIlHm(A)\nDo(HP)\nT\nB\nF\n"
    HP_excludes = ""

    ## --- # # # # # # # # # # # # # Isobaric example # # # # # # # # # # # # # # # #

    # Input parameters
    P = 1000 # Pressure, bar
    T_range = (0+273.15, 1500+273.15) # Temperature range, Kelvin

    # Configure (run build and vertex)
    melt_model = "melt(HP)"
    @time perplex_configure_isobar(perplexdir, scratchdir, composition, elements, P, T_range,
        dataset="hp02ver.dat",
        npoints=100,
        excludes=HP_excludes,
        solution_phases=melt_model*"\n"*HP_solution_phases
    )

    ## --- Query all properties at a single temperature -- results returned as text

    T = 850+273.15
    data_isobaric = perplex_query_point(perplexdir, scratchdir, T)
    @test isa(data_isobaric, String)

    ## --- Query the full isobar -- results returned as dict

    bulk = perplex_query_system(perplexdir, scratchdir, importas=:Tuple)
    @test isa(bulk, NamedTuple)
    @test haskey(bulk, :SIO2)
    if haskey(bulk, :SIO2)
        print("bulk.SIO2: ")
        println(bulk.SIO2)
        @test haskey(bulk, :SIO2) && all(isapprox.(bulk.SIO2, 50.66433039859823, atol=0.1))
    end

    melt = perplex_query_phase(perplexdir, scratchdir, melt_model, importas=:Tuple)
    @test isa(melt, NamedTuple)
    @test haskey(melt, :SIO2)
    if haskey(melt, :SIO2)
        print("melt.SIO2: ")
        println(melt.SIO2)
        @test !isempty(melt.SIO2) && !any(x->x<45, melt.SIO2) && !any(x->x>75, melt.SIO2)
        # @test isapprox(melt.SIO2, [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 66.78282537747364, 66.82016525351406, 66.82117995364602, 66.80298329925418, 66.73744938571255, 66.62063664135016, 66.50000133000002, 66.22747284673613, 66.00795379443232, 65.7150802854759, 65.19438696112262, 64.29384856492115, 63.325731662865834, 62.36298129110562, 61.43282457312982, 60.47865161707871, 59.48121189624237, 58.55474098831869, 57.61577695368922, 56.7251829824451, 55.91185527051578, 55.093105509310554, 54.32996519595824, 53.61523753066627, 52.9413, 52.39319476068052, 52.11132084452833, 51.85100518510051, 51.61326903203858, 51.416005141600515, 51.37846165415398, 51.34318973136206, 51.3116, 51.27793076675845, 51.24849487515052, 51.223459021232784, 51.1924153577246, 51.16727441636279, 51.384574307712846, 50.940294905970504, 50.23598995280201, 50.23410000000001], nans=true)
    end

    modes = perplex_query_modes(perplexdir, scratchdir, importas=:Dict)
    @test isa(modes, Dict)
    @test haskey(modes, "Omph(HP)")
    if haskey(modes, "Omph(HP)")
        print("modes[\"Omph(HP)\"]: ")
        println(modes["Omph(HP)"])
        # @test isapprox(modes["Omph(HP)"],[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 1.5736, 7.33183, 13.3273, 13.874, 13.8044, 13.7504, 13.6605, 13.6055, 13.2465, 12.8556, 12.8012, 12.909, 12.8774, 12.8621, 12.8379, 12.8239, 12.8205, 12.839, 12.8654, 12.8914, 12.9423, 13.0084, 13.1195, 13.2487, 13.391, 13.5401, 13.7082, 13.9396, 14.1879, 14.4729, 14.754, 15.0912, 15.5081, 15.9689, 16.4671, 17.0194, 17.5064, 17.1991, 16.9685, 16.6926, 16.4602, 16.1634, 15.921, 15.659, 15.4497, 15.2485, 15.0301, 14.8809, 14.6926, 15.0711, 9.19562, NaN, NaN], nans=true)
    end

    ## --- # # # # # # # # # # # Geothermal gradient example # # # # # # # # # # # #

    # Input parameters
    P_range = (280, 28000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_surf = 273.15 # Temperature of surface (K)
    geotherm = 0.01 # Geothermal gradient of 0.1 K/bar == about 28.4 K/km
    melt_model = ""

    # Configure (run build and vertex)
    @time perplex_configure_geotherm(perplexdir, scratchdir, composition, elements, P_range, T_surf, geotherm;
        dataset="hp02ver.dat",
        excludes=HP_excludes,
        solution_phases=HP_solution_phases,
        npoints=200,
        index=2
    )

    seismic = perplex_query_seismic(perplexdir, scratchdir, index=2)
    @test isa(seismic, Dict)
    @test haskey(seismic, "T(K)")
    @test isa(seismic["T(K)"], Vector{Float64})

    print("seismic[\"T(K)\"]: ")
    println(seismic["T(K)"])

    ## --- # # # # # # # # # # # Pseudosection example # # # # # # # # # # # # #

    P_range = (1000, 5000) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_range = (400+273.15, 600+273.15) # Temperature range to explore, K
    melt_model = ""

    solution_phases = "Opx(HP)\nO(HP)\nF\n"
    excludes = ""

    # Configure (run build and vertex)
    @time perplex_configure_pseudosection(perplexdir, scratchdir, composition,
        elements, P_range, T_range, dataset="hp02ver.dat", excludes=excludes,
        solution_phases=melt_model*solution_phases, index=1, xnodes=50, ynodes=50)

    # Query modes on diagonal line across P-T space
    modes = perplex_query_modes(perplexdir, scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(modes, Dict) && !isempty(modes)
    @test haskey(modes,"T(K)") && all(extrema(modes["T(K)"]) .≈ T_range)

    phase = perplex_query_phase(perplexdir, scratchdir, "Opx(HP)", P_range, T_range, index=1, npoints=200)
    @test isa(phase, Dict) && !isempty(phase)
    @test haskey(phase,"T(K)") && all(extrema(phase["T(K)"]) .≈ T_range)

    system = perplex_query_system(perplexdir, scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(system, Dict) && !isempty(system)
    @test haskey(system,"T(K)") && all(extrema(system["T(K)"]) .≈ T_range)

    # Query seismic properties on diagonal line across P-T space
    seismic = perplex_query_seismic(perplexdir, scratchdir, P_range, T_range, index=1, npoints=200)
    @test isa(seismic, Dict) && !isempty(seismic)
    @test haskey(seismic,"T(K)") && !isempty(seismic["T(K)"]) && all(extrema(seismic["T(K)"]) .≈ T_range)
    @test haskey(seismic, "rho,kg/m3") && !isempty(seismic["rho,kg/m3"]) && !any(x->x<2700, seismic["rho,kg/m3"]) && !any(x->x>3200, seismic["rho,kg/m3"])


    # Query properties on a manually-specified diagonal P-T line
    P = range(first(P_range), last(P_range), length=16)
    T = range(first(T_range), last(T_range), length=16)

    modes = perplex_query_modes(perplexdir, scratchdir, P, T, index=1)
    @test isa(modes, Dict) && !isempty(modes)
    @test haskey(modes,"T(K)") && all(isapprox.(modes["T(K)"], T, atol=0.1))

    phase = perplex_query_phase(perplexdir, scratchdir, "Opx(HP)", P, T, index=1)
    @test isa(phase, Dict) && !isempty(phase)
    @test haskey(phase,"T(K)") && all(isapprox.(phase["T(K)"], T, atol=0.1))

    system = perplex_query_system(perplexdir, scratchdir, P, T, index=1)
    @test isa(system, Dict) && !isempty(system)
    @test haskey(system,"T(K)") && all(isapprox.(system["T(K)"], T, atol=0.1))
    @test haskey(system, "rho,kg/m3") && !any(x->x<2700, system["rho,kg/m3"]) && !any(x->x>3200, system["rho,kg/m3"])

    seismic = perplex_query_seismic(perplexdir, scratchdir, P, T, index=1)
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

## --- End of File
