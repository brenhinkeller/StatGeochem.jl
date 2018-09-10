# A few PerpleX calculation examples using the Julia-Perplex interface

## --- Import some useful packages

    using StatGeochem
    using Plots

## --- Configure

    # Absolute paths to perplex resources
    perplexdir = "/Users/cbkeller/Applications/perplex-stable/" # Location of executables and solution models to use
    scratchdir = "./scratch/" # Location of directory to store output files

## --- # # # # # # # # # # # # # Initial composition # # # # # # # # # # # # # #

    ## McDonough Pyrolite
    #elements =    [ "SIO2", "TIO2", "AL2O3",  "FEO",  "MNO",  "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    #composition = [45.1242, 0.2005, 4.4623, 8.0723, 0.1354, 37.9043, 3.5598, 0.3610, 0.0291, 0.1511, 0.0440,]

    ## Kelemen (2014) primitive continental basalt. H2O and CO2 are guesses
    #elements =    [ "SIO2", "TIO2", "AL2O3",  "FEO",  "MNO",  "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    #composition = [50.0956, 0.9564, 15.3224, 8.5103, 0.1659, 9.2520, 9.6912, 2.5472, 0.8588, 2.0000, 0.6000,]

    # Kelemen (2014) primitive continental basalt excluding Mn and Ti since most melt models can"t handle them..
    elements =    [ "SIO2", "AL2O3",  "FEO",  "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    composition = [50.0956, 15.3224, 8.5103, 9.2520, 9.6912, 2.5472, 0.8588, 2.0000, 0.6000,]

    ## Average Archean basalt (EarthChem data)
    #elements =    [ "SIO2", "TIO2", "AL2O3",   "FEO",  "MNO",   "MGO",  "CAO", "NA2O",  "K2O",  "H2O",  "CO2",]
    #composition = [49.2054, 0.8401, 12.0551, 11.4018, 0.2198, 12.3997, 9.3113, 1.6549, 0.4630, 1.8935, 0.5555,]

## --- # # # # # # # # # # # Some solution model options # # # # # # # # # # # #

    # Emphasis on phases from Green (2016) -- developed for metabasites, includes what is probably the best (and most expensive) amphibole model. Use with hp11ver.dat
    G_solution_phases = "Augite(G)\nOpx(JH)\ncAmph(G)\noAmph(DP)\nO(JH)\nSp(JH)\nGrt(JH)\nfeldspar_B\nMica(W)\nBio(TCC)\nChl(W)\nCtd(W)\nCrd(W)\nSa(WP)\nSt(W)\nIlm(WPH)\nAtg(PN)\nT\nB\nF\nDo(HP)\nScap\nChum\nNeph(FB)\n"
    G_excludes ="ged\nfanth\ngl\n"

    # Emphasis on phases from White (2014) -- developed for metapelites. Use with hp11ver.dat
    W_solution_phases = "Omph(HP)\nOpx(W)\ncAmph(DP)\noAmph(DP)\nO(JH)\nSp(JH)\nGt(W)\nfeldspar_B\nMica(W)\nBi(W)\nChl(W)\nCtd(W)\nCrd(W)\nSa(WP)\nSt(W) \nIlm(WPH)\nAtg(PN)\nT\nB\nF\nDo(HP)\nScap\nChum\nPu(M)\n"
    W_excludes = "andr\nts\nparg\ngl\nged\nfanth\n"

    # Emphasis on phases from Jennings and Holland (2015) -- developed for mantle melting. Use with hp11ver.dat
    JH_solution_phases = "Cpx(JH)\nOpx(JH)\ncAmph(DP)\noAmph(DP)\nO(JH)\nSp(JH)\nGrt(JH)\nfeldspar_B\nMica(W)\nBio(TCC)\nChl(W)\nCtd(W)\nCrd(W)\nSa(WP)\nSt(W)\nIlm(WPH)\nAtg(PN)\nT\nB\nF\nDo(HP)\nScap\nChum\nNeph(FB)\n"
    JH_excludes = "ts\nparg\ngl\nged\nfanth\n"

    # Emphasis on phases from Holland and Powell -- all phases can be used with hp02ver.dat.
    HP_solution_phases = "Omph(HP)\nOpx(HP)\nGlTrTsPg\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nfeldspar_B\nMica(CF)\nBio(TCC)\nChl(HP)\nCtd(HP)\nSapp(HP)\nSt(HP)\nIlHm(A)\nDo(HP)\nT\nB\nF\n"
    HP_excludes = ""

## --- # # # # # # # # # # # # # Isobaric example # # # # # # # # # # # # # # # #

    # Input parameters
    P = 10000 # Pressure, bar
    T_range = [500+273.15, 1500+273.15] # Temperature range, Kelvin
    melt_model = "melt(G)"

    # Configure (run build and vertex)
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, "hp11ver.dat", melt_model*"\n"*G_solution_phases, G_excludes)

## ---

    # Query all properties at a single temperature -- results returned as text
    T = 1450+273.15
    data_isobaric = perplex_query_isobar(perplexdir, scratchdir, T)
    print(data_isobaric)

## ---

    # Query the full isobar -- results returned as pandas DataFrames
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, melt_model)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

## ---

    # Plot melt composition as a function of melt percent

    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(G) + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    # savefig(h,"MeltComposition.pdf")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Magma SIO2 (wt.%)", ylabel="Wt. % in melt", title="melt(G) + G_solution_phases, $P bar")
    for e in ["AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h,melt["SIO2"], melt[e], label=e)
    end
    # savefig(h,"MeltCompositionvsSiO2.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="$melt_model + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    # savefig(h,"SolidComposition.pdf")

    # Plot modes of all phases as a function of temperature
    h = plot(xlabel="T (C)", ylabel="Weight percent", title("$melt_model + G_solution_phases, $P bar")
    for m in modes["elements"][3:end]
        plot!(h, modes["T(K)"]-273.15, modes[m])
    end
    # savefig(h,"PhaseModes.pdf")

    # Plot modes of all phases as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Weight percent", title="$melt_model + G_solution_phases, $P bar")
    for m in modes["elements"][3:end]
        plot!(h, modes[melt_model], modes[m], label=m)
    end
    # savefig(h,"PhaseModesvsF.pdf")

## --- # # # # # # # # # # # Geothermal gradient example # # # # # # # # # # # #

    # Input parameters
    P_range = [280, 28000] # Pressure range to explore, bar (1-100 km)
    T_surf = 273.15 # Temperature of surface (K)
    geotherm = 0.1 # Geothermal gradient of 0.1 K/bar == about 28.4 K/km
    melt_model = "melt(G)"

    # Configure (run build and vertex)
    @time perplex_configure_geotherm(perplexdir, scratchdir, composition, elements, P_range, T_surf, geotherm, "hp11ver.dat", melt_model*"\n"*G_solution_phases, G_excludes)

    # Query seismic properties along the whole geotherm
    geotherm_sesimic = perplex_query_geotherm_seismic(perplexdir, scratchdir, P_range, 100)

    # Query all properties at a single pressure
    P = 10000
    data_geotherm = perplex_query_geotherm(perplexdir, scratchdir, P)
    print data_geotherm

## --- End of File
