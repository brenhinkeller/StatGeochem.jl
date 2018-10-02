################################################################################
# Use the Julia-PerpleX interface to run the same isobaric Perplex calculation
# with different solution models to compare results

## --- Import some useful packages
    using StatGeochem
    using Plots; gr(); default(fmt = :svg);

    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    end

## --- Configure
    # Absolute paths to perplex resources
    perplexdir = "/Users/cbkeller/Applications/perplex-stable/" # Location of executables and solution models to use
    scratchdir = "./scratch/" # Location of directory to store output files
    h0 = plot(xlabel="T (C)", ylabel="Melt percent")

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
    HP_excludes = "";

## --- # # # # # # # # # # melt(G) + G_solution_phases # # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 1
    print("\nmelt(G) + G_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(G)\n"*G_solution_phases, excludes=G_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(G)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(G) + G")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(G) + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_G.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(G) + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_G.pdf")

## --- # # # # # # # # # # melt(G) + W_solution_phases # # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 2
    print("\nmelt(G) + W_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(G)\n"*W_solution_phases, excludes=W_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(G)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(G) + W")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(G) + W_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_G_W.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(G) + W_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_G_W.pdf")


## --- # # # # # # # # # # melt(G) +JH_solution_phases # # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 3
    print("\nmelt(G) + JH_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(G)\n"*JH_solution_phases, excludes=JH_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(G)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(G) + JH")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(G) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_G_JH.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(G) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_G_JH.pdf")

## --- # # # # # # # # # # pMELTS(G) +HP_solution_phases # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 4
    print("\npMELTS(G) + JH_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp02ver.dat", solution_phases="pMELTS(G)\n"*JH_solution_phases, excludes=JH_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "pMELTS(G)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="pMELTS(G) + JH")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="pMELTS(G) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_pMELTS_JH.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="pMELTS(G) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_pMELTS_JH.pdf")

## --- # # # # # # # # # # # melt(W) + W_solution_phases # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 5
    print("\nmelt(W) + W_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(W)\n"*W_solution_phases, excludes=W_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(W)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(W) + W")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(W) + W_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_W.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(W) + W_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_W.pdf")

## --- # # # # # # # # # # melt(W) + G_solution_phases # # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 6
    print("\nmelt(W) + G_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(W)\n"*G_solution_phases, excludes=G_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(W)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(W) + G")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(W) + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_W_G.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(W) + G_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_W_G.pdf")

## --- # # # # # # # # # # # melt(W) +JH_solution_phases # # # # # # # # # # # #

    # Input parameters
    P = 10000 # bar
    T_range = [500+273.15, 1500+273.15]
    idx = 7
    print("\nmelt(W) + JH_solution_phases\n")
    @time perplex_configure_isobaric(perplexdir, scratchdir, composition, elements, P, T_range, dataset="hp11ver.dat", solution_phases="melt(W)\n"*JH_solution_phases, excludes=JH_excludes, index=idx)

    # Query the full isobar -- results returned as elementified dictionary
    T_range_inc = [floor(Int,T_range[1])+1, ceil(Int,T_range[2])-1]
    npoints = T_range_inc[2] - T_range_inc[1] + 1
    bulk = perplex_query_isobar_system(perplexdir, scratchdir, T_range_inc, npoints, index=idx)             # Get system data for all temperatures. Set include_fluid = "n" to get solid+melt only
    modes = perplex_query_isobar_modes(perplexdir, scratchdir, T_range, npoints, index=idx)                 # || phase modes
    melt = perplex_query_isobar_phase(perplexdir, scratchdir, T_range_inc, npoints, "melt(W)", index=idx)  # || melt data

    # Create dictionary to hold solid composition and fill it using what we know from system and melt
    solid = Dict()
    solid["wt_pct"] = 100 - melt["wt_pct"]
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        solid[e] = (bulk[e] - (melt[e] .* melt["wt_pct"]/100)) ./ (solid["wt_pct"]/100)
    end

    # Add results to melt % vs temperature figure
    plot!(h0, melt["T(K)"]-273.15, melt["wt_pct"], label="melt(W) + JH")

    # Plot melt composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in melt", title="melt(W) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], melt[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box)
    savefig(h,"Perplex_MeltTest_W_JH.pdf")

    # Plot solid composition as a function of melt percent
    h = plot(xlabel="Percent melt", ylabel="Wt. % in solid", title="melt(W) + JH_solution_phases, $P bar")
    for e in ["SIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O"]
        plot!(h, melt["wt_pct"], solid[e], label=e)
    end
    plot!(h,fg_color_legend=:white, framestyle=:box, legend=:topleft)
    savefig(h,"Perplex_SolidTest_W_JH.pdf")

## --- Format melt comparison

    plot!(h0,fg_color_legend=:white,legend=:topleft)
    display(h0)
    savefig(h0,"Perplex_T-F_comparison.pdf")

## --- End of File
