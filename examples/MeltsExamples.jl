## --- Load required packages
    using StatGeochem
    using Plots

## --- # # # # # # # # # # # pMelts equil. batch melting # # # # # # # # # # # #
    # Note: melts_configure requires a working alphamelts installation.
    # If you don't have alphamelts, you'll need to download, install, and
    # configure it. The absolute path to the alphamelts perl script must be input here:
    meltspath = "/usr/local/bin/run_alphamelts.command"

    # The absolute or relative path to a direcory where you want to run MELTS -
    # A number of configuration and output files will be written here.
    scratchdir = "scratch/"

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

## --- Plot melt composition
    h = plot(xlabel="Percent melt",ylabel="Abudance (wt. %) in melt")
    for e in ["SiO2","Al2O3","CaO","MgO","FeO","Na2O","K2O"]
        plot!(h,melt_comp["mass"],melt_comp[e],label=e)
    end
    plot!(h,xlims=(0,100),framestyle=:box,fg_color_legend=:white,legend=:right)
    display(h)

## --- Plot solid composition
    h = plot(xlabel="Percent melt",ylabel="Abudance (wt. %) in solid")
    for e in ["SiO2","Al2O3","CaO","MgO","FeO","Na2O","K2O"]
        plot!(h,100 .- solid_comp["mass"],solid_comp[e],label=e)
    end
    plot!(h,xlims=(0,100),framestyle=:box,fg_color_legend=:white,legend=:right)
    display(h)

## --- Plot phase modes
    h = plot(xlabel="Temperature (C)",ylabel="Abudance (wt. %)")
    for m in modes["elements"][4:end]
        plot!(h,modes["Temperature"],modes[m],label=m)
    end
    plot!(h,ylims=(0,100),framestyle=:box,fg_color_legend=:white,legend=:topleft)
    display(h)

## --- End of File
