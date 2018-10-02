## --- Hafnium isotopes

    # Calculate the initial Hf ratio and epsilon Hf at time t Ma
    function eHf(Hf176_Hf177, Lu176_Hf177, t; eHfOnly::Bool=true)

        # Lutetium decay constant (Soderlund et al., 2004
        lambda = 1.867E-11

        # Present-day CHUR composition (Bouvier et al., 2008)
        CHUR_Hf176_Hf177 = 0.282785
        CHUR_Lu176_Hf177 = 0.0336

        # Calculate initial Hf ratio at time t
        Hf176_Hf177_t = Hf176_Hf177 .- Lu176_Hf177.*(exp.(t .* 10^6*lambda) .- 1)

        # Calculate CHUR Hf ratio at time t
        CHUR_Hf176_Hf177_t = CHUR_Hf176_Hf177 .- CHUR_Lu176_Hf177.*(exp.(t .* 10^6*lambda) .- 1)

        # Calculate corresponding epsilon Hf
        eHf=(Hf176_Hf177_t ./ CHUR_Hf176_Hf177_t .- 1) .* 10^4;

        if eHfOnly
            return eHf
        else
            return (eHf, Hf176_Hf177_t)
        end
    end
    export eHf

    function bin_bsr_eHf(x,Hf176_Hf177,Lu176_Hf177,age,min,max,nbins,x_sigma,Hf176_Hf177_sigma,Lu176_Hf177_sigma,age_sigma,nresamples)
        data = hcat(x,Hf176_Hf177,Lu176_Hf177,age)
        sigma = hcat(x_sigma,Hf176_Hf177_sigma,Lu176_Hf177_sigma,age_sigma)

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(age))
            eHf_resampled = eHf(dbs[:,2], dbs[:,3], dbs[:,4])
            (c,m,s) = binmeans(dbs[:,1], eHf_resampled, min, max, nbins)
            means[:,i] = m;
        end

        m = nanmean(means,dim=2)
        el = m - pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) - m

        return (c, m, el, eu)
    end
    export bin_bsr_eHf

## --- MELTS interface

    # Configure and run MELTS simulation
    function melts_configure(meltspath::String, scratchdir::String, composition::Array{Float64},
        elements::Array, T_range::Array=[1400, 600], P_range::Array=[10000,10000];
        batchstring::String="1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n",
        dT=-10, dP=0, index=1, version="pMELTS",mode="isobaric",fo2path="FMQ",
        fractionatesolids::Bool=false,verbose::Bool=true)

        ############################ Default Settings ###############################
        ##MELTS or pMELTS
        #version = "pMELTS"
        ##Mode (isothermal, isobaric, isentropic, isenthalpic, isochoric, geothermal or PTPath)
        #mode = "isobaric"
        ## Set fO2 constraint, i.e. "IW","COH","FMQ","NNO","HM","None" as a string
        #fo2path = "FMQ"
        ## Fractionate all solids? ("!" for no, "" for yes)
        #fractionatesolids = "!"
        # Mass retained during fractionation
        massin = 0.001
        # Ouptut temperatures in celcius? ("!" for no, "" for yes)
        celciusoutput = ""
        # Save all output? ("!" for no, "" for yes)
        saveall = "!"
        # Fractionate all water? ("!" for no, "" for yes)
        fractionatewater = "!"
        # Fractionate individual phases (specify as strings in cell array, i.e. {"olivine","spinel"})
        fractionate = []
        # Supress individual phases (specify as strings in cell array, i.e. {"leucite"})
        supress = []
        # Coninuous (fractional) melting? ("!" for no, "" for yes)
        continuous = "!"
        # Threshold above which melt is extracted (if fractionation is turned on)
        minf = 0.005
        # Do trace element calculations
        dotrace = "!"
        # Treat water as a trace element
        dotraceh2o = "!"
        # Initial trace compositionT
        tsc = [];
        # Initial trace elements
        telements = []
        # Default global constraints
        Pmax = 90000
        Pmin = 2
        Tmax = 3000
        Tmin = 450
        # Simulation number (for folder, etc)

        ########################## end Default Settings ############################

        # Guess if intention is for calculation to end at Tf or Pf as a min or max
        if T_range[2]<T_range[1]
            Tmin=T_range[2]
        end
        if T_range[2]>T_range[1]
            Tmax=T_range[2]
        end
        if P_range[2]<P_range[1]
            Pmin=P_range[2]
        end
        if P_range[2]>P_range[1]
            Pmax=P_range[2]
        end

        if fractionatesolids
            fractionatesolids = ""
        else
            fractionatesolids = "!"
        end

        # Normalize starting composition
        composition = composition./sum(composition)*100

        # output prefixectory name
        prefix = scratchdir*"out$index/"
        # Ensure directory exists and is empty
        system("rm -rf $prefix; mkdir -p $prefix")

        # Make .melts file containing the starting composition you want to run simulations on
        fp = open(prefix*"sc.melts", "w")
        for i = 1:length(elements)
            write(fp,"Initial Composition: $(elements[i]) $(trunc(composition[i],digits=4))\n")
        end
        for i = 1:length(telements)
            write(fp, "Initial Trace: $(telements[i]) $(trunc(tsc[i],digits=4))\n")
        end

        write(fp, "Initial Temperature: $(trunc(T_range[1],digits=2))\nInitial Pressure: $(trunc(P_range[1],digits=2))\nlog fo2 Path: $fo2path\n")

        for i = 1:length(fractionate)
            write(fp,"Fractionate: $(fractionate[i])\n")
        end
        for i = 1:length(supress)
            write(fp,"Suppress: $(supress[i])\n")
        end

        close(fp)


        # Make melts_env file to specify type of MELTS calculation
        fp = open(prefix*"/melts_env.txt", "w")
        write(fp, "! *************************************\n!  Julia-generated environment file\n! *************************************\n\n"  *
            "! this variable chooses MELTS or pMELTS; for low-pressure use MELTS\n" *
            "ALPHAMELTS_VERSION		$version\n\n" *
            "! do not use this unless fO2 anomalies at the solidus are a problem\n"  *
            "!ALPHAMELTS_ALTERNATIVE_FO2	true\n\n"  *
            "! use this if you want to buffer fO2 for isentropic, isenthalpic or isochoric mode\n! e.g. if you are doing isenthalpic AFC\n"  *
            "!ALPHAMELTS_IMPOSE_FO2		true\n\n"  *
            "! use if you want assimilation and fractional crystallization (AFC)\n"  *
            "!ALPHAMELTS_ASSIMILATE		true\n\n"  *
            "! isothermal, isobaric, isentropic, isenthalpic, isochoric, geothermal or PTPath\n"  *
            "ALPHAMELTS_MODE			$mode\n"  *
            "!ALPHAMELTS_PTPATH_FILE		ptpath.txt\n\n"  *
            "! need to set DELTAP for polybaric paths; DELTAT for isobaric paths\nALPHAMELTS_DELTAP	$(trunc(dP,digits=1))\n"  *
            "ALPHAMELTS_DELTAT	$(trunc(dT,1))\n"  *
            "ALPHAMELTS_MAXP		$(trunc(Pmax,digits=1))\n"  *
            "ALPHAMELTS_MINP		$(trunc(Pmin,digits=1))\n"  *
            "ALPHAMELTS_MAXT		$(trunc(Tmax,digits=1))\n"  *
            "ALPHAMELTS_MINT		$(trunc(Tmin,digits=1))\n\n"  *
            "! this one turns on fractional crystallization for all solids\n! use Fractionate: in the melts file instead for selective fractionation\n"  *
            "$(fractionatesolids)ALPHAMELTS_FRACTIONATE_SOLIDS	true\n"  *
            "$(fractionatesolids)ALPHAMELTS_MASSIN		$massin\n\n"  *
            "! free water is unlikely but can be extracted\n"  *
            "$(fractionatewater)ALPHAMELTS_FRACTIONATE_WATER	true\n"  *
            "$(fractionatewater)ALPHAMELTS_MINW			0.005\n\n"  *
            "! the next one gives an output file that is always updated, even for single calculations\n"  *
            "$(saveall)ALPHAMELTS_SAVE_ALL		true\n"  *
            "!ALPHAMELTS_SKIP_FAILURE		true\n\n"  *
            "! this option converts the output temperature to celcius, like the input\n"  *
            "$(celciusoutput)ALPHAMELTS_CELSIUS_OUTPUT	true\n\n"  *
            "! the next two turn on and off fractional melting\n"  *
            "$(continuous)ALPHAMELTS_CONTINUOUS_MELTING	true\n"  *
            "$(continuous)ALPHAMELTS_MINF			$minf\n"  *
            "$(continuous)ALPHAMELTS_INTEGRATE_FILE	integrate.txt\n\n"  *
            "! the next two options refer to the trace element engine\n"  *
            "$(dotrace)ALPHAMELTS_DO_TRACE		true\n"  *
            "$(dotraceh2o)ALPHAMELTS_DO_TRACE_H2O		true\n")
        close(fp)

        # Make a batch file to run the above .melts file starting from the liquidus
        fp = open(prefix*"/batch.txt", "w")
        write(fp,batchstring)
        close(fp)

        # Run the command
        # Edit the following line(s to make sure you have a correct path to the "run_alphamelts.command" perl script
        if verbose
            system("cd " * prefix  * "; " * meltspath * " -f melts_env.txt -b batch.txt")
        else
            system("cd " * prefix  * "; " * meltspath * " -f melts_env.txt -b batch.txt &>./melts.log")
        end
        return 0
    end
    export melts_configure

    # Get melts results, return as string
    function melts_query(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files

        # Read results and return them if possible
        data = ""
        try
            # Read entire output file as a string
            fp = open(prefix*"Phase_main_tbl.txt", "r")
            data = read(fp,String)
            close(fp)
        catch
            # Return empty string if file doesn't exist
            data = ""
        end
        return data
    end
    export melts_query

    # Get modal phase proportions, return as elementified dictionary
    function melts_query_modes(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"Phase_mass_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_modes

    # Get liquid composition, return as elementified dictionary
    function melts_query_liquid(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"Liquid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_liquid

    # Read solid composition, return as elementified dictionary
    function melts_query_solid(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"Solid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_solid

    # Read system thermodynamic data, return as elementified dictionary
    function melts_query_system(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"System_main_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_system

## -- Perplex interface

    """
    perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::String=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1; dataset::String="hp02ver.dat",
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", index::Int=1)

    Set up a PerpleX calculation for a single bulk composition along a specified
    geothermal gradient and pressure (depth) range. P specified in bar and T_surf
    in Kelvin, with geothermal gradient in units of Kelvin/bar
    """
    function perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1; dataset::String="hp02ver.dat",
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", index::Int=1)

        build = perplexdir * "build" # path to PerpleX build
        vertex = perplexdir * "vertex" # path to PerpleX vertex

        #Configure working directory
        prefix = scratchdir * "out_$index/"
        system("rm -rf $prefix; mkdir -p $prefix")

        # Place required data files
        system("cp $perplexdir$dataset $prefix")
        system("cp $(perplexdir)perplex_option.dat $prefix")
        system("cp $(perplexdir)solution_model.dat $prefix")

        # Create build batch file
        fp = open(prefix*"build.bat", "w")
        # Name, components, and basic options. Holland and Powell (1998) "CORK" fluid equation state.
        elementstring = ""
        for e in elements
            elementstring = elementstring * uppercase(e) * "\n"
        end
        write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n5\n")
        # Pressure gradient details
        write(fp,"3\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(P_range[1])\n$(P_range[2])\ny\n")
        # Whole-rock composition
        for i = 1:length(composition)
            write(fp,"$(composition[i]) ")
        end
        # Solution model
        write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\nGeothermal")
        close(fp)

        # build PerpleX problem definition
        system("cd $prefix; $build < build.bat > build.log")

        # Run PerpleX vertex calculations
        system("cd $prefix; echo $index | $vertex > vertex.log")

        return 0
    end
    export perplex_configure_geotherm

    """
    perplex_configure_isobaric(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::String=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"]
        P::Number=10000, T_range::Array{<:Number}=[500+273.15, 1500+273.15]; dataset::String="hp11ver.dat",
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", index::Int=1)

    Set up a PerpleX calculation for a single bulk composition along a specified
    isobaric temperature gradient. P specified in bar and T_range in Kelvin
    """
    function perplex_configure_isobaric(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P::Number=10000, T_range::Array{<:Number}=[500+273.15, 1500+273.15]; dataset::String="hp11ver.dat",
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", index::Int=1)

        build = perplexdir * "build" # path to PerpleX build
        vertex = perplexdir * "vertex" # path to PerpleX vertex

        #Configure working directory
        prefix = scratchdir * "out_$index/"
        system("rm -rf $prefix; mkdir -p $prefix")

        # Place required data files
        system("cp $perplexdir$dataset $prefix")
        system("cp $(perplexdir)perplex_option.dat $prefix")
        system("cp $(perplexdir)solution_model.dat $prefix")

        # Create build batch file
        fp = open(prefix*"build.bat", "w")
        # Name, components, and basic options. Holland and Powell (1998) "CORK" fluid equation state.
        elementstring = ""
        for e in elements
            elementstring = elementstring * uppercase(e) * "\n"
        end
        write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n5\n")
        # Pressure gradient details
        write(fp,"3\nn\nn\n2\n$(T_range[1])\n$(T_range[2])\n$P\ny\n")
        # Whole-rock composition
        for i = 1:length(composition)
            write(fp,"$(composition[i]) ")
        end
        # Solution model
        write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\nIsobaric")
        close(fp)

        # build PerpleX problem definition
        system("cd $prefix; $build < build.bat > build.log")

        # Run PerpleX vertex calculations
        system("cd $prefix; echo $index | $vertex > vertex.log")

        return 0
    end
    export perplex_configure_isobaric

    # Query perplex results at a single pressure on a geotherm. Results are returned
    # as string read from perplex text file output
    function perplex_query_geotherm(perplexdir::String, scratchdir::String, P::Number; index::Int=1)
        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Sanitize P inputs to avoid PerpleX escape sequence
        if P == 999
            P = 999.001
        end

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n1\n$P\n999\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.txt")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Read results and return them if possible
        data = ""
        try
            # Read entire output file as a string
            fp = open("$(prefix)$(index)_1.txt", "r")
            data = read(fp)
            close(fp)
        catch
            # Return empty string if file doesn't exist
            data = ""
        end
        return data
    end
    export perplex_query_geotherm

    # Query perplex seismic results along a geotherm. Results are returned as
    # a dictionary
    function perplex_query_geotherm_seismic(perplexdir::String, scratchdir::String, P_range::Array{<:Number}=[284.2, 28420], npoints::Int=100; index::Int=1)
        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n3\n1\n$(P_range[1])\n$(P_range[2])\n$npoints\n2\nn\nn\n13\nn\nn\n15\nn\nn\n0\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i 0 $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i 0 $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            # Convert to a dictionary
            data = elementify(data)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export perplex_query_geotherm_seismic

    # Query perplex results at a single temperature on an isobar. Results are
    # returned as string.
    function perplex_query_isobar(perplexdir::String, scratchdir::String, T::Number; index::Int=1)
        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Sanitize T inputs to avoid PerpleX escape sequence
        if T == 999
            T = 999.001
        end

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n1\n$T\n999\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.txt")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Read results and return them if possible
        data = ""
        try
            # Read entire output file as a string
            fp = open("$(prefix)$(index)_1.txt", "r")
            data = read(fp, String)
            close(fp)
        catch
            # Return empty string if file doesn't exist
            data = ""
        end
        return data
    end
    export perplex_query_isobar

    # Query perplex results for a specified phase along an entire isobar.
    # Results are returned as a dictionary
    function perplex_query_isobar_phase(perplexdir::String, scratchdir::String,
        T_range::Array{<:Number}=[773.15,1773.15], npoints::Int=1000, phase="melt(G)"; index::Int=1,
        include_fluid="y", clean_units::Bool=true)

        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n3\n1\n$(T_range[1])\n$(T_range[2])\n$npoints\n36\n2\n$phase\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i 0 $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i 0 $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end
            # Convert to a dictionary
            data = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export perplex_query_isobar_phase

    # Query modal mineralogy along a given isobar. Results are returned as a
    # dictionary
    function perplex_query_isobar_modes(perplexdir::String, scratchdir::String, T_range::Array{<:Number}=[773.15,1773.15], npoints::Int=1000; index::Int=1, include_fluid="y")
        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n3\n1\n$(T_range[1])\n$(T_range[2])\n$npoints\n25\nn\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i 0 $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i 0 $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            # Convert to a dictionary
            data = elementify(data)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export perplex_query_isobar_modes

    # Query calculated system properties along an entire isobar. Results are
    # returned as a dictionary. Set include_fluid = "n" to get solid+melt only.
    function perplex_query_isobar_system(perplexdir::String, scratchdir::String,
        T_range::Array{<:Number}=[773.15,1773.15], npoints::Int=1000; index::Int=1,
        include_fluid="y", clean_units::Bool=true)

        werami = perplexdir * "werami" # path to PerpleX werami
        prefix = scratchdir * "out_$index/" # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n3\n1\n$(T_range[1])\n$(T_range[2])\n$npoints\n36\n1\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i 0 $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i 0 $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        data = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end
            # Convert to a dictionary
            data = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export perplex_query_isobar_system

## --- End of File
