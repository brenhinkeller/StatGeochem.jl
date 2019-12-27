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
        eHf=(Hf176_Hf177_t ./ CHUR_Hf176_Hf177_t .- 1) .* 10^4

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
            means[:,i] = m
        end

        m = nanmean(means,dim=2)
        el = m - pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) - m

        return (c, m, el, eu)
    end
    export bin_bsr_eHf

## --- Calculate Eu*

    # Full four-element log-linear interpolation, using ionic radii
    function eustar(Nd::Number, Sm::Number, Gd::Number, Tb::Number)
        # Ionic radii, in pm [Tb, Gd, Sm, Nd]
        r = [106.3, 107.8, 109.8, 112.3] # or x = [1, 2, 4, 6]

        # Normalize to chondrite
        y = log.([Tb/0.0374, Gd/0.2055, Sm/0.1530, Nd/0.4670])
        notnan = .~ isnan.(y)

        # Make sure we're interpolating and not extrapolating
        if any(notnan[1:2]) && any(notnan[3:4])
            # Fit a straight line through the chondrite-normalized values
            (a,b) = linreg(r[notnan], y[notnan])
            # De-dormalize output for Eu, interpolating at r = 108.7 pm or x = 3
            eu_interp = 0.0580*exp(a + b*108.7)
        else
            eu_interp = NaN
        end
        return eu_interp
    end

    # Simple geometric mean interpolation from Sm and Gd alone
    function eustar(Sm::Number, Gd::Number)
        # Geometric mean in regular space is equal to the arithmetic mean in log space. Fancy that!
        return 0.0580*sqrt(Sm/0.1530 * Gd/0.2055)
    end

    export eustar

## --- Fe oxide conversions

    function feoconversion(FeO::Number=NaN, Fe2O3::Number=NaN, FeOT::Number=NaN, Fe2O3T::Number=NaN)
        # Compiles data from FeO, Fe2O3, FeOT, and Fe2O3T into
        # a single FeOT value.

        # To convert from Fe2O3 wt % to FeO wt %, multiply by
        conversionfactor = (55.845+15.999) / (55.845+1.5*15.999)

        # If FeOT or Fe2O3T already exists, use that
        if isnan(FeOT)
            if isnan(Fe2O3T)
                FeOT=nansum([Fe2O3*conversionfactor, FeO])
            else
                FeOT=Fe2O3T*conversionfactor
            end
         end

        return FeOT
    end
    export feoconversion

## --- Oxide conversions

    """
    ```julia
    dataset = oxideconversion(dataset::Dict; unitratio::Number=10000)
    ```

    Convert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3)...
    If metals are as PPM, set unitratio=10000 (default); if metals are as wt%,
    set unitratio = 1
    """
    function oxideconversion(dataset::Dict; unitratio::Number=10000)
        # Convert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3)...
        # for i=1:length(source)
        #     conversionfactor(i)=mass.percation.(dest[i])./mass.(source[i]);
        # end

        # Array of elements to convert
        source = ["Si","Ti","Al","Fe","Fe","Mg","Ca","Mn","Na","K","P","Cr","Ni","Co","C","S","H"]
        dest = ["SiO2","TiO2","Al2O3","FeOT","Fe2O3T","MgO","CaO","MnO","Na2O","K2O","P2O5","Cr2O3","NiO","CoO","CO2","SO2","H2O"]
        conversionfactor = [2.13932704290547,1.66847584248889,1.88944149488507,1.28648836426407,1.42973254639611,1.65825961736268,1.39919258253823,1.29121895771597,1.34795912485574,1.20459963614796,2.29133490474735,1.46154369861159,1.27258582901258,1.27147688434143,3.66405794688203,1.99806612601372,8.93601190476191]

        # If source field exists, fill in destination from source
        for i = 1:length(source)
            if haskey(dataset, source[i])
                if ~haskey(dataset, dest[i]) # If destination field doesn't exist, make it.
                    dataset[dest[i]] = fill(NaN, size(dataset[source[i]]))
                end
                t = isnan.(dataset[dest[i]]) .& (.~ isnan.(dataset[source[i]]))
                dataset[dest[i]][t] = dataset[source[i]][t] .* (conversionfactor[i] / unitratio)
            end
        end

        return dataset
    end
    export oxideconversion


## --- Chemical Index of Alteration

    # Chemical Index of Alteration as defined by Nesbitt and Young, 1982
    # Note that CaO should be only igneous CaO excluding any Ca from calcite or apatite
    function CIA(Al2O3::Number, CaO::Number, Na2O::Number, K2O::Number)
        A = Al2O3 / 101.96007714
        C = CaO / 56.0774
        N = Na2O / 61.978538564
        K = K2O / 94.19562
        return A / (A + C + N + K) * 100
    end
    export CIA

    # "Weathering Index of Parker" as defined by Parker, 1970
    function WIP(Na2O::Number, MgO::Number, K2O::Number, CaO::Number)
        Na = Na2O / 30.9895
        Mg = MgO / 40.3044
        K = K2O / 47.0980
        Ca = CaO / 56.0774
        # Denominator for each element is a measure of Nicholls' bond strengths
        return (Na/0.35 + Mg/0.9 + K/0.25 + Ca/0.7) * 100
    end
    export WIP

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
        tsc = []
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
        prefix = joinpath(scratchdir, "out$(index)/")
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
            "ALPHAMELTS_DELTAT	$(trunc(dT,digits=1))\n"  *
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
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

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

## -- Perplex interface: 1. Configuration

    """
    ```julia
    perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        \telements::String=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        \tP_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1;
        \tdataset::String="hp02ver.dat", index::Integer=1, npoints::Integer=100,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n", fluid_eos::Integer=5)
    ```

    Set up a PerpleX calculation for a single bulk composition along a specified
    geothermal gradient and pressure (depth) range. P specified in bar and T_surf
    in Kelvin, with geothermal gradient in units of Kelvin/bar
    """
    function perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1;
        dataset::String="hp02ver.dat", index::Integer=1, npoints::Integer=100,
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", fluid_eos::Integer=5)

        build = joinpath(perplexdir, "build")# path to PerpleX build
        vertex = joinpath(perplexdir, "vertex")# path to PerpleX vertex

        #Configure working directory
        prefix = joinpath(scratchdir, "out$(index)/")
        system("rm -rf $prefix; mkdir -p $prefix")

        # Place required data files
        system("cp $(joinpath(perplexdir,dataset)) $prefix")
        system("cp $(joinpath(perplexdir,"perplex_option.dat")) $prefix")
        system("cp $(joinpath(perplexdir,"solution_model.dat")) $prefix")

        # Edit perplex_option.dat to specify number of nodes at which to solve
        system("sed -e \"s/1d_path .*|/1d_path                   $npoints $npoints |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file.
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(uppercase.(elements) .* "\n")
        write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(P_range[1])\n$(P_range[2])\ny\n") # v6.8.7
        # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n5\n3\nn\ny\n2\n1\n$T_surf\n$geotherm\n$(P_range[1])\n$(P_range[2])\ny\n") # v6.8.1

        # Whole-rock composition
        for i = 1:length(composition)
            write(fp,"$(composition[i]) ")
        end
        # Solution model
        if length(excludes) > 0
            write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\nGeothermal")
        else
            write(fp,"\nn\nn\ny\nsolution_model.dat\n$(solution_phases)\nGeothermal")
        end
        close(fp)

        # build PerpleX problem definition
        system("cd $prefix; $build < build.bat > build.log")

        println("Built problem definition")

        # Run PerpleX vertex calculations
        result = system("cd $prefix; echo $index | $vertex > vertex.log")
        return result
    end
    export perplex_configure_geotherm

    """
    ```julia
    perplex_configure_isobar(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        \telements::String=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"]
        \tP::Number=10000, T::Array{<:Number}=[500+273.15, 1500+273.15];
        \tdataset::String="hp11ver.dat", index::Integer=1, npoints::Integer=100,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n", fluid_eos::Integer=5)
    ```

    Set up a PerpleX calculation for a single bulk composition along a specified
    isobaric temperature gradient. P specified in bar and T_range in Kelvin
    """
    function perplex_configure_isobar(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P::Number=10000, T::Array{<:Number}=[500+273.15, 1500+273.15];
        dataset::String="hp11ver.dat", index::Integer=1, npoints::Integer=100,
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", fluid_eos::Integer=5)

        build = joinpath(perplexdir, "build")# path to PerpleX build
        vertex = joinpath(perplexdir, "vertex")# path to PerpleX vertex

        #Configure working directory
        prefix = joinpath(scratchdir, "out$(index)/")
        system("rm -rf $prefix; mkdir -p $prefix")

        # Place required data files
        system("cp $(joinpath(perplexdir,dataset)) $prefix")
        system("cp $(joinpath(perplexdir,"perplex_option.dat")) $prefix")
        system("cp $(joinpath(perplexdir,"solution_model.dat")) $prefix")

        # Edit perplex_option.dat to specify number of nodes at which to solve
        system("sed -e \"s/1d_path .*|/1d_path                   $npoints $npoints |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file
        # Options based on Perplex v6.8.7
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(uppercase.(elements) .* "\n")
        write(fp,"$index\n$dataset\nperplex_option.dat\nn\n3\nn\nn\nn\n$elementstring\n$fluid_eos\nn\nn\n2\n$(T[1])\n$(T[2])\n$P\ny\n") # v6.8.7
        # write(fp,"$index\n$dataset\nperplex_option.dat\nn\nn\nn\nn\n$elementstring\n$fluid_eos\n3\nn\nn\n2\n$(T[1])\n$(T[2])\n$P\ny\n") # v6.8.1

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
        result = system("cd $prefix; echo $index | $vertex > vertex.log")
        return result
    end
    export perplex_configure_isobar

    """
    ```julia
    perplex_configure_pseudosection(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        \telements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        \tP::Array{<:Number}=[280, 28000], T::Array{<:Number}=[273.15, 1500+273.15];
        \tdataset::String="hp11ver.dat", index::Integer=1, xnodes::Integer=42, ynodes::Integer=42,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n", fluid_eos::Number=5)
    ```

    Set up a PerpleX calculation for a single bulk composition across an entire
    2d P-T space. P specified in bar and T in Kelvin
    """
    function perplex_configure_pseudosection(perplexdir::String, scratchdir::String, composition::Array{<:Number},
        elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
        P::Array{<:Number}=[280, 28000], T::Array{<:Number}=[273.15, 1500+273.15];
        dataset::String="hp11ver.dat", index::Integer=1, xnodes::Integer=42, ynodes::Integer=42,
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n", fluid_eos::Number=5)

        build = joinpath(perplexdir, "build")# path to PerpleX build
        vertex = joinpath(perplexdir, "vertex")# path to PerpleX vertex

        #Configure working directory
        prefix = joinpath(scratchdir, "out$(index)/")
        system("rm -rf $prefix; mkdir -p $prefix")

        # Place required data files
        system("cp $(joinpath(perplexdir,dataset)) $prefix")
        system("cp $(joinpath(perplexdir,"perplex_option.dat")) $prefix")
        system("cp $(joinpath(perplexdir,"solution_model.dat")) $prefix")

        # Edit data files to specify number of nodes at which to solve
        system("sed -e \"s/x_nodes .*|/x_nodes                   $xnodes $xnodes |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/y_nodes .*|/y_nodes                   $ynodes $ynodes |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file
        # Options based on Perplex v6.8.7
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(uppercase.(elements) .* "\n")
        write(fp,"$index\n$dataset\nperplex_option.dat\nn\n2\nn\nn\nn\n$elementstring\n$fluid_eos\nn\n2\n$(T[1])\n$(T[2])\n$(P[1])\n$(P[2])\ny\n") # v6.8.7

        # Whole-rock composition
        for i = 1:length(composition)
            write(fp,"$(composition[i]) ")
        end
        # Solution model
        write(fp,"\nn\ny\nn\n$excludes\ny\nsolution_model.dat\n$solution_phases\nPseudosection")
        close(fp)

        # build PerpleX problem definition
        system("cd $prefix; $build < build.bat > build.log")

        # Run PerpleX vertex calculations
        result = system("cd $prefix; echo $index | $vertex > vertex.log")
        return result
    end
    export perplex_configure_pseudosection

## -- Perplex interface: 2. 0d queries

    """
    ```julia
    perplex_query_point(perplexdir::String, scratchdir::String, indvar::Number; index::Integer=1)
    ```

    Query perplex results at a single temperature on an isobar or single pressure
    on a geotherm. Results are returned as a string.
    """
    function perplex_query_point(perplexdir::String, scratchdir::String, indvar::Number; index::Integer=1)
        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Sanitize T inputs to avoid PerpleX escape sequence
        if indvar == 999
            indvar = 999.001
        end

        # Create werami batch file
        # Options based on Perplex v6.7.2
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n1\n$indvar\n999\n0\n")
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
            @warn "$(prefix)$(index)_1.txt could not be parsed, perplex may not have run"

        end
        return data
    end
    """
    ```julia
    perplex_query_point(perplexdir::String, scratchdir::String, P::Number, T::Number; index::Integer=1)
    ```

    Query perplex results at a single P,T point in a pseudosection.
    Results are returned as a string.
    """
    function perplex_query_point(perplexdir::String, scratchdir::String, P::Number, T::Number; index::Integer=1)
        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Sanitize T inputs to avoid PerpleX escape sequence
        if P == 99
            P = 99.001
        end
        if T == 99
            T = 99.001
        end

        # Create werami batch file
        # Options based on Perplex v6.7.2
        fp = open(prefix*"werami.bat", "w")
        write(fp,"$index\n1\n$T\n$P\n99\n99\n0\n")
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
            @warn "$(prefix)$(index)_1.txt could not be parsed, perplex may not have run"

        end
        return data
    end
    export perplex_query_point

## --- Perplex interface: 3. 1d queries

    # We'll need this for when perplex messes up
    molarmass = Dict("SIO2"=>60.083, "TIO2"=>79.8651, "AL2O3"=>101.96007714, "FE2O3"=>159.6874, "FEO"=>71.8442, "MGO"=>40.304, "CAO"=>56.0774, "MNO"=>70.9370443, "NA2O"=>61.978538564, "K2O"=>94.19562, "H2O"=>18.015, "CO2"=>44.009, "P2O5"=>141.942523997)

    """
    ```julia
    perplex_query_seismic(perplexdir::String, scratchdir::String;
        \tdof::Integer=1, index::Integer=1, include_fluid="n")
    ```

    Query perplex seismic results along a previously configured 1-d path (dof=1,
    isobar or geotherm) or 2-d grid / pseudosection (dof=2).
    Results are returned as a dictionary.
    """
    function perplex_query_seismic(perplexdir::String, scratchdir::String;
        dof::Integer=1, index::Integer=1, include_fluid::String="n")
        # Query a pre-defined path (isobar or geotherm)

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        if dof == 1
            # v6.7.8 1d path
            write(fp,"$index\n3\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\n0\n")
        elseif dof == 2
            # v6.7.8 2d grid
            write(fp,"$index\n2\n2\nn\n$include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\nn\n1\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

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
    """
    ```julia
    perplex_query_seismic(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        \tindex::Integer=1, npoints::Integer=200, include_fluid="n")
    ```

    Query perplex seismic results along a specified P-T path using a pre-computed
    pseudosection. Results are returned as a dictionary.
    """
    function perplex_query_seismic(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        index::Integer=1, npoints::Integer=200, include_fluid="n")
        # Query a new path from a pseudosection

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        # v6.7.8 pseudosection
        write(fp,"$index\n3\nn\n$(T[1])\n$(P[1])\n$(T[2])\n$(P[2])\n$npoints\n2\nn
                $include_fluid\n13\nn\n$include_fluid\n15\nn\n$include_fluid\n0\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

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
    export perplex_query_seismic


    """
    ```julia
    perplex_query_phase(perplexdir::String, scratchdir::String, phase::String;
        \tdof::Integer=1, index::Integer=1, include_fluid="y", clean_units::Bool=true)
    ```

    Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
    along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d
    grid / pseudosection (dof=2). Results are returned as a dictionary.
    """
    function perplex_query_phase(perplexdir::String, scratchdir::String, phase::String;
        dof::Integer=1, index::Integer=1, include_fluid="y", clean_units::Bool=true)
        # Query a pre-defined path (isobar or geotherm)

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        if dof == 1
            # v6.7.8, 1d path
            write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\n0\n")
        elseif dof == 2
            # v6.7.8, 2d grid
            write(fp,"$index\n2\n36\n2\n$phase\n$include_fluid\nn\n1\n0\n") # v6.7.8
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]

            # Renormalize weight percentages
            t = contains.(elements,"wt%")
            total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            # Check if perplex is messing up and outputting mole proportions
            if nanmean(total_weight) < 50
                @warn "Perplex seems to be reporting mole fractions instead of weight percentages, attempting to correct"
                for col = findall(t)
                    data[2:end,col] .*= molarmass[replace(elements[col], ",wt%" => "")]
                end
                total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            end
            data[2:end,t] .*= 100 ./ total_weight

            # Clean up element names
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end

            # Convert to a dictionary
            result = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    """
    ```julia
    perplex_query_phase(perplexdir::String, scratchdir::String, phase::String, P::Array{<:Number}, T::Array{<:Number};
        \tindex::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true)
    ```

    Query all perplex-calculated properties for a specified phase (e.g. "Melt(G)")
    along a specified P-T path using a pre-computed pseudosection. Results are
    returned as a dictionary.
    """
    function perplex_query_phase(perplexdir::String, scratchdir::String, phase::String, P::Array{<:Number}, T::Array{<:Number};
        index::Integer=1, npoints::Integer=200, include_fluid="y", clean_units::Bool=true)
        # Query a new path from a pseudosection

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        # v6.7.8 pseudosection
        write(fp,"$index\n3\nn\n$(T[1])\n$(P[1])\n$(T[2])\n$(P[2])\n$npoints\n36\n2\n$phase\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]

            # Renormalize weight percentages
            t = contains.(elements,"wt%")
            total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            # Check if perplex is messing up and outputting mole proportions
            if nanmean(total_weight) < 50
                @warn "Perplex seems to be reporting mole fractions instead of weight percentages, attempting to correct"
                for col = findall(t)
                    data[2:end,col] .*= molarmass[replace(elements[col], ",wt%" => "")]
                end
                total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            end
            data[2:end,t] .*= 100 ./ total_weight

            # Clean up element names
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end

            # Convert to a dictionary
            result = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    export perplex_query_phase


    """
    ```julia
    perplex_query_modes(perplexdir::String, scratchdir::String;
        \tdof::Integer=1, index::Integer=1, include_fluid="y")
    ```

    Query modal mineralogy (mass proportions) along a previously configured 1-d
    path (dof=1, isobar or geotherm) or 2-d grid / pseudosection (dof=2).
    Results are returned as a dictionary.
    """
    function perplex_query_modes(perplexdir::String, scratchdir::String;
        dof::Integer=1, index::Integer=1, include_fluid="y")
        # Query a pre-defined path (isobar or geotherm)

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        if dof == 1
            # v6.7.8 1d path
            write(fp,"$index\n3\n25\nn\n$include_fluid\n0\n")
        elseif dof == 2
            # v6.7.8 2d grid
            write(fp,"$index\n2\n25\nn\n$include_fluid\nn\n1\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            # Convert to a dictionary
            result = elementify(data)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    """
    ```julia
    perplex_query_modes(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        \tindex::Integer=1, npoints::Integer=200, include_fluid="y")
    ```

    Query modal mineralogy (mass proportions) along a specified P-T path using a
    pre-computed pseudosection. Results are returned as a dictionary.
    """
    function perplex_query_modes(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        index::Integer=1, npoints::Integer=200, include_fluid="y")
        # Query a new path from a pseudosection

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        # v6.7.8 pseudosection
        write(fp,"$index\n3\nn\n$(T[1])\n$(P[1])\n$(T[2])\n$(P[2])\n$npoints\n25\nn\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn"t already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            # Convert to a dictionary
            result = elementify(data)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    export perplex_query_modes


    """
    ```julia
    perplex_query_system(perplexdir::String, scratchdir::String;
        \tindex::Integer=1, include_fluid="y", clean_units::Bool=true)
    ```?

    Query all perplex-calculated properties for the system (with or without fluid)
    along a previously configured 1-d path (dof=1, isobar or geotherm) or 2-d
    grid / pseudosection (dof=2). Results are returned as a dictionary.
    Set include_fluid="n" to return solid+melt only.
    """
    function perplex_query_system(perplexdir::String, scratchdir::String;
        index::Integer=1, include_fluid="y", clean_units::Bool=true)
        # Query a pre-defined path (isobar or geotherm)

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        if dof == 1
            # v6.7.8, 1d path
            write(fp,"$index\n3\n36\n1\n$include_fluid\n0\n")
        elseif dof == 2
            # v6.7.8, 2d grid
            write(fp,"$index\n2\n36\n1\n$include_fluid\nn\n1\n0\n")
        else
            error("Expecting dof = 1 (path) or 2 (grid/pseudosection) degrees of freedom")
        end
        close(fp)

        # Make sure there isn't already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]

            # Renormalize weight percentages
            t = contains.(elements,"wt%")
            total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            data[2:end,t] .*= 100 ./ total_weight

            # Clean up element names
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end

            # Convert to a dictionary
            result = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    """
    ```julia
    function perplex_query_system(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        \tindex::Integer=1, npoints::Integer=200, include_fluid="y",clean_units::Bool=true)
    ```

    Query all perplex-calculated properties for the system (with or without fluid)
    along a specified P-T path using a pre-computed pseudosection. Results are
    returned as a dictionary. Set include_fluid="n" to return solid+melt only.
    """
    function perplex_query_system(perplexdir::String, scratchdir::String, P::Array{<:Number}, T::Array{<:Number};
        index::Integer=1, npoints::Integer=200, include_fluid="y",clean_units::Bool=true)
        # Query a new path from a pseudosection

        werami = joinpath(perplexdir, "werami")# path to PerpleX werami
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Create werami batch file
        fp = open(prefix*"werami.bat", "w")
        # v6.7.8 pseudosection
        write(fp,"$index\n3\nn\n$(T[1])\n$(P[1])\n$(T[2])\n$(P[2])\n$npoints\n36\n1\n$include_fluid\n0\n")
        close(fp)

        # Make sure there isn't already an output
        system("rm -f $(prefix)$(index)_1.tab*")

        # Extract Perplex results with werami
        system("cd $prefix; $werami < werami.bat > werami.log")

        # Ignore initial and trailing whitespace
        system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(index)_1.tab")
        # Merge delimiters
        system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(index)_1.tab")

        # Read results and return them if possible
        result = Dict()
        try
            # Read data as an Array{Any}
            data = readdlm("$(prefix)$(index)_1.tab", ' ', skipstart=8)
            elements = data[1,:]

            # Renormalize weight percentages
            t = contains.(elements,"wt%")
            total_weight = nansum(Float64.(data[2:end,t]),dim=2)
            data[2:end,t] .*= 100 ./ total_weight

            # Clean up element names
            if clean_units
                elements = elements .|> x -> replace(x, ",%" => "_pct") # substutue _pct for ,% in column names
                elements = elements .|> x -> replace(x, ",wt%" => "") # Remove units on major oxides
            end

            # Convert to a dictionary
            result = elementify(data,elements)
        catch
            # Return empty dictionary if file doesn't exist
            @warn "$(prefix)$(index)_1.tab could not be parsed, perplex may not have run"
        end
        return result
    end
    export perplex_query_system


## --- End of File
