## --- Calculate Eu*

    """
    ```julia
    eustar(Nd::Number, Sm::Number, Gd::Number, Tb::Number)
    ```
    Calculate expected europium concentration, Eu*, based on abundance of
    adjacent rare earths.
    Full four-element log-linear interpolation, using ionic radii
    """
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

    """
    ```julia
    eustar(Sm::Number, Gd::Number)
    ```
    Calculate expected europium concentration, Eu*, based on abundance of
    adjacent rare earths.
    Simple geometric mean interpolation from Sm and Gd alone
    """
    function eustar(Sm::Number, Gd::Number)
        # Geometric mean in regular space is equal to the arithmetic mean in log space. Fancy that!
        return 0.0580*sqrt(Sm/0.1530 * Gd/0.2055)
    end

    export eustar

## --- Fe oxide conversions

    """
    ```julia
    feoconversion(FeO::Number=NaN, Fe2O3::Number=NaN, FeOT::Number=NaN, Fe2O3T::Number=NaN)
    ```
    Compiles data from FeO, Fe2O3, FeOT, and Fe2O3T into  a single FeOT value.
    """
    function feoconversion(FeO::Number=NaN, Fe2O3::Number=NaN, FeOT::Number=NaN, Fe2O3T::Number=NaN)

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
    Convert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3, ...).

    If metals are as PPM, set unitratio=10000 (default); if metals are as wt%,
    set unitratio = 1
    """
    function oxideconversion(dataset::Dict; unitratio::Number=10000)
        result = copy(dataset)
        return oxideconversion!(result)
    end
    export oxideconversion
    """
    ```julia
    dataset = oxideconversion!(dataset::Dict; unitratio::Number=10000)
    ```
    Convert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3, ...) in place.

    If metals are as PPM, set unitratio=10000 (default); if metals are as wt%,
    set unitratio = 1
    """
    function oxideconversion!(dataset::Dict; unitratio::Number=10000)
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
    export oxideconversion!


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

    """
    ```julia
    melts_configure(meltspath::String, scratchdir::String, composition::Array{Float64},
        \telements::Array,
        \tT_range::Array=[1400, 600],
        \tP_range::Array=[10000,10000];)
    ```
    Configure and run a MELTS simulation using alphaMELTS. Optional keyword arguments and defaults:

    `batchstring::String="1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n"`

    `dT = -10`

    `dP = 0`

    `index = 1`

    `version = "pMELTS"`

    `mode = "isobaric"`

    `fo2path = "FMQ"`
    Oxygen fugacity buffer to follow, e.g., `FMQ` or `NNO+1`

    `fractionatesolids::Bool = false`
    Fractionate all solids

    `suppress::Array{String} = []`
    Supress individual phases (specify as strings in array, i.e. `["leucite"]`)

    `verbose::Bool = true`
    Print verbose MELTS output to terminal (else, write it to `melts.log`)
    """
    function melts_configure(meltspath::String, scratchdir::String, composition::Array{Float64},
        elements::Array, T_range::Array=[1400, 600], P_range::Array=[10000,10000];
        batchstring::String="1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n",
        dT=-10, dP=0, index=1, version="pMELTS",mode="isobaric",fo2path="FMQ",
        fractionatesolids::Bool=false, suppress::Array{String}=String[], verbose::Bool=true)

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
        for i = 1:length(suppress)
            write(fp,"Suppress: $(suppress[i])\n")
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

    """
    ```julia
    melts_query_modes(scratchdir::String; index=1)
    ```
    Read all phase proportions from `Phase_main_tbl.txt` in specified MELTS run directory
    Returns an elementified dictionary
    """
    function melts_query(scratchdir::String; index=1)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        melts = Dict()
        if isfile(prefix*"/Phase_main_tbl.txt")
            data = readdlm(prefix*"/Phase_main_tbl.txt", ' ', skipblanks=false)
            pos = findall(all(isempty.(data), dims=2) |> vec)
            melts["minerals"] = Array{String}(undef, length(pos)-1)
            for i=1:(length(pos)-1)
                name = data[pos[i]+1,1]
                melts[name] = elementify(data[pos[i]+2:pos[i+1]-1,:], skipnameless=true)
                melts["minerals"][i] = name
            end
        end
        return melts
    end
    export melts_query

    """
    ```julia
    melts_query_modes(scratchdir::String; index=1)
    ```
    Read modal phase proportions from `Phase_mass_tbl.txt` in specified MELTS run
    Returns an elementified dictionary
    """
    function melts_query_modes(scratchdir::String; index=1)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Phase_mass_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Phase_mass_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, floatout=true, skipnameless=true)
        else
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_modes

    """
    ```julia
    melts_clean_modes(scratchdir::String; index=1)
    ```
    Read and parse / clean-up modal phase proportions from specified MELTS run directory
    Returns an elementified dictionary
    """
    function melts_clean_modes(scratchdir::String; index=1)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Phase_mass_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Phase_mass_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, floatout=true, skipnameless=true)

            # Start by transferring over all the non-redundant elements
            modes = Dict()
            for e in data["elements"]
                m = replace(e, r"_.*" => s"")
                if haskey(modes, m)
                    modes[m] .+= data[e]
                else
                    modes[m] = copy(data[e])
                end
            end

            # Add the sum of all solids
            modes["solids"] = zeros(size(data["Temperature"]))
            for e in data["elements"][4:end]
                if !contains(e, "water") && !contains(e, "liquid")
                    modes["solids"] .+= data[e]
                end
            end

            # Get full mineral compositions, add feldspar and oxides
            melts = melts_query(scratchdir, index=index)
            if containsi(melts["minerals"],"feldspar")
                modes["anorthite"] = zeros(size(modes["Temperature"]))
                modes["albite"] = zeros(size(modes["Temperature"]))
                modes["orthoclase"] = zeros(size(modes["Temperature"]))
            end
            An_Ca = (238.12507+40.0784) / (15.999+40.0784)
            Ab_Na = (239.22853+22.98977*2) / (15.999+22.98977*2)
            Or_K  = (239.22853+39.09831*2) / (15.999+39.09831*2)
            if containsi(melts["minerals"],"rhm_oxide")
                modes["ilmenite"] = zeros(size(modes["Temperature"]))
                modes["magnetite"] = zeros(size(modes["Temperature"]))
                modes["hematite"] = zeros(size(modes["Temperature"]))
            end
            for m in melts["minerals"]
                if containsi(m,"feldspar")
                    t = vec(findclosest(melts[m]["Temperature"],modes["Temperature"]))
                    AnAbOr = [melts[m]["CaO"]*An_Ca melts[m]["Na2O"]*Ab_Na melts[m]["K2O"]*Or_K] |> x -> x ./ sum(x, dims=2)
                    modes["anorthite"][t] .+= AnAbOr[:,1] .*  melts[m]["mass"]
                    modes["albite"][t] .+= AnAbOr[:,2] .*  melts[m]["mass"]
                    modes["orthoclase"][t] .+= AnAbOr[:,3] .*  melts[m]["mass"]
                elseif containsi(m,"rhm_oxide")
                    t = vec(findclosest(melts[m]["Temperature"],modes["Temperature"]))
                    if  haskey(melts[m],"MnO")
                        Ilmenite = (melts[m]["TiO2"] + melts[m]["MnO"]+(melts[m]["TiO2"]*(71.8444/79.8768) - melts[m]["MnO"]*(71.8444/70.9374))) / 100
                        Magnetite = (melts[m]["FeO"] - (melts[m]["TiO2"])*71.8444/79.8768) * (1+159.6882/71.8444)/100
                    else
                        Ilmenite = (melts[m]["TiO2"] + melts[m]["TiO2"]*71.8444/79.8768) / 100
                        Magnetite = (melts[m]["FeO"] - melts[m]["TiO2"]*71.8444/79.8768) * (1+159.6882/71.8444)/100
                    end
                    Magnetite[Magnetite.<0] .= 0
                    Hematite = (melts[m]["Fe2O3"] - Magnetite*100*159.6882/231.5326)/100
                    modes["ilmenite"][t] .+= melts[m]["mass"] .* Ilmenite
                    modes["magnetite"][t] .+= melts[m]["mass"] .* Magnetite
                    modes["hematite"][t] .+= melts[m]["mass"] .* Hematite
                end
            end
            minerals = sort(collect(keys(modes)))
            modes["elements"] = ["Pressure","Temperature","mass","solids","liquid"] ∪ minerals[.!containsi.(minerals, "feldspar") .& .!containsi.(minerals, "rhm")]
        else
            # Return empty dictionary if file doesn't exist
            modes = Dict()
        end
        return modes
    end
    export melts_clean_modes

    """
    ```julia
    melts_query_liquid(scratchdir::String; index=1)
    ```
    Read liquid composition from `Liquid_comp_tbl.txt` in specified MELTS run directory
    Returns an elementified dictionary
    """
    function melts_query_liquid(scratchdir::String; index=1)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Liquid_comp_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Liquid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, floatout=true, skipnameless=true)
        else
            # Return empty dictionary if file doesn't exist
            data = Dict()
        end
        return data
    end
    export melts_query_liquid

    """
    ```julia
    melts_query_solid(scratchdir::String; index=1)
    ```
    Read solid composition from `Solid_comp_tbl.txt` in specified MELTS run directory
    Returns an elementified dictionary
    """
    function melts_query_solid(scratchdir::String; index=1)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Solid_comp_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Solid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, floatout=true, skipnameless=true)
        else
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
        if isfile(prefix*"/System_main_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"System_main_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, floatout=true, skipnameless=true)
        else
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
        \tdataset::String="hp02ver.dat",
        \tindex::Integer=1,
        \tnpoints::Integer=100,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n",
        \tmode_basis::String="vol",  #["vol", "wt", "mol"]
        \tcomposition_basis::String="wt",  #["vol", "wt", "mol"]
        \tfluid_eos::Integer=5)
    ```

    Set up a PerpleX calculation for a single bulk composition along a specified
    geothermal gradient and pressure (depth) range. P specified in bar and T_surf
    in Kelvin, with geothermal gradient in units of Kelvin/bar
    """
    function perplex_configure_geotherm(perplexdir::String, scratchdir::String, composition::Array{<:Number},
            elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
            P_range::Array{<:Number}=[280,28000], T_surf::Number=273.15, geotherm::Number=0.1;
            dataset::String="hp02ver.dat",
            index::Integer=1,
            npoints::Integer=100,
            solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
            excludes::String="ts\nparg\ngl\nged\nfanth\ng\n",
            mode_basis::String="vol",
            composition_basis::String="wt",
            fluid_eos::Integer=5
        )

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

        # Edit perplex_option.dat to output all seismic properties 
        #println("editing perplex options ")
        system("sed -e \"s/seismic_output .*|/seismic_output                   all |/\" -i.backup $(prefix)perplex_option.dat")

        # Specify whether we want volume or weight percentages
        system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file.
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(elements .* "\n")
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
        \tdataset::String="hp11ver.dat",
        \tindex::Integer=1,
        \tnpoints::Integer=100,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n",
        \tmode_basis::String="vol",  #["vol", "wt", "mol"]
        \tcomposition_basis::String="wt",  #["vol", "wt", "mol"]
        \tfluid_eos::Integer=5)
    ```

    Set up a PerpleX calculation for a single bulk composition along a specified
    isobaric temperature gradient. P specified in bar and T_range in Kelvin
    """
    function perplex_configure_isobar(perplexdir::String, scratchdir::String, composition::Array{<:Number},
            elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
            P::Number=10000, T::Array{<:Number}=[500+273.15, 1500+273.15];
            dataset::String="hp11ver.dat",
            index::Integer=1,
            npoints::Integer=100,
            solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
            excludes::String="ts\nparg\ngl\nged\nfanth\ng\n",
            mode_basis::String="vol",
            composition_basis::String="wt",
            fluid_eos::Integer=5
        )

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

        # Specify whether we want volume or weight percentages
        system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file
        # Options based on Perplex v6.8.7
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(elements .* "\n")
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
        \tdataset::String="hp11ver.dat",
        \tindex::Integer=1,
        \txnodes::Integer=42,
        \tynodes::Integer=42,
        \tsolution_phases::String="O(HP)\\nOpx(HP)\\nOmph(GHP)\\nGt(HP)\\noAmph(DP)\\ncAmph(DP)\\nT\\nB\\nChl(HP)\\nBio(TCC)\\nMica(CF)\\nCtd(HP)\\nIlHm(A)\\nSp(HP)\\nSapp(HP)\\nSt(HP)\\nfeldspar_B\\nDo(HP)\\nF\\n",
        \texcludes::String="ts\\nparg\\ngl\\nged\\nfanth\\ng\\n",
        \tmode_basis::String="vol", #["vol", "wt", "mol"]
        \tcomposition_basis::String="wt", #["wt", "mol"]
        \tfluid_eos::Number=5)
    ```

    Set up a PerpleX calculation for a single bulk composition across an entire
    2d P-T space. P specified in bar and T in Kelvin.
    """
    function perplex_configure_pseudosection(perplexdir::String, scratchdir::String, composition::Array{<:Number},
            elements::Array{String}=["SIO2","TIO2","AL2O3","FEO","MGO","CAO","NA2O","K2O","H2O"],
            P::Array{<:Number}=[280, 28000], T::Array{<:Number}=[273.15, 1500+273.15];
            dataset::String="hp11ver.dat",
            index::Integer=1,
            xnodes::Integer=42,
            ynodes::Integer=42,
            solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
            excludes::String="ts\nparg\ngl\nged\nfanth\ng\n",
            mode_basis::String="vol",
            composition_basis::String="wt",
            fluid_eos::Number=5
        )

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

        # Specify whether we want volume or weight percentages
        system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
        system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

        # Create build batch file
        # Options based on Perplex v6.8.7
        fp = open(prefix*"build.bat", "w")

        # Name, components, and basic options. P-T conditions.
        # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
        elementstring = join(elements .* "\n")
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

    # # We'll need this for when perplex messes up
    # molarmass = Dict("SIO2"=>60.083, "TIO2"=>79.8651, "AL2O3"=>101.96007714, "FE2O3"=>159.6874, "FEO"=>71.8442, "MGO"=>40.304, "CAO"=>56.0774, "MNO"=>70.9370443, "NA2O"=>61.978538564, "K2O"=>94.19562, "H2O"=>18.015, "CO2"=>44.009, "P2O5"=>141.942523997)

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
            write(fp,"$index\n3\n36\n2\n$phase\n$include_fluid\n5\n0\n")
            # If a named phase (e.g. feldspar) has multiple immiscible phases, average them (5)
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
                @warn "Perplex seems to be reporting mole fractions instead of weight percentages"
                # Attempt to change back to weight percentages
                # for col = findall(t)
                #     data[2:end,col] .*= molarmass[replace(elements[col], ",wt%" => "")]
                # end
                # total_weight = nansum(Float64.(data[2:end,t]),dim=2)
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
                @warn "Perplex seems to be reporting mole fractions instead of weight percentages"
                # , attempting to correct
                # for col = findall(t)
                #     data[2:end,col] .*= molarmass[replace(elements[col], ",wt%" => "")]
                # end
                # total_weight = nansum(Float64.(data[2:end,t]),dim=2)
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

    Currently returns vol % 
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
            # Convert to a dictionary. 
            # Perplex sometimes returns duplicates of a single solution model, sum them.
            result = elementify(data, sumduplicates=true)
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
            result = elementify(data, sumduplicates=true)
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
        index::Integer=1, include_fluid="y", clean_units::Bool=true, dof::Integer=1)
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

    # Translate between perplex names and germ names
    function germ_perplex_name_matches(germ_name, perplex_name)
        # Feldspar
        if germ_name == "Albite"
            any(perplex_name .== ["ab", "abh"])
        elseif germ_name == "Anorthite"
            perplex_name == "an"
        elseif germ_name == "Orthoclase"
            any(perplex_name .== ["mic", "Kf", "San", "San(TH)"])
        # Amphibole
        elseif germ_name == "Amphibole"
            any(lowercase(perplex_name) .== ["gl", "fgl", "rieb", "anth", "fanth", "cumm", "grun", "tr", "ftr", "ged", "parg", "ts"]) ||
            any(contains.(perplex_name, ["Amph", "GlTrTs", "Act(", "Anth"]))
        # Mica
        elseif germ_name == "Biotite"
            any(perplex_name .== ["ann"]) ||
            any(contains.(perplex_name, ["Bi(", "Bio("]))
        elseif germ_name == "Phlogopite"
            any(lowercase(perplex_name) .== ["naph", "phl"])
        # Pyroxene
        elseif germ_name == "Clinopyroxene"
            any(lowercase(perplex_name) .== ["di", "hed", "acm", "jd"]) ||
            any(contains.(perplex_name, ["Augite", "Cpx", "Omph"]))
        elseif germ_name == "Orthopyroxene"
            any(lowercase(perplex_name) .== ["en", "fs"]) ||
            contains(perplex_name, "Opx")
        # Cordierite
        elseif germ_name == "Cordierite"
            any(lowercase(perplex_name) .== ["crd", "fcrd", "hcrd", "mncrd"]) ||
            contains(perplex_name, "Crd")
        # Garnet
        elseif germ_name == "Garnet"
            any(lowercase(perplex_name) .== ["py", "spss", "alm", "andr", "gr"]) ||
            any(contains.(perplex_name, ["Grt", "Gt(", "Maj"]))
        # Oxides
        elseif germ_name == "Ilmenite"
            perplex_name == "ilm" || any(contains.(perplex_name, ["Ilm", "IlHm", "IlGk"]))
        elseif germ_name == "Magnetite"
            perplex_name == "mt"
        elseif germ_name == "Rutile"
            perplex_name == "ru"
        # Feldspathoids
        elseif germ_name == "Leucite"
            perplex_name == "lc"
        elseif germ_name == "Nepheline"
            perplex_name == "ne" || contains(perplex_name, "Neph")
        # Olivine
        elseif germ_name == "Olivine"
            any(lowercase(perplex_name) .== ["fo", "fa"]) ||
            any(contains.(perplex_name, ["O(", "Ol("]))
        # Spinel
        elseif germ_name == "Spinel"
            any(lowercase(perplex_name) .== ["sp", "usp"]) ||
            contains(perplex_name, "Sp(")
        # Accessories
        elseif germ_name == "Sphene"
            perplex_name == "sph"
        elseif germ_name == "Zircon"
            perplex_name == "zrc"
        elseif germ_name == "Baddeleyite"
            perplex_name == "bdy"
        else
            false
        end
    end
    export germ_perplex_name_matches

    function perplex_phase_is_fluid(phase_name)
        any(phase_name .== ["F", "WADDAH", "H2O"]) ||
        any(contains.(phase_name, ["Aq_", "F(", "Fluid"]))
    end
    export perplex_phase_is_fluid

    function perplex_phase_is_melt(phase_name)
        any(phase_name .== ["h2oL", "abL", "anL", "diL", "enL", "faL", "kspL", "qL", "silL"]) ||
        any(contains.(phase_name, ["liq", "melt", "LIQ", "MELTS"]))
    end
    export perplex_phase_is_melt

    function perplex_phase_is_solid(phase_name)
        !perplex_phase_is_fluid(phase_name) && !perplex_phase_is_melt(phase_name) &&
        !any(contains.(phase_name, ["P(", "T(", "Pressure", "Temperature", "elements", "minerals"]))
    end
    export perplex_phase_is_solid

    function perplex_expand_name(name)
        abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm","oilm","lime", "mft", "mt", "mang", "oxide", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi",)
        full_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote(ordered)", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite-bredigite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermaks pyroxene", "Diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "mg-tschermak", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite(Na-free)", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "fe-carpholite", "fe-sapphirine(793)", "mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite(14Ang)", "clinochlore(ordered)", "daphnite", "Fe-sudoite", "Mn-chlorite", "Sudoite", "antigorite", "chrysotile", "Fe-talc", "Kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)","lime", "magnesioferrite", "magnetite", "manganosite", "nickel", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "Forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica",)
        t = name .== abbreviations
        if any(t)
            full_names[findfirst(t)]
        else
            name
        end
    end
    export perplex_expand_name

    function perplex_abbreviate_name(name)
        abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm", "oilm", "lime", "mft", "mt", "mang", "oxide", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi",)
        full_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote(ordered)", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite-bredigite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermaks pyroxene", "Diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "mg-tschermak", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite(Na-free)", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "fe-carpholite", "fe-sapphirine(793)", "mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite(14Ang)", "clinochlore(ordered)", "daphnite", "Fe-sudoite", "Mn-chlorite", "Sudoite", "antigorite", "chrysotile", "Fe-talc", "Kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)", "lime", "magnesioferrite", "magnetite", "manganosite", "nickel", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "Forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica",)
        t = name .== full_names
        if any(t)
            abbreviations[findfirst(t)]
        else
            name
        end
    end
    export perplex_abbreviate_name

    function perplex_common_name(name)
        abbreviations = ("ak", "alm", "and", "andr", "chum", "cz", "crd", "ep", "fa", "fctd", "fcrd", "fep", "fosm", "fst", "fo", "geh", "gr", "hcrd", "tpz", "ky", "larn", "law", "merw", "mctd", "mst", "mnctd", "mncrd", "mnst", "mont", "osm1", "osm2", "phA", "pump", "py", "rnk", "sill", "spss", "sph", "spu", "teph", "ty", "vsv", "zrc", "zo", "acm", "cats", "di", "en", "fs", "hed", "jd", "mgts", "pswo", "pxmn", "rhod", "wo", "anth", "cumm", "fanth", "fgl", "ftr", "ged", "gl", "grun", "parg", "rieb", "tr", "ts", "deer", "fcar", "fspr", "mcar", "spr4", "spr7", "ann", "cel", "east", "fcel", "ma", "mnbi", "mu", "naph", "pa", "phl", "afchl", "ames", "clin", "daph", "fsud", "mnchl", "sud", "atg", "chr", "fta", "kao", "pre", "prl", "ta", "tats", "ab", "anl", "an", "coe", "crst", "heu", "abh", "kals", "lmt", "lc", "me", "mic", "ne", "q", "san", "stlb", "stv", "trd", "wrk", "bdy", "cor", "geik", "hem", "herc", "ilm", "oilm", "lime", "mft", "mt", "mang", "oxide", "per", "pnt", "ru", "sp", "usp", "br", "dsp", "gth", "ank", "arag", "cc", "dol", "mag", "rhc", "sid", "diam", "gph", "iron", "Ni", "CO2", "CO", "H2", "CH4", "O2", "H2O", "abL", "anL", "diL", "enL", "faL", "fliq", "foL", "h2oL", "hliq", "kspL", "mliq", "qL", "silL", "H+", "Cl-", "OH-", "Na+", "K+", "Ca++", "Mg++", "Fe++", "Al+++", "CO3", "AlOH3", "AlOH4-", "KOH", "HCL", "KCL", "NaCl", "CaCl2", "CaCl+", "MgCl2", "MgCl", "FeCl2", "aqSi", "Augite(G)", "Cpx(JH)", "Cpx(l)", "Cpx(h)", "Cpx(stx)", "Cpx(stx7)", "Omph(HP)", "Cpx(HP)", "Cpx(m)", "Cpx(stx8)", "Omph(GHP)", "cAmph(G)", "Cumm", "Gl", "Tr", "GlTrTsPg", "Amph(DHP)", "Amph(DPW)", "Ca-Amph(D)", "Na-Amph(D)", "Act(M)", "GlTrTsMr", "cAmph(DP)", "melt(G)", "melt(W)", "melt(HP)", "pMELTS(G)", "mMELTS(G)", "LIQ(NK)", "LIQ(EF)", "Chl(W)", "Chl(HP)", "Chl(LWV)", "O(JH)", "O(SG)", "O(HP)", "O(HPK)", "O(stx)", "O(stx7)", "Ol(m)", "O(stx8)", "Sp(JH)", "GaHcSp", "Sp(JR)", "Sp(GS)", "Sp(HP)", "Sp(stx)", "CrSp", "Sp(stx7)", "Sp(WPC)", "Sp(stx8)", "Pl(JH)", "Pl(h)", "Pl(stx8)", "Kf", "San", "San(TH)", "Grt(JH)", "Gt(W)", "CrGt", "Gt(MPF)", "Gt(B)", "Gt(GCT)", "Gt(HP)", "Gt(EWHP)", "Gt(WPH)", "Gt(stx)", "Gt(stx8)", "Gt(WPPH)", "ZrGt(KP)", "Maj", "Opx(JH)", "Opx(W)", "Opx(HP)", "CrOpx(HP)", "Opx(stx)", "Opx(stx8)", "Mica(W)", "Pheng(HP)", "MaPa", "Mica(CF)", "Mica(CHA1)", "Mica(CHA)", "Mica+(CHA)", "Mica(M)", "Mica(SGH)", "Ctd(W)", "Ctd(HP)", "Ctd(SGH)", "St(W)", "St(HP)", "Bi(W)", "Bio(TCC)", "Bio(WPH)", "Bio(HP)", "Crd(W)", "hCrd", "Sa(WP)", "Sapp(HP)", "Sapp(KWP)", "Sapp(TP)", "Osm(HP)", "F", "F(salt)", "COH-Fluid", "Aq_solven0", "WADDAH", "T", "Scap", "Carp", "Carp(M)", "Carp(SGH)", "Sud(Livi)", "Sud", "Sud(M)", "Anth", "o-Amph", "oAmph(DP)", "feldspar", "feldspar_B", "Pl(I1,HP)", "Fsp(C1)", "Do(HP)", "M(HP)", "Do(AE)", "Cc(AE)", "oCcM(HP)", "Carb(M)", "oCcM(EF)", "dis(EF)", "IlHm(A)", "IlGkPy", "Ilm(WPH)", "Ilm(WPH0)", "Neph(FB)", "Chum", "Atg(PN)", "B", "Pu(M)", "Stlp(M)", "Wus",)
        common_names = ("akermanite", "almandine", "andalusite", "andradite", "clinohumite", "clinozoisite", "cordierite", "epidote", "fayalite", "Fe-chloritoid", "Fe-cordierite", "Fe-epidote", "Fe-osumilite", "Fe-staurolite", "forsterite", "gehlenite", "grossular", "hydrous cordierite", "hydroxy-topaz", "kyanite", "larnite", "lawsonite", "merwinite", "Mg-chloritoid", "Mg-staurolite", "Mn-chloritoid", "Mn-cordierite", "Mn-staurolite", "monticellite", "osumilite(1)", "osumilite(2)", "phase A", "pumpellyite", "pyrope", "rankinite", "sillimanite", "spessartine", "sphene", "spurrite", "tephroite", "tilleyite", "vesuvianite", "zircon", "zoisite", "acmite", "Ca-tschermakite", "diopside", "enstatite", "ferrosilite", "hedenbergite", "jadeite", "Mg-tschermakite", "pseudowollastonite", "pyroxmangite", "rhodonite", "wollastonite", "anthophyllite", "cummingtonite", "Fe-anthophyllite", "Fe-glaucophane", "ferroactinolite", "gedrite", "glaucophane", "grunerite", "pargasite", "riebeckite", "tremolite", "tschermakite", "deerite", "Fe-carpholite", "Fe-sapphirine(793)", "Mg-carpholite", "sapphirine(442)", "sapphirine(793)", "annite", "celadonite", "eastonite", "Fe-celadonite", "margarite", "Mn-biotite", "muscovite", "Na-phlogopite", "paragonite", "phlogopite", "Al-free chlorite", "amesite", "clinochlore", "daphnite", "Fe-sudoite", "Mn-chlorite", "sudoite", "antigorite", "chrysotile", "Fe-talc", "kaolinite", "prehnite", "pyrophyllite", "talc", "tschermak-talc", "albite", "analcite", "anorthite", "coesite", "cristobalite", "heulandite", "highalbite", "kalsilite", "laumontite", "leucite", "meionite", "microcline", "nepheline", "quartz", "sanidine", "stilbite", "stishovite", "tridymite", "wairakite", "baddeleyite", "corundum", "geikielite", "hematite", "hercynite", "ilmenite", "ilmenite(ordered)", "lime", "magnesioferrite", "magnetite", "manganosite", "nickel", "periclase", "pyrophanite", "rutile", "spinel", "ulvospinel", "brucite", "diaspore", "goethite", "ankerite", "aragonite", "calcite", "dolomite", "magnesite", "rhodochrosite", "siderite", "diamond", "graphite", "iron", "nickel", "carbon dioxide", "carbon monoxide", "hydrogen", "methane", "oxygen", "water fluid", "albite liquid", "anorthite liquid", "diopside liquid", "enstatite liquid", "fayalite liquid", "Fe-liquid (in KFMASH)", "forsterite liquid", "H2O liquid", "H2O liquid (in KFMASH)", "K-feldspar liquid", "Mg liquid (in KFMASH)", "Silica liquid", "Sillimanite liquid", "H+(aq)", "Cl(aq)", "OH(aq)", "Na+(aq)", "K+(aq)", "Ca2+(aq)", "Mg2+(aq)", "Fe2+(aq)", "Al3+(aq)", "CO3--(aq)", "Al(OH)3(aq)", "Al(OH)4----(aq)", "KOH(aq)", "HCl(aq)", "KCl(aq)", "NaCl(aq)", "CaCl(aq)", "CaCl+(aq)", "MgCl2(aq)", "MgCl+(aq)", "FeCl(aq)", "Aqueous silica", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinopyroxene", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "clinoamphibole", "melt", "melt", "melt", "melt", "melt", "melt", "melt", "chlorite", "chlorite", "chlorite", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "olivine", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "spinel", "plagioclase", "plagioclase", "plagioclase", "k-feldspar", "k-feldspar", "k-feldspar", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "garnet", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "orthopyroxene", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "white mica", "chloritoid", "chloritoid", "chloritoid", "staurolite", "staurolite", "biotite", "biotite", "biotite", "biotite", "cordierite", "cordierite", "sapphirine", "sapphirine", "sapphirine", "sapphirine", "osumilite", "fluid", "fluid", "fluid", "fluid", "fluid", "talc", "scapolite", "carpholite", "carpholite", "carpholite", "sudoite", "sudoite", "sudoite", "orthoamphibole", "orthoamphibole", "orthoamphibole", "ternary feldspar", "ternary feldspar", "ternary feldspar", "ternary feldspar", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "calcite", "ilmenite", "ilmenite", "ilmenite", "ilmenite", "nepheline", "clinohumite", "serpentine", "brucite", "pumpellyite", "stilpnomelane", "wüstite",)
        t = name .== abbreviations
        if any(t)
            common_names[findfirst(t)]
        elseif contains(name,"anorthite")
            "anorthite"
        elseif contains(name,"albite")
            "albite"
        elseif contains(name,"orthoclase")
            "orthoclase"
        else
            name
        end
    end
    export perplex_common_name


## -- Zircon saturation calculations

    """
    ```julia
    M = tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate zircon saturation M-value based on major element concentrations
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    """
    function tzircM(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MnO::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/55.8667
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation ratios
        normconst = nansum([Na, K, Ca, Al, Si, Ti, Fe, Mg, Mn, P])
        K = K / normconst
        Na = Na / normconst
        Ca = Ca / normconst
        Al = Al / normconst
        Si = Si / normconst

        return (Na + K + 2*Ca)/(Al * Si)
    end
    function tzircM(SiO2::AbstractArray, TiO2::AbstractArray, Al2O3::AbstractArray, FeOT::AbstractArray, MnO::AbstractArray, MgO::AbstractArray, CaO::AbstractArray, Na2O::AbstractArray, K2O::AbstractArray, P2O5::AbstractArray)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/55.8667
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation ratios
        normconst = nansum([Na K Ca Al Si Ti Fe Mg Mn P], dim=2)
        K .= K ./ normconst
        Na .= Na ./ normconst
        Ca .= Ca ./ normconst
        Al .= Al ./ normconst
        Si .= Si ./ normconst

        return (Na + K + 2*Ca)./(Al .* Si)
    end
    export tzircM

    """
    ```julia
    ZrSat = tzircZr(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
    ```
    Calculate zircon saturation Zr concentration for a given temperature (in C)
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    """
    function tzircZr(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
        M = tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        # Boehnke, Watson, et al., 2013
        Zr = @. max(496000. /(exp(10108. /(T+273.15) -0.32 -1.16*M)), 0)
        return Zr
    end
    export tzircZr

    """
    ```julia
    T = tzirc(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, Zr)
    ```
    Calculate zircon saturation temperature in degrees Celsius
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    """
    function tzirc(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, Zr)
        M = tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        # Boehnke, Watson, et al., 2013
        T = @. 10108. / (0.32 + 1.16*M + log(496000. / Zr)) - 273.15
        return T
    end
    export tzirc


    function tspheneC(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MnO::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/55.8667
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation ratios
        normconst = nansum([Na, K, Ca, Al, Si, Ti, Fe, Mg, Mn, P])
        K = K / normconst
        Na = Na / normconst
        Ca = Ca / normconst
        Al = Al / normconst
        Si = Si / normconst

        eCa = Ca - Al/2 + Na/2 + K/2
        C = (10 * eCa) / (Al * Si)
    end
    function tspheneC(SiO2::AbstractArray, TiO2::AbstractArray, Al2O3::AbstractArray, FeOT::AbstractArray, MnO::AbstractArray, MgO::AbstractArray, CaO::AbstractArray, Na2O::AbstractArray, K2O::AbstractArray, P2O5::AbstractArray)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/55.8667
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation ratios
        normconst = nansum([Na K Ca Al Si Ti Fe Mg Mn P], dim=2)
        K .= K ./ normconst
        Na .= Na ./ normconst
        Ca .= Ca ./ normconst
        Al .= Al ./ normconst
        Si .= Si ./ normconst

        eCa = Ca - Al/2 + Na/2 + K/2
        C = (10 * eCa) ./ (Al .* Si)
    end

    """
    ```julia
    TiSat = tspheneTi(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
    ```
    Calculate sphene saturation Ti concentration for a given temperature (in C)
    Following the sphene saturation calibration of Ayers et al., 2018
    (10.1130/abs/2018AM-320568)
    """
    function tspheneTi(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
        C = tspheneC(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        TiO2 = @. max(0.79*C - 7993/(T+273.15) + 7.88, 0)
    end
    export tspheneTi


    """
    ```julia
    T = tsphene(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate sphene saturation temperature in degrees Celsius
    Following the sphene saturation calibration of Ayers et al., 2018
    (10.1130/abs/2018AM-320568)
    """
    function tsphene(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        C = tspheneC(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        T = @. 7993/(0.79*C - TiO2 + 7.88) - 273.15
    end
    export tsphene


## --- End of File
