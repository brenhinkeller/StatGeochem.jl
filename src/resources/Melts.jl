## --- MELTS interface

    """
    ```julia
    melts_configure(meltspath::String, scratchdir::String, composition::Collection{Number},
        \telements::Collection{String},
        \tT_range=(1400, 600),
        \tP_range=(10000,10000);)
    ```
    Configure and run a MELTS simulation using alphaMELTS.
    Optional keyword arguments and defaults include:

        batchstring::String = "1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n"

    A string defining the sequence of options that would be entered to produce
    the desired calculation if running alphaMELTS at the command line. The
    default string specifies a batch calculation starting at the liquidus.

        dT = -10

    The temperature step, in degrees, between each step of the MELTS calculation

        dP = 0

    The pressure step, in bar, between each step of the MELTS calculation

        index = 1

    An optional variable used to specify a unique suffix for the run directory name

        version::String = "pMELTS"

    A string specifying the desired version of MELTS. Options include `MELTS` and `pMELTS`.

        mode::String = "isobaric"

    A string specifying the desired calculation mode for MELTS. Options include
    `isothermal`, `isobaric`, `isentropic`, `isenthalpic`, `isochoric`,
    `geothermal` and `PTPath`.

        fo2path::String = "FMQ"

    A string specifying the oxygen fugacity buffer to follow, e.g., `FMQ` or `NNO+1`.
    Available buffers include `IW`,`COH`,`FMQ`,`NNO`,`HM`, and `None`

        fractionatesolids::Bool = false

    Fractionate all solids? default is `false`

        suppress::Collection{String} = String[]

    Supress individual phases (specify as strings in array, i.e. `["leucite"]`)

        verbose::Bool = true

    Print verbose MELTS output to terminal (else, write it to `melts.log`)
    """
    function melts_configure(meltspath::String, scratchdir::String, composition::Collection{Number},
        elements::Collection{String}, T_range::Collection{Number}=(1400, 600), P_range::Collection{Number}=(10000,10000);
        batchstring::String="1\nsc.melts\n10\n1\n3\n1\nliquid\n1\n1.0\n0\n10\n0\n4\n0\n",
        dT=-10, dP=0, index=1, version="pMELTS",mode="isobaric",fo2path="FMQ",
        fractionatesolids::Bool=false, suppress::Collection{String}=String[], verbose::Bool=true)

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
        fractionate = String[]
        # Coninuous (fractional) melting? ("!" for no, "" for yes)
        continuous = "!"
        # Threshold above which melt is extracted (if fractionation is turned on)
        minf = 0.005
        # Do trace element calculations
        dotrace = "!"
        # Treat water as a trace element
        dotraceh2o = "!"
        # Initial trace compositionT
        tsc = Float64[]
        # Initial trace elements
        telements = String[]
        # Default global constraints
        Pmax = 90000
        Pmin = 2
        Tmax = 3000
        Tmin = 450
        # Simulation number (for folder, etc)

        ########################## end Default Settings ############################

        # Guess if intention is for calculation to end at Tf or Pf as a min or max
        if last(T_range)<first(T_range)
            Tmin=last(T_range)
        end
        if last(T_range)>first(T_range)
            Tmax=last(T_range)
        end
        if last(P_range)<first(P_range)
            Pmin=last(P_range)
        end
        if last(P_range)>first(P_range)
            Pmax=last(P_range)
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
        for i ∈ eachindex(elements)
            write(fp,"Initial Composition: $(elements[i]) $(trunc(composition[i],digits=4))\n")
        end
        for i ∈ eachindex(telements)
            write(fp, "Initial Trace: $(telements[i]) $(trunc(tsc[i],digits=4))\n")
        end

        write(fp, "Initial Temperature: $(trunc(first(T_range),digits=2))\nInitial Pressure: $(trunc(first(P_range),digits=2))\nlog fo2 Path: $fo2path\n")

        for i ∈ eachindex(fractionate)
            write(fp,"Fractionate: $(fractionate[i])\n")
        end
        for i ∈ eachindex(suppress)
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
    melts_query(scratchdir::String; index=1)
    ```
    Read all phase proportions from `Phase_main_tbl.txt` in specified MELTS run directory
    Returns an elementified dictionary
    """
    function melts_query(scratchdir::String; index=1, importas=:Dict)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        if importas==:Dict
            melts = Dict{String, Union{Vector{String}, Dict}}()
        else
            melts = Dict{String, Union{Vector{String}, NamedTuple}}()
        end
        if isfile(prefix*"/Phase_main_tbl.txt")
            data = readdlm(prefix*"/Phase_main_tbl.txt", ' ', skipblanks=false)
            pos = findall(all(isempty.(data), dims=2) |> vec)
            melts["minerals"] = Array{String}(undef, length(pos)-1)
            for i=1:(length(pos)-1)
                name = data[pos[i]+1,1]
                melts[name] = elementify(data[pos[i]+2:pos[i+1]-1,:], skipnameless=true, importas=importas)
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
    function melts_query_modes(scratchdir::String; index=1, importas=:Dict)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Phase_mass_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Phase_mass_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, standardize=true, skipnameless=true, importas=importas)
        else
            # Return empty dictionary if file doesn't exist
            data = importas==:Dict ? Dict() : ()
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
            data = elementify(data, standardize=true, skipnameless=true, importas=:Dict)

            # Start by transferring over all the non-redundant elements
            modes = typeof(data)()
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
                    Ilmenite = Vector{Float64}(undef, length(t))
                    Magnetite = Vector{Float64}(undef, length(t))
                    if  haskey(melts[m],"MnO")
                        Ilmenite .= (melts[m]["TiO2"] + melts[m]["MnO"]+(melts[m]["TiO2"]*(71.8444/79.8768) - melts[m]["MnO"]*(71.8444/70.9374))) / 100
                        Magnetite .= (melts[m]["FeO"] - (melts[m]["TiO2"])*71.8444/79.8768) * (1+159.6882/71.8444)/100
                    else
                        Ilmenite .= (melts[m]["TiO2"] + melts[m]["TiO2"]*71.8444/79.8768) / 100
                        Magnetite .= (melts[m]["FeO"] - melts[m]["TiO2"]*71.8444/79.8768) * (1+159.6882/71.8444)/100
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
    function melts_query_liquid(scratchdir::String; index=1, importas=:Dict)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Liquid_comp_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Liquid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, standardize=true, skipnameless=true, importas=importas)
        else
            # Return empty dictionary if file doesn't exist
            data = importas==:Dict ? Dict() : ()
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
    function melts_query_solid(scratchdir::String; index=1, importas=:Dict)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/Solid_comp_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"Solid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, standardize=true, skipnameless=true, importas=importas)
        else
            # Return empty dictionary if file doesn't exist
            data = importas==:Dict ? Dict() : ()
        end
        return data
    end
    export melts_query_solid

    """
    ```julia
    melts_query_system(scratchdir::String; index=1, importas=:Dict)
    ```
    Read system thermodynamic and composition data from `System_main_tbl.txt` in
    specified MELTS run directory. Returns an elementified dictionary or tuple.
    """
    function melts_query_system(scratchdir::String; index=1, importas=:Dict)
        prefix = joinpath(scratchdir, "out$(index)/") # path to data files

        # Read results and return them if possible
        if isfile(prefix*"/System_main_tbl.txt")
            # Read data as an Array{Any}
            data = readdlm(prefix*"System_main_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data, standardize=true, skipnameless=true, importas=importas)
        else
            # Return empty dictionary if file doesn't exist
            data = importas==:Dict ? Dict() : ()
        end
        return data
    end
    export melts_query_system

## --- End of File