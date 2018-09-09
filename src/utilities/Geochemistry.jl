## --- Hafnium isotopes

    # Calculate the initial Hf ratio and epsilon Hf at time t Ma
    function eHf(Hf176_Hf177, Lu176_Hf177, t; eHfOnly=true)

        # Lutetium decay constant (Soderlund et al., 2004
        lambda = 1.867E-11

        # Present-day CHUR composition (Bouvier et al., 2008)
        CHUR_Hf176_Hf177 = 0.282785
        CHUR_Lu176_Hf177 = 0.0336

        # Calculate initial Hf ratio at time t
        Hf176_Hf177_t = Hf176_Hf177 - Lu176_Hf177.*(exp.(t*10^6*lambda) - 1)

        # Calculate CHUR Hf ratio at time t
        CHUR_Hf176_Hf177_t = CHUR_Hf176_Hf177 - CHUR_Lu176_Hf177.*(exp.(t .* 10^6 .* lambda) - 1);

        # Calculate corresponding epsilon Hf
        eHf=(Hf176_Hf177_t ./ CHUR_Hf176_Hf177_t - 1) * 10^4;

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

        means = Array{Float64}(nbins,nresamples)
        c = Array{Float64}(nbins)
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
        fractionatesolids=false,verbose=true)

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
            write(fp,"Initial Composition: $(elements[i]) $(trunc(composition[i],4))\n")
        end
        for i = 1:length(telements)
            write(fp, "Initial Trace: $(telements[i]) $(trunc(tsc[i],4))\n")
        end

        write(fp, "Initial Temperature: $(trunc(T_range[1],2))\nInitial Pressure: $(trunc(P_range[1],2))\nlog fo2 Path: $fo2path\n")

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
            "! need to set DELTAP for polybaric paths; DELTAT for isobaric paths\nALPHAMELTS_DELTAP	$(trunc(dP,1))\n"  *
            "ALPHAMELTS_DELTAT	$(trunc(dT,1))\n"  *
            "ALPHAMELTS_MAXP		$(trunc(Pmax,1))\n"  *
            "ALPHAMELTS_MINP		$(trunc(Pmin,1))\n"  *
            "ALPHAMELTS_MAXT		$(trunc(Tmax,1))\n"  *
            "ALPHAMELTS_MINT		$(trunc(Tmin,1))\n\n"  *
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
        data = ""

        # Read results and return them if possible
        try
            # Returns results as text string
            fp = open(prefix*"Phase_main_tbl.txt", "r")
            data = read(fp,String)
            close(fp)
        end
        return data
    end
    export melts_query

    # Get modal phase proportions, return as elementified dictionary
    function melts_query_modes(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files
        data = Dict()

        # Read results and return them if possible
        try
            data = readdlm(prefix*"Phase_mass_tbl.txt", ' ', skipstart=1)
            data = elementify(data,floatout=true)
        end
        return data
    end
    export melts_query_modes

    # Get liquid composition, return as elementified dictionary
    function melts_query_liquid(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files
        data = Dict()

        # Read results and return them if possible
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"Liquid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        end
        return data
    end
    export melts_query_liquid

    # Read solid composition, return as elementified dictionary
    function melts_query_solid(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files
        data = Dict()

        # Read results and return them if possible
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"Solid_comp_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        end
        return data
    end
    export melts_query_solid

    # Read system thermodynamic data, return as elementified dictionary
    function melts_query_system(scratchdir::String; index=1)
        prefix = scratchdir*"out$index/" # path to data files
        data = Dict()

        # Read results and return them if possible
        try
            # Read data as an Array{Any}
            data = readdlm(prefix*"System_main_tbl.txt", ' ', skipstart=1)
            # Convert to a dictionary
            data = elementify(data,floatout=true)
        end
        return data
    end
    export melts_query_system

## --- End of File
