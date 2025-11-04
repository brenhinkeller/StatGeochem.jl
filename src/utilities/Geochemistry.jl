## --- Calculate Eu*

    """
    ```julia
    eustar(Nd::Number, Sm::Number, Gd::Number, Tb::Number)
    ```
    Calculate expected europium concentration, Eu*, based on abundance of
    adjacent rare earths.

    Full four-element log-linear interpolation, assuming 3+ ionic radii and the
    chondritic abundances of Sun and McDonough 1989 (doi: 10.1144/gsl.sp.1989.042.01.19)
    """
    function eustar(Nd::Number, Sm::Number, Gd::Number, Tb::Number)
        # Ionic radii, in pm [Tb, Gd, Sm, Nd]
        r = [106.3, 107.8, 109.8, 112.3] # or x = [1, 2, 4, 6]

        # Normalize to chondrite
        y = log.([Tb/0.0374, Gd/0.2055, Sm/0.1530, Nd/0.4670])
        notnan = .!isnan.(y)

        # Make sure we're interpolating and not extrapolating
        if any(view(notnan, 1:2)) && any(view(notnan, 3:4))
            # Fit a straight line through the chondrite-normalized values
            x = r[notnan]
            (a,b) = hcat(fill!(similar(x), 1), x) \ y[notnan]
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

    Simple geometric mean interpolation from Sm and Gd alone, assuming the chondritic 
    abundances of Sun and McDonough 1989 (doi: 10.1144/gsl.sp.1989.042.01.19), that is
    Eu* = `0.0580*sqrt(Sm/0.1530 * Gd/0.2055)`
    """
    function eustar(Sm::Number, Gd::Number)
        # Geometric mean in regular space is equal to the arithmetic mean in log space. Fancy that!
        return 0.0580*sqrt(Sm/0.1530 * Gd/0.2055)
    end

    export eustar

## --- CIPW norm

    """
    ```julia
    cipw_norm(SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Returns
    ```
    quartz, orthoclase, plagioclase, corundum, nepheline, diopside, orthopyroxene, olivine, magnetite, ilmenite, apatite
    ```
    """
    function cipw_norm(SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O, K2O, P2O5)
        SiO2  /= 60.0843
        TiO2  /= 79.8988
        Al2O3 /= 101.9613
        Fe2O3 /= 159.6922
        FeO  /= 71.8464
        MnO  /= 70.9374
        MgO  /= 40.3044
        CaO  /= 56.0794
        Na2O /= 61.9789
        K2O  /= 94.1960
        P2O5 /= 141.9445

        FeO = nanadd(FeO, MnO)
        CaO -= 3.333333333333333 * P2O5
        apatite = 0.6666666666666666 * P2O5
        # P2O5 = 0
        FeO -= TiO2
        ilmenite = TiO2
        FeO -= Fe2O3
        magnetite = Fe2O3
        # Fe2O3 = 0
        Al2O3 -= K2O
        orthoclase = K2O
        # K2O = 0
        Al2O3 -= Na2O
        albite = Na2O
        if CaO > Al2O3
            CaO -= Al2O3
            anorthite = Al2O3
            Al2O3 = 0
        else
            Al2O3 -= CaO
            anorthite = CaO
            CaO = 0
        end
        if Al2O3 > 0
            corundum = Al2O3
            Al2O3 = 0
        else
            corundum = 0
        end
        Mg′ = MgO / (MgO + FeO)
        FMO = FeO + MgO
        FMO_weight = (Mg′*40.3044)+((1-Mg′)*71.8464)
        if CaO > 0
            FMO -= CaO
            diopside = CaO
        else
            diopside = 0
        end
        orthopyroxene = FMO
        pSi1 = 6orthoclase + 6albite + 2anorthite + 2diopside + orthopyroxene
        if pSi1 < SiO2
            quartz = SiO2 - pSi1
            nepheline = 0
            olivine = 0
        else
            quartz = 0
            pSi2 = 6orthoclase + 6albite + 2anorthite + 2diopside
            pSi3 = SiO2 - pSi2
            if FMO > 2pSi3
                orthopyroxene = 0
                olivine = FMO
                FMO = 0
                pSi4 = 6orthoclase + 2anorthite + 2diopside + 0.5olivine
                pSi5 = SiO2 - pSi4
                Albite = (pSi5-(2*Na2O))/4
                nepheline = Na2O-Albite
            else
                nepheline = 0
                orthopyroxene = 2pSi3 - FMO
                olivine = FMO - pSi3
            end
        end
        orthoclase *= 2
        nepheline *= 2
        albite *= 2
        An′ = anorthite/(anorthite+albite)
        plag_weight = (An′*278.2093)+((1-An′)*262.2230)
        plagioclase = albite+anorthite

        quartz *= 60.0843
        orthoclase *= 278.3315
        plagioclase *= plag_weight
        corundum *= 101.9613
        nepheline *= 142.0544
        diopside *= (172.248 + FMO_weight)
        orthopyroxene *= (60.0843 + FMO_weight)
        olivine *= (60.0843 + 2FMO_weight)
        magnetite *= 231.5386
        ilmenite *= 151.7452
        apatite *= 504.3152

        return (quartz=quartz, orthoclase=orthoclase, plagioclase=plagioclase,
            corundum=corundum, nepheline=nepheline, diopside=diopside,
            orthopyroxene=orthopyroxene, olivine=olivine, magnetite=magnetite,
            ilmenite=ilmenite, apatite=apatite)
    end
    # export cipw_norm

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
                if isnan(Fe2O3)
                    FeOT = FeO
                elseif isnan(FeO)
                    FeOT = Fe2O3*conversionfactor
                else
                    FeOT = Fe2O3*conversionfactor + FeO
                end
            else
                FeOT=Fe2O3T*conversionfactor
            end
         end

        return FeOT
    end
    export feoconversion

## --- Oxide conversions

    function fillifnan!(dest::AbstractArray, source::AbstractArray)
        @inbounds for i in eachindex(dest, source)
            if isnan(dest[i]) && !isnan(source[i])
                dest[i] = source[i]
            end
        end
        return dest
    end
    function fillifnan!(dest::AbstractArray, source::AbstractArray, factor::Number)
        @inbounds for i in eachindex(dest, source)
            if isnan(dest[i]) && !isnan(source[i])
                dest[i] = source[i] * factor
            end
        end
        return dest
    end

    function nannegative!(a::AbstractArray)
        @inbounds for i in eachindex(a)
            if a[i] < 0
                a[i] = NaN
            end
        end
        return a
    end

    """
    ```julia
    converted_dataset = oxideconversion(dataset::Union{Dict,NamedTuple}; unitratio::Number=10000)
    ```
    As `oxideconversion!`, but returning a copy rather than modifying in-place
    """
    oxideconversion(ds::Union{Dict,NamedTuple}; kwargs...) = oxideconversion!(deepcopy(ds); kwargs...)
    export oxideconversion
    
    """
    ```julia
    oxideconversion!(dataset::Dict; unitratio::Number=10000)
    ```
    Convert major elements (Ti, Al, etc.) into corresponding oxides (TiO2, Al2O3, ...) in place if extant.

    If metals are expected as PPM, set unitratio=10000 (default); if metals are as wt%,
    set unitratio = 1

    See also `oxideconversion`, c.f. `metalconversion!`
    """
    function oxideconversion!(dataset::NamedTuple; unitratio::Number=10000)
        # List of elements to convert
        source = (:Si, :Ti, :Al, :Fe, :Fe, :Mg, :Ca, :Mn, :Li, :Na, :K, :P, :Cr, :Ni, :Co, :S, :H)
        dest = (:SiO2, :TiO2, :Al2O3, :FeOT, :Fe2O3T, :MgO, :CaO, :MnO, :Li2O, :Na2O, :K2O, :P2O5, :Cr2O3, :NiO, :CoO, :SO3, :H2O)
        conversionfactor = (2.13932704290547,1.66847584248889,1.88944149488507,1.28648836426407,1.42973254639611,1.65825961736268,1.39919258253823,1.29121895771597,2.1526657060518732,1.34795912485574,1.20459963614796,2.29133490474735,1.46154369861159,1.27258582901258,1.27147688434143,2.4970991890205863,8.93601190476191)
        @assert eachindex(source) == eachindex(dest) == eachindex(conversionfactor)

        # If source field exists, fill in destination from source
        for i ∈ eachindex(source)
            if haskey(dataset, source[i])
                if haskey(dataset, dest[i]) # If destination field doesn't exist, make it.
                    oxide, metal = dataset[dest[i]], dataset[source[i]]
                    fillifnan!(oxide, metal, conversionfactor[i]/unitratio)
                end
            end
        end
        return dataset
    end
    oxideconversion!(ds::Dict; kwargs...) = (oxideconversion!(TupleDataset(ds); kwargs...); ds)
    export oxideconversion!


    """
    ```julia
    converted_dataset = metalconversion(dataset::Union{Dict,NamedTuple}; unitratio::Number=10000)
    ```
    As `metalconversion!`, but returning a copy rather than modifying in-place
    """
    metalconversion(ds::Union{Dict,NamedTuple}; kwargs...) = metalconversion!(copy(ds); kwargs...)
    export metalconversion

    """
    ```julia
    dataset = metalconversion!(dataset::Union{Dict,NamedTuple}; unitratio::Number=10000)
    ```
    Convert minor element oxides (MnO, Cr2O3, NiO, ...) into corresponding metals (Mn, Cr, Ni, ...) in place if extant.

    If metals are expected as parts per million (ppm), set unitratio=10000 (default); if metals are as wt%, set unitratio = 1

    See also `metalconversion`, c.f. `oxideconversion!`
    """
    function metalconversion!(dataset::NamedTuple; unitratio::Number=10000)
        # List of elements to convert
        dest = (:Mn, :P, :Cr, :Ni, :Co, :Sr, :Ba, :Li, :S,)
        source = (:MnO, :P2O5, :Cr2O3, :NiO, :CoO, :SrO, :BaO, :Li2O, :SO3)
        conversionfactor = (0.7744619872751028, 0.4364268173666496, 0.6842080746199798, 0.785801615263874, 0.786486968277016, 0.8455993051534453, 0.8956541815613328, 0.46454031259412965, 0.4004646689233921)

        # If source field exists, fill in destination from source
        for i ∈ eachindex(source)
            if haskey(dataset, source[i])
                if haskey(dataset, dest[i]) # If destination field doesn't exist, make it.
                    metal, oxide = dataset[dest[i]], dataset[source[i]]
                    fillifnan!(metal, oxide, conversionfactor[i]*unitratio)
                end
            end
        end
        return dataset
    end
    metalconversion!(ds::Dict; kwargs...) = (metalconversion!(TupleDataset(ds); kwargs...); ds)
    export metalconversion!



    """
    ```julia
    carbonateconversion!(dataset::NamedTuple)
    ```
    Convert carbonates (CaCO3, MgCO3) into corresponding metal oxides and CO2 if extant, in place,
    as well as synchonizing TIC, TOC, TC, C and CO2. All are assumed to be reported in the same units,
    (likely wt. %) except for C, which is assumed to be equivalent to unitratio * TC, 
    """
    function carbonateconversion!(ds::NamedTuple; unitratio=10000)
        # Calculate CO2 if both CaCO3 and MgCO3 are reported
        if haskey(ds, :CaCO3) && haskey(ds, :MgCO3) && haskey(ds, :CO2)
            fillifnan!(ds.CO2, ds.CaCO3*0.43971009048182363 .+ ds.MgCO3*0.5219717006867268)
        end

        # Populate oxides and CO2 from carbonates and TIC
        source = (:CaCO3, :CaCO3, :MgCO3, :MgCO3, :TIC,)
        dest = (:CaO, :CO2, :MgO, :CO2, :CO2)
        conversionfactor = (0.5602899095181764, 0.43971009048182363, 0.4780282993132732, 0.5219717006867268, 3.664057946882025)
        for i in eachindex(source)
            if haskey(ds, source[i])
                if haskey(ds, dest[i])
                    d, s = ds[dest[i]], ds[source[i]]
                    fillifnan!(d, s, conversionfactor[i])
                end
            end
        end

        # Fill TC from C and TIC from CO2
        if haskey(ds,:TC) && haskey(ds, :C)
            fillifnan!(ds.TC, ds.C, 1e-4)
        end
        if haskey(ds,:TIC) && haskey(ds, :CO2)
            fillifnan!(ds.TIC, ds.CO2, 0.27292144788565975)
        end

        # Synchronise TOC, TIC, TC
        if haskey(ds, :TC) && haskey(ds, :TOC) && haskey(ds, :TIC)
            fillifnan!(ds.TC, ds.TOC + ds.TIC)
            fillifnan!(ds.TOC, ds.TC - ds.TIC)
            nannegative!(ds.TOC)
            fillifnan!(ds.TIC, ds.TC - ds.TOC)
            nannegative!(ds.TIC)
            if haskey(ds, :CO2)
                # If we have new TIC values, fill CO2 again
                fillifnan!(ds.CO2, ds.TIC, 3.664057946882025)
            end
        end

        # Fill C from any available source
        if haskey(ds,:TC) && haskey(ds, :C)
            fillifnan!(ds.C, ds.TC, 1e4)
        end
        if haskey(ds,:TOC) && haskey(ds,:TIC) && haskey(ds, :C)
            fillifnan!(ds.C, ds.TOC + ds.TIC, 1e4)
        end
        if haskey(ds,:TOC) && haskey(ds,:CO2) && haskey(ds, :C)
            fillifnan!(ds.C, ds.TOC + ds.CO2/3.664057946882025, 1e4)
        end
        if haskey(ds,:TOC) && haskey(ds, :C)
            fillifnan!(ds.C, ds.TOC, 1e4)
        end
        if haskey(ds,:TIC) && haskey(ds, :C)
            fillifnan!(ds.C, ds.TIC, 1e4)
        end
        if haskey(ds,:CO2) && haskey(ds, :C)
            fillifnan!(ds.C, ds.CO2, 1e4/3.664057946882025)
        end

        return ds
    end
    carbonateconversion!(ds::Dict) = (carbonateconversion!(TupleDataset(ds)); ds)
    export carbonateconversion!

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

## -- Zircon saturation calculations

    """
    ```julia
    M = Boehnke_tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate zircon saturation M-value based on major element concentrations
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    (doi: 10.1016/j.chemgeo.2013.05.028)
    """
    function Boehnke_tzircM(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MnO::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation fractions
        normconst = nansum((Na, K, Ca, Al, Si, Ti, Fe, Mg, Mn, P))
        K, Na, Ca, Al, Si = (K, Na, Ca, Al, Si) ./ normconst

        M = (Na + K + 2*Ca)/(Al * Si)
        return M
    end
    function Boehnke_tzircM(SiO2::AbstractArray, TiO2::AbstractArray, Al2O3::AbstractArray, FeOT::AbstractArray, MnO::AbstractArray, MgO::AbstractArray, CaO::AbstractArray, Na2O::AbstractArray, K2O::AbstractArray, P2O5::AbstractArray)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation fractions
        normconst = nansum([Na K Ca Al Si Ti Fe Mg Mn P], dim=2)
        K .= K ./ normconst
        Na .= Na ./ normconst
        Ca .= Ca ./ normconst
        Al .= Al ./ normconst
        Si .= Si ./ normconst

        M = (Na + K + 2*Ca)./(Al .* Si)
        return M
    end
    export Boehnke_tzircM
    tzircM = Boehnke_tzircM
    export tzircM

    """
    ```julia
    ZrSat = Boehnke_tzircZr(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
    ```
    Calculate zircon saturation Zr concentration for a given temperature (in C)
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    (doi: 10.1016/j.chemgeo.2013.05.028)
    """
    function Boehnke_tzircZr(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
        M = Boehnke_tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        # Boehnke, Watson, et al., 2013
        ZrSat = @. max(496000. /(exp(10108. /(T+273.15) -0.32 -1.16*M)), 0)
        return ZrSat
    end
    export Boehnke_tzircZr
    tzircZr = Boehnke_tzircZr
    export tzircZr

    """
    ```julia
    T = Boehnke_tzirc(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, Zr)
    ```
    Calculate zircon saturation temperature in degrees Celsius
    Following the zircon saturation calibration of Boehnke, Watson, et al., 2013
    (doi: 10.1016/j.chemgeo.2013.05.028)
    """
    function Boehnke_tzirc(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, Zr)
        M = Boehnke_tzircM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        # Boehnke, Watson, et al., 2013
        TC = @. 10108. / (0.32 + 1.16*M + log(496000. / Zr)) - 273.15
        return TC
    end
    export Boehnke_tzirc
    tzirc = Boehnke_tzirc
    export tzirc


## --- Sphene saturation calculations

    function Ayers_tspheneC(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MnO::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

        # Normalize cation fractions
        normconst = nansum((Na, K, Ca, Al, Si, Ti, Fe, Mg, Mn, P))
        K, Na, Ca, Al, Si = (K, Na, Ca, Al, Si) ./ normconst

        eCa = Ca - Al/2 + Na/2 + K/2
        return (10 * eCa) / (Al * Si)
    end

    function Ayers_tspheneM(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MnO::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        Mn = MnO/70.9374
        P = P2O5/70.9723

         # Normalize cation fractions
        normconst = nansum((Na, K, Ca, Al, Si, Ti, Fe, Mg, Mn, P))
        K, Na, Ca, Al, Si = (K, Na, Ca, Al, Si) ./ normconst

        M = (Na + K + (2 * Ca))/(Al * Si)
        return M
    end

     """
    ```julia
    TiO2Sat = Ayers_tspheneTiO2(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
    ```
    Calculate sphene saturation TiO2 concentration (in wt. %) for a given temperature
    (in C) following the sphene saturation calibration of Ayers et al., 2022
    (doi: 10.1007/s00410-022-01902-z)
    """
    function Ayers_tspheneTiO2(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, TC)
        M = Ayers_tspheneM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        TiO2 = max((0.978*M)+(0.0048*(TC+273.15))-5.90, 0)
        return TiO2
    end
    export Ayers_tspheneTiO2

    """
    ```julia
    TC = Ayers_tsphene(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate sphene saturation temperature in degrees Celsius
    Following the sphene saturation calibration of Ayers et al., 2022
    (doi: 10.1007/s00410-022-01902-z)
    """
    function Ayers_tsphene(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        M = Ayers_tspheneM(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        TC = (TiO2 + 5.90 - (0.978*M))/0.0048 - 273.15
        return TC
    end
    export Ayers_tsphene

    """
    ```julia
    TiO2Sat = Ayers_tspheneTiO2_18(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, T)
    ```
    Calculate sphene saturation TiO2 concentration (in wt. %) for a given temperature
    (in C) following the sphene saturation calibration of Ayers et al., 2018
    (doi: 10.1130/abs/2018AM-320568)
    """
    function Ayers_tspheneTiO2_18(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5, TC)
        C = Ayers_tspheneC(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        TiO2 = max(0.79*C - 7993/(TC+273.15) + 7.88, 0)
        return TiO2
    end
    export Ayers_tspheneTiO2_18


    """
    ```julia
    TC = Ayers_tsphene_18(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate sphene saturation temperature in degrees Celsius
    Following the sphene saturation calibration of Ayers et al., 2018
    (doi: 10.1130/abs/2018AM-320568)
    """
    function Ayers_tsphene_18(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        C = Ayers_tspheneC(SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O, P2O5)
        TC = 7993/(0.79*C - TiO2 + 7.88) - 273.15
        return TC
    end
    export Ayers_tsphene_18

## --- Rutile saturation calculations

    function FM(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        #Cations
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        P = P2O5/70.9723

        # Normalize cation fractions
        normconst = nansum((Na, K, Ca, Al, Si, Ti, Fe, Mg, P))
        Na, K, Ca, Mg, Fe, Al, Si = (Na, K, Ca, Mg, Fe, Al, Si) ./ normconst

        return (Na + K + 2(Ca + Mg + Fe)) / (Al * Si)
    end


    """
    ```julia
    TC = Hayden_trutile(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate rutile saturation temperature in degrees Celcius
    following the rutile saturation model of Hayden and Watson, 2007
    (doi: 10.1016/j.epsl.2007.04.020)
    """
    function Hayden_trutile(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number)
        TK = 5305.0 / (7.95 - log10(TiO2 * 10000 * 47.867/(47.867+15.999*2)) + 0.124*FM(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5))
        TC = TK - 273.15
    end

    """
    ```julia
    TC = Hayden_trutile(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5, TC)
    ```
    Calculate the TiO2 concentration in weight percent required for rutile
    saturation at temperature `TC` degrees Celcius, following the rutile
    saturation model of Hayden and Watson, 2007
    (doi: 10.1016/j.epsl.2007.04.020)
    """
    function Hayden_trutileTiO2(SiO2::Number, TiO2::Number, Al2O3::Number, FeOT::Number, MgO::Number, CaO::Number, Na2O::Number, K2O::Number, P2O5::Number, TC::Number)
        TK = TC + 273.15
        return exp10(7.95 - 5315.0/TK + 0.124*FM(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5)) * (47.867+15.999*2)/47.867 * 1e-5
    end


## --- Monazite and xenotime saturation calculations

    """
    ```julia
    LREEmolwt(La, Ce, Pr, Nd, Sm, Gd)
    ```
    Returns the average molecular weight of the LREE considered in the
    REEt value from the monazite saturation model of Montel 1993
    (doi: 10.1016/0009-2541(93)90250-M)
    """
    function LREEmolwt(La, Ce, Pr, Nd, Sm, Gd)
        # All as PPM
        nansum((138.905477La, 140.1161Ce, 140.907662Pr, 144.2423Nd, 150.362Sm, 157.253Gd)) / nansum((La, Ce, Pr, Nd, Sm, Gd))
    end

    """
    ```julia
    LREEt(La, Ce, Pr, Nd, Sm, Gd)
    ```
    Returns the sum of the LREE concentrations divided by their respective molar masses.
    If REE are input in parts per million by weight (ppmw), the result is in units of moles per megagram.
    This is equivalent to the REEt value from the monazite saturation model of Montel 1993
    (doi: 10.1016/0009-2541(93)90250-M)
    """
    function LREEt(La::T, Ce::T, Pr::T, Nd::T, Sm::T, Gd::T) where T <: Number
        # All as PPM
        nansum((La/138.905477, Ce/140.1161, Pr/140.907662, Nd/144.2423, Sm/150.362, Gd/157.253))
    end


    function Montel_tmonaziteD(SiO2::T, TiO2::T, Al2O3::T, FeOT::T, MgO::T, CaO::T, Na2O::T, K2O::T, Li2O::T, H2O::T) where T <: Number
        #Cations
        H = H2O/9.0075
        Li = Li2O/14.9395
        Na = Na2O/30.9895
        K = K2O/47.0827
        Ca = CaO/56.0774
        Al = Al2O3/50.9806
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        # Anions
        # O = 0.5H + 0.5Li + 0.5Na + 0.5K + Ca + Mg + Fe + 1.5Al + 2Si + 2Ti

        # Calculate cation fractions
        normconst = nansum((H, Li, Na, K, Ca, Al, Si, Ti, Fe, Mg))
        Li, Na, K, Ca, Al, Si = (Li, Na, K, Ca, Al, Si) ./ normconst

        # Note that the paper incorrectly describes this equation as being
        # written in terms of atomic percent ("at.%"), but in fact it appears
        # to require molar cation fractions, just as does the analagous M-value
        # equation found in zircon saturation papers
        D = (Na + K + Li + 2Ca) / (Al * (Al + Si))
        return D
    end

    """
    ```julia
    REEt = Montel_tmonaziteREE(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, T)
    ```
    Calculate monazite saturation REEt value (in [ppm/mol.wt.]) for a given
    temperature (in C) following the monazite saturation model of Montel 1993
    (doi: 10.1016/0009-2541(93)90250-M), where:

    D = (Na + K + Li + 2Ca) / Al * 1/(Al + Si)) # all as molar cation fractions (not at. %!)
    ln(REEt) = 9.50 + 2.34D + 0.3879√H2O - 13318/T # H2O as wt.%
    REEt = Σ REEᵢ(ppm) / at. weight (g/mol)
    """
    function Montel_tmonaziteREE(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, T)
        D = Montel_tmonaziteD(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O) # input as wt. %
        REEt = exp(9.50 + 2.34*D + 0.3879*sqrt(H2O) - 13318/(T+272.15))
        return REEt
    end
    export Montel_tmonaziteREE

    """
    ```julia
    TC = Montel_tmonazite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, La, Ce, Pr, Nd, Sm, Gd)
    ```
    Calculate monazite saturation temperature in degrees Celcius
    following the monazite saturation model of Montel 1993
    (doi: 10.1016/0009-2541(93)90250-M)
    """
    function Montel_tmonazite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O, La, Ce, Pr, Nd, Sm, Gd)
        D = Montel_tmonaziteD(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, Li2O, H2O) # input as wt. %
        REEt = LREEt(La, Ce, Pr, Nd, Sm, Gd) # input in PPM
        TC = 13318/(9.50 + 2.34*D + 0.3879*sqrt(H2O) - log(REEt)) - 272.15
        return TC
    end
    export Montel_tmonazite


    """
    ```julia
    LREEt = Rusiecka_tmonaziteREE(P_ppm, TC)
    ```
    Calculate the LREEt (mol/Megagram) value required for monazite saturation at a
    temperature of `TC` degrees celcius and `P` ppmw phosphorous present,
    following the solubility model of Rusiecka & Baker, 2019
    (doi: 10.2138/am-2019-6931)
    """
    function Rusiecka_tmonaziteREE(P_ppm, TC)
        #[LREE][P] (mol^2/100g^2) = 10^(2.3055 - 1.029e4/T)
        Kₛₚ = 10^(2.3055 - 1.029e4/(TC + 273.15))
        LREE_μmol_g = Kₛₚ/(P_ppm/10000/30.97)*10000
        return LREE_μmol_g
    end

    """
    ```julia
    LREEt = Rusiecka_txenotimeY(P_ppm, TC)
    ```
    Calculate the Y (ppmw) concentration required for xenotime saturation at a
    temperature of `TC` degrees celcius and `P` ppmw phosphorous present,
    following the solubility model of Rusiecka & Baker, 2019
    (doi: 10.2138/am-2019-6931)
    """
    function Rusiecka_txenotimeY(P_ppm, TC)
        # [Y][P] (mol^2/100g^2) = 10^(3.6932 - 1.1469e4/T)
        Kₛₚ = 10^(3.6932 - 1.1469e4/(TC + 273.15))
        Y_ppm = Kₛₚ/(P_ppm/10000/30.97)*10000*88.905
        return Y_ppm
    end

## --- Apatite saturation calculations

    """
    ```julia
    P2O5 = Harrison_tapatiteP2O5(SiO2, Al2O3, CaO, Na2O, K2O, T)
    ```
    Calculate `P2O5` concentration (in wt.%) required for apatite saturation at a
    given `T` (in C) following the apatite saturation model of Harrison and Watson
    1984 (doi: 10.1016/0016-7037(84)90403-4) with the correction of Bea et al. 1992
    (doi: 10.1016/0024-4937(92)90033-U) where applicable
    """
    function Harrison_tapatiteP2O5(SiO2::T, Al2O3::T, CaO::T, Na2O::T, K2O::T, TC::T) where T <: Number
        TK = TC + 273.16
        ASI = (Al2O3/50.9806)/(CaO/56.0774 + Na2O/30.9895 + K2O/47.0827)
        P2O5sat = 52.5525567/exp( (8400 + 2.64e4(SiO2/100 - 0.5))/TK - (3.1 + 12.4(SiO2/100 - 0.5)) )
        return max(P2O5sat, P2O5sat * (ASI-1) * 6429/TK)
    end
    function Harrison_tapatiteP2O5(SiO2::T, TC::T) where T <: Number
        TK = TC + 273.16
        P2O5sat = 52.5525567/exp( (8400 + 2.64e4(SiO2/100 - 0.5))/TK - (3.1 + 12.4(SiO2/100 - 0.5)) )
        return P2O5sat
    end
    export Harrison_tapatiteP2O5

    """
    As `Harrison_tapatiteP2O5`, but returns saturation phosphorus concentration in PPM P
    """
    Harrison_tapatiteP(x...) = Harrison_tapatiteP2O5(x...) * 10_000 / 2.2913349
    export Harrison_tapatiteP

    """
    ```julia
    TC = Harrison_tapatite(SiO2, P2O5)
    ```
    Calculate apatite saturation temperature in degrees Celcius
    following the apatite saturation model of Harrison and Watson 1984
    (doi: 10.1016/0016-7037(84)90403-4)
    """
    function Harrison_tapatite(SiO2::T, P2O5::T) where T <: Number
        TK = (8400 + 2.64e4(SiO2/100 - 0.5)) / (log(52.5525567/P2O5) + (3.1 + 12.4(SiO2/100 - 0.5)))
        return TK - 273.16
    end
    export Harrison_tapatite

    """
    ```julia
    P2O5 = Tollari_tapatiteP2O5(SiO2, CaO, T)
    ```
    Calculate `P2O5` concentration (in wt.%) required for apatite saturation at a
    given `T` (in C) following the apatite saturation model of Tollari et al. 2006
    (doi: 10.1016/j.gca.2005.11.024)
    """
    function Tollari_tapatiteP2O5(SiO2::T, CaO::T, TC::T) where T <: Number
        # Using conversions from Tollari et al.
        SiO2ₘ = 1.11 * SiO2
        CaOₘ = 1.18 * CaO
        TK = TC+273.15
        P2O5satₘ = exp(TK * (-0.8579/(139.0 - SiO2ₘ) + 0.0165)  -  10/3*log(CaOₘ))
        return P2O5satₘ / 0.47
    end

    """
    ```julia
    TC = Tollari_tapatite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5)
    ```
    Calculate apatite saturation temperature in degrees Celcius
    following the apatite saturation model of Tollari et al. 2006
    (doi: 10.1016/j.gca.2005.11.024)
    """
    function Tollari_tapatite(SiO2, TiO2, Al2O3, FeOT, MgO, CaO, Na2O, K2O, P2O5)
        # Cations
        Na2 = Na2O/61.97854
        K2 = K2O/94.19562
        Ca = CaO/56.0774
        Al2 = Al2O3/101.960077
        Si = SiO2/60.0843
        Ti = TiO2/79.865
        Fe = FeOT/71.8444
        Mg = MgO/24.3050
        P2 = P2O5/141.94252

        # Normalize to mole percent oxides
        normconst = nansum((Na2, K2, Ca, Al2, Si, Ti, Fe, Mg, P2))
        CaOₘ = Ca / normconst * 100
        SiO2ₘ = Si / normconst * 100
        P2O5ₘ = P2 / normconst * 100
        TC = (log(P2O5ₘ) + 10/3*log(CaOₘ)) / (-0.8579/(139.0 - SiO2ₘ) + 0.0165) - 273.15
        return TC
    end

## --- Ti-in-zircon thermometry

    """
    ```julia
    Ti = Ferry_Ti_in_zircon(TC::Number, aSiO2::Number, aTiO2::Number)
    ```
    Parts per million by weight of titanium in zircon at temperature `TC` degrees
    Celsius given `aSiO2` silica activity and `aTiO2` titanium activity, following
    the equations of Ferry and Watson, 2007.
    (doi: 10.1007/s00410-007-0201-0)
    """
    function Ferry_Ti_in_zircon(TC::Number, aSiO2::Number, aTiO2::Number)
        exp10(5.711 - 4800.0/(TC+273.15) - log10(aSiO2) +log10(aTiO2))
    end

    """
    ```julia
    TC = Ferry_Ti_in_zirconT(TC::Number, aSiO2::Number, aTiO2::Number)
    ```
    Calculate titanium-in-zircon temperature in degrees Celcius `TC`
    given `Ti` parts per million by weight of titanium in zircon,
    `aSiO2` silica activity and `aTiO2` titanium activity, following
    the equations of Ferry and Watson, 2007.
    (doi: 10.1007/s00410-007-0201-0)
    """
    function Ferry_Ti_in_zirconT(Ti::Number, aSiO2::Number, aTiO2::Number)
        1 / ((5.711) - log10(aSiO2) + log10(aTiO2) - log10(Ti)) * (4800.0) .- 273.15
    end

    """
    ```julia
    Ti = Crisp_Ti_in_zircon(TC::Number, Pbar::Number, aSiO2::Number, aTiO2::Number)
    ```
    Parts per million by weight of titanium in zircon at temperature `TC` degrees
    Celsius and pressure `Pbar` bar given `aSiO2` silica activity and `aTiO2`
    titanium activity, following the equations of Crisp et al., 2023.
    (doi: 10.1016/j.gca.2023.04.031)
    """
    function Crisp_Ti_in_zircon(TC::Number, Pbar::Number, aSiO2::Number, aTiO2::Number)
        T = TC+273.15
        P = Pbar*1e-4
        f = 1.0/(1.0+10.0^(0.775P - 3.3713))
        exp10(5.84 - 4800.0/T - 0.12*P - 0.0056*P^3 - log10(aSiO2)*f +log10(aTiO2)) / f
    end

    """
    ```julia
    Zr = Ferry_Zr_in_rutile(TC::Number, aSiO2::Number)
    ```
    Parts per million by weight of zirconium in rutile at temperature `TC`
    degrees Celsius given `aSiO2` silica activity, following the
    equations of Ferry and Watson, 2007.
    (doi: 10.1007/s00410-007-0201-0)
    """
    function Ferry_Zr_in_rutile(TC::Number, aSiO2::Number)
        exp10(7.420 - 4530.0/(TC+273.15) - log10(aSiO2))
    end

    # calculate the temperature of rutile saturation in degrees Celsius
    """
    ```julia
    TC = Ferry_Zr_in_rutileT(Zr::Number, aSiO2::Number)
    ```
    Calculate zirconium-in-rutile temperature in degrees Celcius
    given `Zr` parts per million by weight of zirconium in rutile and 
    `aSiO2` silica activity, following the equations of Ferry and Watson, 2007.
    (doi: 10.1007/s00410-007-0201-0)
    """
    function Ferry_Zr_in_rutileT(Zr::Number, aSiO2::Number)
        1 / ((7.420) - log10(aSiO2) - log10(Zr)) * (4530.0) .- 273.15
    end


## --- End of File
