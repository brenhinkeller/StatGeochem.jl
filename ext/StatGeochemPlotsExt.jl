module StatGeochemPlotsExt

    using StatGeochem, Plots

    StatGeochem.mapplot(args...; kwargs...) = mapplot!(plot(), args...; kwargs...)
    function StatGeochem.mapplot!(h, args...; file="world2048.jpg", seriestype=:scatter, kwargs...)
        if !isfile(file)
            file = joinpath(StatGeochem.moduleresourcepath, "maps", file)
        end
        img = reverse!(StatGeochem.load(file), dims=1)
        x = range(-180, 180, length=size(img,2))
        y = range(-90, 90, length=size(img,1))
        plot!(h, x, y, img,
            framestyle=:box,
            yflip=false,
            xlabel="Longitude",
            ylabel="Latitude",
            xlims=(-180, 180),
            ylims=(-90, 90),
            size=(800,400),
        )
        plot!(h, args...; seriestype, kwargs...)
        return h
    end

    export mapplot, mapplot!


    chondrite = (
        La = 0.367,Ce = 0.957,Pr = 0.137,Nd = 0.711,Sm = 0.231,Eu = 0.087,Gd = 0.306,
        Tb = 0.058,Dy = 0.381,Ho = 0.085,Er = 0.249,Tm = 0.036,Yb = 0.248,Lu = 0.038,
    )
    PAAS = (
        La = 38,Ce = 80,Pr = 8.9,Nd = 32,Sm = 5.6,Eu = 1.1,Gd = 4.7,
        Tb = 0.77,Dy = 4.4,Ho = 1.0,Er = 2.9,Tm = 0.40,Yb = 2.8,Lu = 0.43,
    )

    """
    Construct a normalized multi-element diagram (spider diagram) from the rare earth 
    elements in `data`. 

    Use `spidergram` to create a new plot object, and `spidergram!` to add to an existing 
    one:
    ```julia
    spidergram(data; [normalizer], kwargs...)              # Create a new spider diagram
    spidergram!(plotobj, data; [normalizer], kwargs...)    # Add to the plot `plotobj`
    ```

    Specify a `normalizer` set of REE values (in ppm) to normalize `data.` Chondrite and 
    PAAS (post-Archean Australian shale) from Taylor and McLennan (1985) are pre-defined 
    respectively as `chondrite` and `PAAS`. If no values are specified, `data` will be 
    normalized to chondrite values.

    Values in `data` and `normalizer` may be passed as a dictonary, named tuple, or an 
    array. All arrays should be in element order:

        La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu

    Dictonaries and NamedTuples should be organized by element:

        NamedTuple with 14 elements:
        La  = Float64 0.367
        Ce  = Float64 0.957
        Pr  = Float64 0.137
        ⋮   = ⋮

    Or

        Dict{String, Float64} with 14 entries:
        "La" => 20.78
        "Ce" => 35.61
        "Pr" => 2.344
        ⋮    => ⋮

    If `data` is not passed as an array, `normalizer` must be a named tuple in element 
    order.

    """
    function StatGeochem.spidergram(data; normalizer=chondrite, markershape=:circle, kwargs...)
        h = Plots.plot(
            fg_color_legend=:white,
            framestyle=:box,
            grid=false,
            yaxis=:log10,
            xticks=(1:15, ["La","Ce","Pr","Nd","","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                "Yb","Lu"]),
            yminorticks=log.(1:10),
        )

        if normalizer==chondrite
            ylabel!("Chondrite Normalized")
            ylims!(10^0, 10^3)
            yticks!(10.0.^(0:3), ["1", "10", "100", "1000"])
        elseif normalizer==PAAS 
            ylabel!("PAAS Normalized")
            ylims!(10^-1, 10^2)
            yticks!(10.0.^(-1:2), ["0.1", "1", "10", "100",])
        end

        spidergram!(h, data; normalizer=normalizer, markershape=markershape, kwargs...)
    end
    
    function StatGeochem.spidergram!(h, data::Dict; normalizer::NamedTuple=chondrite, 
            markershape=:circle, kwargs...
        )

        REEindex = NamedTuple{keys(normalizer)}(i for i in collect([1:4; 6:15]))
        Key = keytype(data)
        
        x = collect(values(REEindex))
        y = [(haskey(data, Key(k)) ? data[Key(k)]/normalizer[Symbol(k)] : NaN) for k in keys(normalizer)]
        _spidergram!(h, x, y; markershape=markershape, kwargs...)
    end
    
    function StatGeochem.spidergram!(h, data::NamedTuple; normalizer::NamedTuple=chondrite, 
            markershape=:circle, kwargs...
        )
        
        REEindex = NamedTuple{keys(normalizer)}(i for i in collect([1:4; 6:15]))
    
        x = [REEindex[Symbol(k)] for k in keys(data)]
        y = [data[k]/normalizer[k] for k in keys(data)]
        _spidergram!(h, x, y; markershape=markershape, kwargs...)
    end
    
    StatGeochem.spidergram!(h, data::AbstractArray; normalizer=chondrite, 
            markershape=:circle, kwargs...) = 
        _spidergram!(h, collect([1:4; 6:15]), data ./ collect(values(normalizer)); 
            markershape=markershape, kwargs...
        )
    
    function _spidergram!(h, x::AbstractArray, y::AbstractArray; kwargs...,)
        Plots.plot!(h, x[.!isnan.(y)], y[.!isnan.(y)]; kwargs...)
        return h
    end
    
    export spidergram, spidergram!
end
