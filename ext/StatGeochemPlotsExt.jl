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


    taylormclennan = (
        La = 0.367,Ce = 0.957,Pr = 0.137,Nd = 0.711,Sm = 0.231,Eu = 0.087,Gd = 0.306,
        Tb = 0.058,Dy = 0.381,Ho = 0.085,Er = 0.249,Tm = 0.036,Yb = 0.248,Lu = 0.038,
    )

    """
    Construct a `chondrite` normalized multi-element diagram (spider diagram) from the 
    rare earth elements in `data`. 

    Use `spidergram` to create a new plot object, and `spidergram!` to add to an existing 
    one:
    ```julia
    spidergram(data; [chondrite], kwargs...)              # Create a new spider diagram
    spidergram!(plotobj, data; [chondrite], kwargs...)    # Add to the plot `plotobj`
    ```

    If no chondrite values are specified, `data` will be normalized to the values reported
    by Taylor and McLennan (1985).

    Values in `data` and `chondrite` may be passed as a dictonary, named tuple, or an 
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

    If `data` is not passed as an array, `chondrite` must be a named tuple in element 
    order.
    """
    StatGeochem.spidergram(data; chondrite=taylormclennan, kwargs...) = 
        spidergram!(plot(), data; chondrite=chondrite, kwargs...)
    
    function StatGeochem.spidergram!(h, data::Dict; chondrite::NamedTuple=taylormclennan, kwargs...)
        Key = keytype(data)
        data = [data[Key(k)] for k in keys(chondrite)]
        _spidergram!(h, data, values(chondrite); kwargs...)
    end
    function StatGeochem.spidergram!(h, data::NamedTuple; chondrite::NamedTuple=taylormclennan, kwargs...)
        data = [data[Symbol(k)] for k in keys(chondrite)]
        _spidergram!(h, data, values(chondrite); kwargs...)
    end

    StatGeochem.spidergram!(h, data::AbstractArray; chondrite=taylormclennan, kwargs...) = 
        _spidergram!(h, data, collect(values(chondrite)); kwargs...)

    function _spidergram!(h, data::AbstractArray, chondrite; kwargs...,)
        Plots.plot!(h, 
            ylabel="Chondrite Normalized",
            fg_color_legend=:white,
            framestyle=:box,
            grid=false,
            yaxis=:log10,
            ylims=(10^0, 10^3),
            yticks=(10.0.^(0:3), ("1", "10", "100", "1000")),
            xticks=(1:15, ["La","Ce","Pr","Nd","","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
                "Yb","Lu"]),
            yminorticks=log.(1:10),
        )
        Plots.plot!(h, collect([1:4; 6:15]), data ./ chondrite; kwargs...)

        return h
    end

    export spidergram, spidergram!
end
