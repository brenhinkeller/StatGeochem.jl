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
end
