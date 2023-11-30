module StatGeochemPlotsExt

    using StatGeochem, Plots

    mapplot(args...; kwargs...) = mapplot!(plot(), args...; kwargs...)
    function mapplot!(h, args...; file="world2048.jpg", seriestype=:scatter, kwargs...)
        if !isfile(file)
            file = joinpath(moduleresourcepath, "maps", file)
        end
        img = reverse!(load(file), dims=1)
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
    end

    export mapplot, mapplot!
end
