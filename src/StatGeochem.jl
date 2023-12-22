__precompile__()

module StatGeochem

    using Reexport
    @reexport using NaNStatistics
    @reexport using StatGeochemBase

    # Vectorization and parallelization tools
    using LoopVectorization: @turbo
    using Polyester: @batch

    # General requirements
    using DelimitedFiles, Random, Downloads
    using ProgressMeter: @showprogress, Progress, update!, next!
    const Collection{T} = Union{AbstractArray{<:T}, NTuple{N,T}} where N
    include("utilities/System.jl")
    include("utilities/Import.jl")
    include("utilities/Resampling.jl")
    include("utilities/Changepoint.jl")

    include("resources/Chemistry.jl")
    include("utilities/Geochronology.jl")
    include("utilities/Geochemistry.jl")
    include("utilities/GIS.jl")
    include("utilities/Etc.jl")

    # Resources
    using LazyArtifacts
    resourcepath = joinpath(homedir(),"resources")
    moduleresourcepath = joinpath(Base.source_dir(),"resources")
    export resourcepath, moduleresourcepath

    using FileIO: load
    using HDF5: h5read
    using Colors: Color, RGBX, RGB, N0f8
    include("resources/tc1/tc1.jl")
    include("resources/Crust1.jl")
    include("resources/Litho1.jl")
    include("resources/Geography.jl")
    include("resources/Elevation.jl")
    include("resources/Seafloorage.jl")
    include("resources/PartitionCoefficients/PartitionCoefficients.jl")

    # Custom pretty printing for some types
    include("utilities/Display.jl")

    # Methods exclusively for package extensions
    mapplot() = @warn "`mapplot` requires Plots.jl"
    mapplot!() = @warn "`mapplot!` requires Plots.jl"
    spidergram() = @warn "`spidergram` requires Plots.jl"
    spidergram!() = @warn "`spidergram!` requires Plots.jl"
    export mapplot, mapplot!, spidergram, spidergram!

end # module
