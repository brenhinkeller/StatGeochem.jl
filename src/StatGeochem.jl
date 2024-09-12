__precompile__()

module StatGeochem

    using Reexport
    @reexport using NaNStatistics
    @reexport using StatGeochemBase

    # Vectorization and parallelization tools
    using LoopVectorization: @turbo
    using Polyester: @batch

    # General requirements
    using DelimitedFiles, Random, Downloads, LazyArtifacts
    using ProgressMeter: @showprogress, Progress, update!, next!
    const Collection{T} = Union{DenseArray{<:T}, AbstractRange{<:T}, NTuple{N,T}} where N
    include("utilities/System.jl")
    include("utilities/Resampling.jl")
    include("utilities/Changepoint.jl")

    include("resources/Chemistry.jl")
    include("utilities/Geochronology.jl")
    include("utilities/Geochemistry.jl")
    include("utilities/GIS.jl")
    include("utilities/Etc.jl")

    # Resources
    resourcepath = joinpath(homedir(),"resources")
    moduleresourcepath = joinpath(Base.source_dir(),"resources")
    export resourcepath, moduleresourcepath

    using Perple_X_jll
    include("resources/Perplex.jl")

    using FileIO: load
    using HDF5: h5read
    using Colors: Color, RGBX, RGB, N0f8
    include("resources/tc1/tc1.jl")
    include("resources/Crust1.jl")
    include("resources/Litho1.jl")
    include("resources/Geology.jl")
    include("resources/Geography.jl")
    include("resources/Seafloorage.jl")
    include("resources/PartitionCoefficients/PartitionCoefficients.jl")

    # Custom pretty printing for some types
    include("utilities/Display.jl")

    # Functions for which methods will be added in package extensions
    function mapplot end
    function mapplot! end
    function spidergram end
    function spidergram! end
    export mapplot, mapplot!, spidergram, spidergram!

end # module
