__precompile__()

module StatGeochem

    using Reexport
    @reexport using NaNStatistics
    @reexport using StatGeochemBase

    # AVX vectorziation tools
    using LoopVectorization

    # General requirements
    using Statistics, DelimitedFiles, SpecialFunctions, Random
    using StatsBase: percentile, mean, std, ProbabilityWeights
    using ProgressMeter: @showprogress, Progress, update!
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
    resourcepath = joinpath(homedir(),"resources")
    moduleresourcepath = joinpath(Base.source_dir(),"resources")
    export resourcepath, moduleresourcepath

    using FileIO, ImageIO, HDF5
    using Colors: Color, RGBX, RGB, N0f8
    include("resources/tc1/tc1.jl")
    include("resources/Crust1.jl")
    include("resources/Elevation.jl")
    include("resources/Seafloorage.jl")
    include("resources/PartitionCoefficients/PartitionCoefficients.jl")

    # Custom pretty printing for some types
    include("utilities/Display.jl")

end # module
