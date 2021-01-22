# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                StatGeochem.jl                                 #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#   Last modified by C. Brenhin Keller 2018-07-12                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

__precompile__()

module StatGeochem

    # Backwards compatibility
    using Compat
    # Forwards compatibility
    VERSION >= v"0.7" && using Statistics, DelimitedFiles, SpecialFunctions

    # AVX vectorziation tools
    using LoopVectorization
    using SIMDPirates: vifelse, verf
    using VectorizationBase: SVec

    import SpecialFunctions.erf
    erf(x::SVec) = verf(x)

    # General requirements
    using Random
    using StatsBase: percentile, mean, std, ProbabilityWeights
    using ProgressMeter: @showprogress, Progress, update!
    using Interpolations
    include("utilities/System.jl")
    include("utilities/Math.jl")
    include("utilities/Import.jl")
    include("utilities/ArrayStats.jl")
    include("utilities/Resampling.jl")
    include("utilities/Changepoint.jl")

    using IndirectArrays: IndirectArray
    using Colors: Color, RGBX, RGB, N0f8
    include("resources/Colormaps.jl")

    include("resources/Chemistry.jl")
    include("utilities/Geochronology.jl")
    include("utilities/Geochemistry.jl")
    include("utilities/GIS.jl")
    include("utilities/Etc.jl")

    # Resources
    resourcepath = joinpath(homedir(),"resources")
    moduleresourcepath = joinpath(Base.source_dir(),"resources")
    export resourcepath, moduleresourcepath

    using FileIO, HDF5
    include("resources/tc1/tc1.jl")
    include("resources/Crust1.jl")
    include("resources/Elevation.jl")
    include("resources/Seafloorage.jl")
    include("resources/PartitionCoefficients/PartitionCoefficients.jl")

    # Custom pretty printing for some types
    include("utilities/Display.jl")

end # module
