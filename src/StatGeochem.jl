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
    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    else
        # Other compatibility not covered by Compat.jl:
        # New syntax for trunc
        import Base.trunc
        trunc(x; digits::Int=0) = trunc(x,digits)
    end

    # Functions
    include("utilities/System.jl");
    include("utilities/Math.jl");
    include("utilities/Import.jl");

    using StatsBase: percentile, mean, std, ProbabilityWeights
    using ProgressMeter: @showprogress, Progress, update!
    using Interpolations
    include("utilities/ArrayStats.jl");
    include("utilities/Resampling.jl");

    using IndirectArrays: IndirectArray
    using Colors: Colorant, ColorTypes, RGBX, RGB, N0f8
    include("utilities/Colormaps.jl");

    include("utilities/Geochemistry.jl");
    include("utilities/GIS.jl");
    include("utilities/Etc.jl");

    # Resources
    resourcepath = joinpath(homedir(),"resources")
    moduleresourcepath = joinpath(Base.source_dir(),"resources")
    export resourcepath, moduleresourcepath

    using FileIO, HDF5
    include("resources/tc1/tc1.jl")
    include("resources/Crust1.jl")
    include("resources/Elevation.jl")
    include("resources/Seafloorage.jl")


end # module
