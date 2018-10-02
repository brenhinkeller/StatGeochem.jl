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
    end

    include("utilities/System.jl");
    include("utilities/Math.jl");
    include("utilities/Import.jl");

    using StatsBase: percentile
    using Interpolations: interpolate, Gridded, Linear
    include("utilities/ArrayStats.jl");
    include("utilities/Resampling.jl");

    using Colors: ColorTypes, RGB, N0f8
    include("utilities/Colormaps.jl");

    include("utilities/Geochemistry.jl");
    include("utilities/GIS.jl");
    include("utilities/Etc.jl");

    using HDF5
    include("Resources.jl")

end # module
