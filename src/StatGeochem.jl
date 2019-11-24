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

    include("utilities/System.jl");
    include("utilities/Math.jl");
    include("utilities/Import.jl");

    using StatsBase: percentile
    using ProgressMeter: @showprogress, Progress, update!
    using Interpolations
    include("utilities/ArrayStats.jl");
    include("utilities/Resampling.jl");

    using IndirectArrays: IndirectArray
    using Colors: Colorant, ColorTypes, RGB4, RGB, N0f8
    include("utilities/Colormaps.jl");

    include("utilities/Geochemistry.jl");
    include("utilities/GIS.jl");
    include("utilities/Etc.jl");

    using FileIO, HDF5
    tc1_550 = readdlm("resources/tc1_550.csv", ',')
    tc1_1300 = readdlm("resources/tc1_1300.csv", ',')
    tc1_age = readdlm("resources/tc1_age.csv", ',', Int)
    include("Resources.jl")

end # module
