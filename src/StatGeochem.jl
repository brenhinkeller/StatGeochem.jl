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

    include("utilities/Math.jl");

    using StatsBase: percentile
    using Interpolations: interpolate, Gridded, Linear
    include("utilities/ArrayStats.jl");

    using Images: ColorTypes, N0f8
    include("utilities/Colormaps.jl");

    include("utilities/Etc.jl");

    include("utilities/Resources.jl");

end # module
