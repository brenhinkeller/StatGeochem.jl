using StatGeochem
using Test, Statistics, StatsBase, Downloads, Perple_X_jll

#TODO: move back to Resources when done testing
@testset "Other Resources" begin include("testResources.jl") end
@testset "Perplex" begin include("testPerplex.jl") end

# Utilities
@testset "Import" begin include("testImport.jl") end
@testset "Resampling" begin include("testResampling.jl") end
@testset "Changepoint" begin include("testChangepoint.jl") end
@testset "Geochemistry" begin include("testGeochemistry.jl") end

using ImageIO, FileIO
@testset "Other Utilities" begin include("testUtilities.jl") end

# Resources
@testset "Crust 1.0" begin include("testCrust1.jl") end
@testset "Litho 1.0" begin include("testLitho1.jl") end
