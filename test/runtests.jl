using StatGeochem
using Test, Statistics, StatsBase

# Utilities
@testset "Import" begin include("testImport.jl") end
@testset "Resampling" begin include("testResampling.jl") end
@testset "Changepoint" begin include("testChangepoint.jl") end
@testset "Other Utilities" begin include("testUtilities.jl") end

# Resources
@testset "Crust 1.0" begin include("testCrust1.jl") end
@testset "Other Resources" begin include("testResources.jl") end
