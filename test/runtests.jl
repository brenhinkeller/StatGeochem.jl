using StatGeochem
using Test, Statistics, StatsBase

# Utilities
@testset "Math" begin include("testMath.jl") end
@testset "ArrayStats" begin include("testArrayStats.jl") end
@testset "Import" begin include("testImport.jl") end
@testset "Changepoint" begin include("testChangepoint.jl") end
@testset "Other Utilities" begin include("testUtilities.jl") end

# Resources
@testset "Crust 1.0" begin include("testCrust1.jl") end
@testset "Other Resources" begin include("testResources.jl") end
