using StatGeochem
using Test, Statistics, StatsBase

@testset "ArrayStats" begin include("testArrayStats.jl") end
@testset "Import" begin include("testImport.jl") end
@testset "Other Utilities" begin include("testUtilities.jl") end
@testset "Crust 1.0" begin include("testCrust1.jl") end
@testset "Other Resources" begin include("testResources.jl") end
