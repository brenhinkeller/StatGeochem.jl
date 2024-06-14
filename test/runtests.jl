using StatGeochem
using Test, Statistics, StatsBase, Downloads

@testset "All tests" begin
    # Utilities
    @testset "Import" begin include("testImport.jl") end
    @testset "Resampling" begin include("testResampling.jl") end
    @testset "Changepoint" begin include("testChangepoint.jl") end
    @testset "Geochemistry" begin include("testGeochemistry.jl") end

    using FileIO
    @testset "Other Utilities" begin include("testUtilities.jl") end

    # Resources
    @testset "Crust 1.0" begin include("testCrust1.jl") end
    @testset "Litho 1.0" begin include("testLitho1.jl") end
    @testset "Other Resources" begin include("testResources.jl") end

    using Plots
    @testset "Package Extensions" begin include("testExtensions.jl") end
end
