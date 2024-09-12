using StatGeochem
using Test, Statistics, Downloads

@testset "All tests" begin
    # Utilities
    @testset "Resampling" begin include("testResampling.jl") end
    @testset "Changepoint" begin include("testChangepoint.jl") end
    @testset "Geochemistry" begin include("testGeochemistry.jl") end

    using FileIO
    @testset "Other Utilities" begin include("testUtilities.jl") end

    # Resources
    @testset "Perplex" begin include("testPerplex.jl") end
    @testset "Crust 1.0" begin include("testCrust1.jl") end
    @testset "Litho 1.0" begin include("testLitho1.jl") end
    @testset "Other Resources" begin include("testResources.jl") end

    using Plots
    @testset "Package Extensions" begin include("testExtensions.jl") end
end
