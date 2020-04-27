## --- ArrayStats.jl

    # Simple functions
    A = [1:10.0..., NaN]
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3,NaN]) == 1.0
    @test nanmad([1,2,3,NaN]) == 1.0
    @test nanaad([1,2,3,NaN]) ≈ 2/3
    @test nanmedian([1,2,3,NaN]) == 2.0
    @test pctile([0:100...,NaN],99) == 99.0

    # Binning
    @test nanmean([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83]
    @test nanmean(1:100, reshape(1:300,100,3), 0, 100, 3) ==
            [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]

    # Weighted Means
    @test MSWD([0,1,2],[1,1,1]) == 1.0
    @test awmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.5, 1.6666666666666667)
    @test gwmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.6454972243679028, 1.6666666666666667)


## --- Math.jl

    @test inpolygon([-1,0,1,0],[0,1,0,-1],[0,0])
    @test all( arcdistance(0,100,[30,0,0],[100,100,95]) .≈ [30,0,5] )
