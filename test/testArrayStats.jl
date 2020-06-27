## --- ArrayStats.jl

    # Weighted Means
    @test MSWD([0,1,2],[1,1,1]) == 1.0
    @test awmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.5, 5/3)
    @test gwmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.6454972243679028, 5/3)

    # Summary statistics: simple cases, Float64
    A = [1:10.0..., NaN]
    @test StatGeochem.nanmin(1.,2.) == 1.
    @test StatGeochem.nanmax(1.,2.) == 2.
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

    # Arrays containing only NaNs should yield NaN
    A = fill(NaN,10)
    @test nansum(A) == 0
    @test isnan(nanmean(A))
    @test isnan(nanrange(A))
    @test isnan(nanminimum(A))
    @test isnan(nanmaximum(A))
    @test all(isnan.(nanextrema(A)))
    @test isnan(nanstd(A))
    @test isnan(nanaad(A))

    # Summary statistics: simple cases, Int64
    A = collect(1:10)
    @test StatGeochem.nanmin(1,2) == 1
    @test StatGeochem.nanmax(1,2) == 2
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3]) == 1.0
    @test nanmad([1,2,3]) == 1.0
    @test nanaad([1,2,3]) ≈ 2/3
    @test nanmedian([1,2,3]) == 2.0
    @test pctile([0:100...],99) == 99.0

    # Summary statistics: simple cases, ranges
    A = 1:10
    @test StatGeochem.nanmin(1,2) == 1
    @test StatGeochem.nanmax(1,2) == 2
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3) == 1.0
    @test nanmad(1:3) == 1.0
    @test nanaad(1:3) ≈ 2/3
    @test nanmedian(1:3) == 2.0
    @test pctile(0:100,99) == 99.0

    A = 1:10.
    @test StatGeochem.nanmin(1,2) == 1
    @test StatGeochem.nanmax(1,2) == 2
    @test nansum(A) == 55.0
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3.) == 1.0
    @test nanmad(1:3.) == 1.0
    @test nanaad(1:3.) ≈ 2/3
    @test nanmedian(1:3.) == 2.0
    @test pctile(0:100.,99) == 99.0

    # Summary statistics: dimensional tests, Int64
    A = reshape(1:300,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0,100,1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3,100,1)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)

    # Summary statistics: dimensional tests, Float64
    A = reshape(1:300.,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanmad(A, dims=1) == [25.0 25.0 25.0]
    @test nanmad(A, dims=2) == fill(100.0,100,1)
    @test nanaad(A, dims=1) == [25.0 25.0 25.0]
    @test nanaad(A, dims=2) ≈ fill(200/3,100,1)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)

    # Summary statistics: binning
    @test nanmean([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83]
    @test nanmean(1:100, reshape(1:300,100,3), 0, 100, 3) ==
                [17.0 117.0 217.0; 50.0 150.0 250.0; 83.0 183.0 283.0]
    @test nanmedian([1:100..., 1],[1:100..., NaN],0,100,3) == [17, 50, 83.5]
    @test nanmedian(1:100, reshape(1:300,100,3), 0, 100, 3) ==
                [17.0 117.0 217.0; 50.0 150.0 250.0; 83.5 183.5 283.5]

    # Standardization
    @test standardize!(collect(1:10.)) ≈ ((1:10) .- mean(1:10)) / std(1:10)

    # Sorting and counting
    A = rand(1:100.,100); B = sort(A)
    @test A[1:count_unique!(A)] == unique(B)

    # Interpolation
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]
    @test movmean(collect(1:10.),5) == movmean(1:10,5)
    @test movmean(1:10,4) == [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0]
    @test movmean(repeat(1:10,1,10),4) == repeat([2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0],1,10)

    # Integration
    @test trapz(1:10,fill(1,10)) == 9
    @test trapz(collect(1:10.),ones(10)) == 9
    @test midpointintegrate(1:10,ones(10)) == 10

## ---
