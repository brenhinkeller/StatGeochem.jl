## --- Distribution functions

    @test normpdf.(0, 1,-1:1) ≈ [0.24197072451914337, 0.3989422804014327, 0.24197072451914337]
    @test normpdf_ll.(0,1,-5:5) == -(-5:5).^2/2
    @test normpdf_ll(0,1,-5:5) ≈ sum(normpdf_ll.(0,1,-5:5))
    @test normcdf(0,1,-2:2) ≈ normcdf.(0,1,-2:2) ≈ [0.02275013194817921, 0.15865525393145707, 0.5, 0.8413447460685429, 0.9772498680518208]

    @test normproduct(0,1,0,1) === normpdf(0,sqrt(2),0) === 0.28209479177387814
    @test normproduct_ll(0,1,0,1) === normpdf_ll(0,1,0) === 0.0

    @test [-2,0,2] ≈ norm_quantile.([0.022750131948, 0.5, 0.977249868052])
    @test norm_quantile.(0:0.25:1) ≈ [-Inf, -0.6744897501960818, 0.0, 0.6744897501960818, Inf]

    @test isapprox(norm_width(390682215445)/2, 7, atol=1e-5)


## -- Geometry functions

    @test inpolygon([-1,0,1,0],[0,1,0,-1],[0,0])
    @test !inpolygon([-1,0,1,0],[0,1,0,-1],[0,10])
    @test inpolygon([-1,1,1,-1],[1,1,-1,-1],(0,0))
    @test !inpolygon([-1,1,1,-1],[1,1,-1,-1],(1.1,1))

    @test all( arcdistance(0,100,[30,0,0],[100,100,95]) .≈ [30,0,5] )

## --- Weighted mean functions

    @test awmean(1:10,ones(10)) == (5.5, 0.31622776601683794, 9.166666666666666)
    @test gwmean(1:10,ones(10)) == (5.5, 0.9574271077563381, 9.166666666666666)

## ---
