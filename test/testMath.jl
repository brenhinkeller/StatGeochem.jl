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
    @test inpolygon([-1,0,1,0],[0,1,0,-1],prevfloat.([0.5,0.5]))
    @test !inpolygon([-1,0,1,0],[0,1,0,-1],nextfloat.([0.5,0.5]))
    @test inpolygon([-1,1,1,-1],[1,1,-1,-1],(0,0))
    @test !inpolygon([-1,1,1,-1],[1,1,-1,-1],(1.1,1))

    i,j = find_grid_inpolygon(-1.5:1/3:1.5, -1.5:1/3:1.5, [-.75,.75,.75,-.75],[.75,.75,-.75,-.75])
    @test sort([i j], dims=2) == [4 4; 4 5; 4 6; 4 7; 4 5; 5 5; 5 6; 5 7; 4 6; 5 6; 6 6; 6 7; 4 7; 5 7; 6 7; 7 7]

    @test all( arcdistance(0,100,[30,0,0],[100,100,95]) .≈ [30,0,5] )

## --- Weighted mean functions

    @test awmean(1:10,ones(10)) == (5.5, 0.31622776601683794, 9.166666666666666)
    @test gwmean(1:10,ones(10)) == (5.5, 0.9574271077563381, 9.166666666666666)


## --- Silly functions

    @test isapprox(StatGeochem.inv_sqrt(5.0), 1/sqrt(5.0), atol=1e-6)
    @test isapprox(StatGeochem.inv_sqrt(5f0), 1/sqrt(5f0), atol=1e-6)

## ---
