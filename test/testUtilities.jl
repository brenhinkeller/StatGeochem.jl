## --- Geochronology.jl

    @test awmean(1:10,ones(10)) == (5.5, 0.31622776601683794, 9.166666666666666)
    @test gwmean(1:10,ones(10)) == (5.5, 0.9574271077563381, 9.166666666666666)

## --- Math.jl

    @test normpdf.(0,1,-1:1) ≈ [0.24197072451914337, 0.3989422804014327, 0.24197072451914337]
    @test normpdf_ll.(0,1,-5:5) == -(-5:5).^2/2
    @test normpdf_ll(0,1,-5:5) == sum(normpdf_ll.(0,1,-5:5))
    @test normcdf(0,1,-2:2) == normcdf.(0,1,-2:2) == [0.02275013194817921, 0.15865525393145707, 0.5, 0.8413447460685429, 0.9772498680518208]

    @test norm_quantile.(0:0.25:1) == [-Inf, -0.6744897501960818, 0.0, 0.6744897501960818, Inf]
    @test norm_width.(0:0.25:1) ≈ [-Inf, -1.6832424671458286, -0.8614545985909146, -0.3600247395854099, 0.0]
    @test normproduct(0,1,0,1) === normpdf(0,sqrt(2),0) === 0.28209479177387814
    @test normproduct_ll(0,1,0,1) === normpdf_ll(0,1,0) === 0.0

    @test inpolygon([-1,0,1,0],[0,1,0,-1],[0,0])
    @test all( arcdistance(0,100,[30,0,0],[100,100,95]) .≈ [30,0,5] )

## --- Resampling.jl

    @test bsr(1:10,fill(0.5,10),1000,0.5)::Array{Float64} |> length == 1000
    @test bsr(1:10,fill(0.5,10),1000,fill(0.5,10))::Array{Float64} |> length == 1000
    @test bsr!(Array{Int64}(undef,1000),1:10,fill(0,10),1000,0.5) |> unique ⊆ 1:10
    @test bsr!(Array{Int64}(undef,1000),1:10,fill(0,10),1000,fill(0.5,1000)) |> unique ⊆ 1:10

    @test bincounts(1:100,0,100,10) == (5:10:95, fill(10,10))
    @test binmeans(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    @test binmedians(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(1.1720982147414096,10))

    @test randsample(1:10,1000)::Array{Int64} |> length == 1000
    @test unique(randsample(1:10,1000)) ⊆ 1:10

## ---
