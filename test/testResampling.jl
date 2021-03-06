## --- Resampling.jl

    index = Array{Int64}(undef,1000)
    resampled = Array{Int64}(undef,1000)
    @test bsr!(resampled,index,1:10,0,0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,0,fill(0.5,1000)) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),fill(0.5,1000)) |> unique ⊆ 1:10

    resampled = Array{Float64}(undef,1000)
    # Gaussian
    bsr!(resampled, index, ones(10), 0.5, 0.5)
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(resampled, index, ones(10), 0.5, fill(0.5,1000))
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(resampled, index, ones(10), fill(0.5,10), 0.5)
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(resampled, index, ones(10), fill(0.5,10), fill(0.5,1000))
    @test isapprox(mean(resampled), 1, atol=0.3)
    # Other distributions
    bsr!(uniform, resampled, index, ones(10), 0.5, 0.5)
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(triangular, resampled, index, ones(10), 0.5, fill(0.5,1000))
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(triangular, resampled, index, ones(10), fill(0.5,10), 0.5)
    @test isapprox(mean(resampled), 1, atol=0.3)
    bsr!(uniform, resampled, index, ones(10), fill(0.5,10), fill(0.5,1000))
    @test isapprox(mean(resampled), 1, atol=0.3)

    @test bsresample(1:10,fill(0.5,10),1000,0.5)::Array{Float64} |> length == 1000
    @test bsresample(1:10,fill(0.5,10),1000,fill(0.5,10))::Array{Float64} |> length == 1000

    d = Dict{String,Vector{Float64}}()
    d["x"] = 1:10;  d["y"] = 2:2:20
    d["x_sigma"] = d["y_sigma"] = fill(0.5,10)
    dbs = bsresample(d, 1000, ["x","y"], 0.5)
    @test isapprox(mean(dbs["x"]), 5.5, atol=0.5)
    @test isapprox(std(dbs["x"]), 3.03, atol=0.5)
    @test isapprox(mean(dbs["y"]), 11, atol=1)
    @test isapprox(std(dbs["y"]), 6.06, atol=1)

    @test bincounts(1:100, 0, 100, 10) == (5:10:95, fill(10,10))
    @test binmeans(1:100, 1:100, 0, 100, 10) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    @test binmeans(1:100, 1:100, 0, 100, 10, ones(100)) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    @test binmedians(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(1.1720982147414096,10))

    @test randsample(1:10., 1000, rand(1000))::Array{Float64} |> length == 1000
    @test randsample(1:10, 1000)::Array{Int64} |> length == 1000
    @test unique(randsample(1:10, 1000)) ⊆ 1:10

    dr = randsample(d, 1000, ["x","y"], 0.5)
    @test isapprox(mean(dr["x"]), 5.5, atol=0.5)
    @test isapprox(std(dr["x"]), 3.03, atol=0.5)
    @test unique(dr["x"]) ⊆ 1:10
    @test isapprox(mean(dr["y"]), 11, atol=1)
    @test isapprox(std(dr["y"]), 6.06, atol=1)
    @test unique(dr["y"]) ⊆ 2:2:20

## --- Invweight

    @test invweight(0:10, 0:10, 0:10) ≈ [13.092772378121769, 13.759663290331229, 14.079390874013654, 14.244556812410089, 14.327747696132253, 14.354508911206949, 14.331218676773712, 14.251150311763046, 14.088257739618454, 13.76917452212827, 13.101581593462868]
    @test invweight_location(0:10, 0:10) ≈ [2.3478642777118957, 2.950680392065609, 3.2200556525889006, 3.348975353894235, 3.4103053158720016, 3.4297605864726126, 3.413776296513463, 3.3555688532471923, 3.2289225181937002, 2.9601916238626513, 2.3566734930529947]
    @test invweight_age(0:10) ≈ [10.744908100409873, 10.808982898265619, 10.859335221424754, 10.895581458515855, 10.91744238026025, 10.924748324734335, 10.91744238026025, 10.895581458515855, 10.859335221424754, 10.808982898265619, 10.744908100409873]

## --- bin_bsr

    x = 0:100; y = 0:100
    xmin = 0; xmax = 100; nbins = 5
    (c,m,e) = bin_bsr(x, y, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.4)
    @test isapprox(e, [1.17, 1.21, 1.23, 1.26, 1.28], atol=0.4)

    # Upper and lower CIs
    (c,m,el,eu) = bin_bsr(nanbinmean!, x, y, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.4)
    @test isapprox(el, [2.29, 2.38, 2.41, 2.49, 2.51], atol=0.5)
    @test isapprox(eu, [2.3, 2.37, 2.42, 2.51, 2.51], atol=0.5)

    # Medians, upper and lower CIs
    (c,m,el,eu) = bin_bsr(nanbinmedian!, x, y, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.01, 29.91, 49.9, 69.88, 89.79], atol=1)
    @test isapprox(el, [4.01, 3.91, 3.9, 3.88, 3.79], atol=2)
    @test isapprox(eu, [3.99, 4.09, 4.1, 4.12, 4.21], atol=2)

    # with weights
    w = ones(101)
    (c,m,e) = bin_bsr(x, y, xmin, xmax, nbins, w, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.4)
    @test isapprox(e, [1.17, 1.21, 1.23, 1.26, 1.28], atol=0.4)

    # with 2-D array (matrix) of y data
    y = repeat(0:100, 1, 4)
    y_sigma = ones(101,4)
    (c,m,e) = bin_bsr(x, y, xmin, xmax, nbins, x_sigma=ones(101), y_sigma=y_sigma)
    @test c == 10.0:20.0:90.0
    @test isapprox(m, repeat([10.04, 29.94, 49.94, 69.92, 89.83], 1, 4), atol=0.5)
    @test isapprox(e, repeat([1.17, 1.21, 1.23, 1.26, 1.28], 1, 4), atol=0.5)


## -- bin_bsr_ratios

    x = 0:100; num = 0:100; denom=reverse(num)
    xmin = 0; xmax = 100; nbins = 5
    (c,m,el,eu) = bin_bsr_ratios(x, num, denom, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.33, 8.99], atol=0.4)
    @test isapprox(el, [0.03, 0.05, 0.09, 0.26, 2.11], atol=0.5)
    @test isapprox(eu, [0.03, 0.05, 0.1, 0.29, 3.03], atol=0.5)

    # With weights
    (c,m,el,eu) = bin_bsr_ratios(x, num, denom, xmin, xmax, nbins, ones(101), x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.33, 8.99], atol=0.4)
    @test isapprox(el, [0.03, 0.05, 0.09, 0.26, 2.11], atol=0.5)
    @test isapprox(eu, [0.03, 0.05, 0.1, 0.29, 3.03], atol=0.5)

## --- Monte Carlo interpolation/fitting

    (c,m) = mcfit(0:11, ones(12), 0:11, ones(12), 1, 10, 10)
    @test c == 1:10
    @test isapprox(m, [1.15, 2.02, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.98, 9.85], atol = 0.25)

## --- Downsampling

    @test downsample(1:100,10) == 10:10:100
    @test downsample(reshape(1:100,10,10),2) == [12 32 52 72 92; 14 34 54 74 94; 16 36 56 76 96; 18 38 58 78 98; 20 40 60 80 100]

## ---
