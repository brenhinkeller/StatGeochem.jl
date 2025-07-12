## --- Resampling.jl

    index = Array{Int64}(undef,1000)
    resampled = Array{Int64}(undef,1000)
    @test bsr!(resampled,index,1:10,0,0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,0,fill(0.5,1000)) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),fill(0.5,1000)) |> unique ⊆ 1:10

    index = Array{Int64}(undef,1000)
    resampled = Array{Int64}(undef,1000,10)
    data = repeat(1:10,1,10)
    @test vec(bsr!(resampled,index,data,zeros(10),0.5)) |> unique ⊆ 1:10
    @test vec(bsr!(resampled,index,data,zeros(10),fill(0.5,10))) |> unique ⊆ 1:10

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

    dt = TupleDataset(d)
    @test isa(dt, NamedTuple)
    dbs = bsresample(dt, 1000, (:x, :y), 0.5)
    @test isapprox(mean(dbs[:x]), 5.5, atol=0.5)
    @test isapprox(std(dbs[:x]), 3.03, atol=0.5)
    @test isapprox(mean(dbs[:y]), 11, atol=1)
    @test isapprox(std(dbs[:y]), 6.06, atol=1)

    @test bincounts(1:100, 0.0:10:100.0) == (5:10:95, fill(10,10))
    @test bincounts(1:100, 0, 100, 10) == (5:10:95, fill(10,10))
    c,N = bincounts(1:100, 0, 100, 10; relbinwidth=3)
    @test c == 5:10:95
    @test N == [20, 30, 30, 30, 30, 30, 30, 30, 30, 20]

    @test binmeans(1:100, 1:100, 0.0:10:100.0) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    @test binmeans(1:100, 1:100, 0, 100, 10) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    c,m,e = binmeans(1:100, 1:100, 0, 100, 10; relbinwidth=3)
    @test c == 5:10:95
    @test m == [10.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 90.5]
    @test e ≈ [1.3228756555322954, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.3228756555322954]
    @test binmeans(1:100, 1:100, 0, 100, 10, ones(100)) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    c,m,e = binmeans(1:100, 1:100, 0, 100, 10, ones(100); relbinwidth=3)
    @test c == 5:10:95
    @test m == [10.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 90.5]
    @test e ≈ [1.3228756555322954, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.607275126832159, 1.3228756555322954]

    @test binmedians(1:100,1:100,0.0:10:100.0) == (5:10:95, 5.5:10:95.5, fill(1.1720982147414096,10))
    @test binmedians(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(1.1720982147414096,10))
    c,m,e = binmedians(1:100, 1:100, 0, 100, 10; relbinwidth=3)
    @test c == 5:10:95
    @test m == [10.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 90.5]
    @test e ≈ [1.6575971917205938, 2.030133659392898, 2.030133659392898, 2.030133659392898, 2.030133659392898, 2.030133659392898, 2.030133659392898, 2.030133659392898, 2.030133659392898, 1.6575971917205938]

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

    @test invweight(0:10,0:10,0:10) ≈ [13.092772378121769, 13.759663290331229, 14.079390874013654, 14.244556812410089, 14.327747696132253, 14.354508911206949, 14.331218676773712, 14.251150311763046, 14.088257739618454, 13.76917452212827, 13.101581593462868]
    @test invweight_location(0:10,0:10) ≈ [2.3478642777118957, 2.950680392065609, 3.2200556525889006, 3.348975353894235, 3.4103053158720016, 3.4297605864726126, 3.413776296513463, 3.3555688532471923, 3.2289225181937002, 2.9601916238626513, 2.3566734930529947]
    @test invweight_age(0:10) ≈ [10.744908100409873, 10.808982898265619, 10.859335221424754, 10.895581458515855, 10.91744238026025, 10.924748324734335, 10.91744238026025, 10.895581458515855, 10.859335221424754, 10.808982898265619, 10.744908100409873]
    @test invweight(3:7, 3:7, 3:7, spatialscale=1:3, agescale=20:10:40) ≈ [6.455459194232295 6.4954174826726305 6.50970747669241; 7.221911822727031 7.261870111167367 7.276160105187145; 7.829876642142075 7.86983493058241 7.8841249246021885;;; 6.795239175454396 6.815587062014095 6.822796065519029; 7.814257428828995 7.8346053153886945 7.841814318893627; 8.464965591963413 8.485313478523112 8.492522482028043;;; 6.866607765314396 6.880327963335503 6.88516052627444; 7.9787064574491104 7.992426655470217 7.9972592184091535; 8.67354089449909 8.687261092520199 8.692093655459136;;; 6.796516961665698 6.816864848225396 6.82407385173033; 7.81538067459392 7.835728561153617 7.84293756465855; 8.46552963870305 8.485877525262747 8.49308652876768;;; 6.456852342594404 6.496810631034739 6.511100625054519; 7.223828686779423 7.263786975219759 7.278076969239537; 7.831650416268091 7.871608704708425 7.885898698728204]

    # Test NaN-ful cases
    lat, lon, age = [0:10..., NaN], [0:10..., NaN], [0:10..., NaN]
    @test invweight(lat, lon, age) ≈ [13.092772378121769, 13.759663290331229, 14.079390874013654, 14.244556812410089, 14.327747696132253, 14.354508911206949, 14.331218676773712, 14.251150311763046, 14.088257739618454, 13.76917452212827, 13.101581593462868, Inf]
    @test invweight_location(lat, lon) ≈ [2.3478642777118957, 2.950680392065609, 3.2200556525889006, 3.348975353894235, 3.4103053158720016, 3.4297605864726126, 3.413776296513463, 3.3555688532471923, 3.2289225181937002, 2.9601916238626513, 2.3566734930529947, Inf]
    @test invweight_age(age) ≈ [10.744908100409873, 10.808982898265619, 10.859335221424754, 10.895581458515855, 10.91744238026025, 10.924748324734335, 10.91744238026025, 10.895581458515855, 10.859335221424754, 10.808982898265619, 10.744908100409873, Inf]

    @test resamplingprobability(invweight(lat, lon, age)) ≈ [0.17709468716990867, 0.16997032242046156, 0.1667541264064975, 0.16513990970091127, 0.16433863731114293, 0.1640825308282905, 0.16430537469842987, 0.16507611788372856, 0.16666666666666666, 0.16987285789479709, 0.17699668840192856, 0.0]
    @test resamplingprobability(invweight_location(lat, lon)) ≈ [0.2157181923483138, 0.1795609083825818, 0.16704894086177344, 0.1616579035210186, 0.15921356093554181, 0.158453529370978, 0.1590774311384848, 0.16139152082816496, 0.16666666666666666, 0.17908729309480173, 0.21508527491483678, 0.0]
    @test resamplingprobability(invweight_age(age)) ≈ [0.16814313324769609, 0.1673131616109112, 0.16666666666666666, 0.16620436987522427, 0.16592678603901323, 0.16583422381194865, 0.16592678603901323, 0.16620436987522427, 0.16666666666666666, 0.1673131616109112, 0.16814313324769609, 0.0]

## --- bin_bsr

    x = 0:100; y = 0:100
    xmin = 0; xmax = 100; nbins = 5
    step = (xmax-xmin)/nbins
    (c,m,e) = bin_bsr(x, y, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.5)
    @test isapprox(e, [1.17, 1.21, 1.23, 1.26, 1.28], atol=0.5)

    (c,m,e) = bin_bsr(x, y, xmin:step:xmax, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.5)
    @test isapprox(e, [1.17, 1.21, 1.23, 1.26, 1.28], atol=0.5)

    # Upper and lower CIs
    (c,m,el,eu) = bin_bsr(nanbinmean!, x, y, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.5)
    @test isapprox(el, [2.29, 2.38, 2.41, 2.49, 2.51], atol=1.0)
    @test isapprox(eu, [2.3, 2.37, 2.42, 2.51, 2.51], atol=1.0)

    # Upper and lower CIs
    (c,m,el,eu) = bin_bsr(nanbinmean!, x, y, xmin:step:xmax, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.5)
    @test isapprox(el, [2.29, 2.38, 2.41, 2.49, 2.51], atol=1.0)
    @test isapprox(eu, [2.3, 2.37, 2.42, 2.51, 2.51], atol=1.0)

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
    @test isapprox(m, [10.04, 29.94, 49.94, 69.92, 89.83], atol=0.5)
    @test isapprox(e, [1.17, 1.21, 1.23, 1.26, 1.28], atol=0.5)

    # with 2-D array (matrix) of y data
    y = repeat(0:100, 1, 4)
    y_sigma = ones(101,4)
    (c,m,e) = bin_bsr(x, y, xmin, xmax, nbins, x_sigma=ones(101), y_sigma=y_sigma)
    @test c == 10.0:20.0:90.0
    @test isapprox(m, repeat([10.04, 29.94, 49.94, 69.92, 89.83], 1, 4), atol=0.5)
    @test isapprox(e, repeat([1.17, 1.21, 1.23, 1.26, 1.28], 1, 4), atol=0.5)

    (c,m,el,eu) = bin_bsr(x, y, xmin, xmax, nbins, x_sigma=ones(101), y_sigma=y_sigma, sem=:CI)
    @test c == 10.0:20.0:90.0
    @test isapprox(m, repeat([10.04, 29.94, 49.94, 69.92, 89.83], 1, 4), atol=0.5)
    @test isapprox(el, repeat([2.29, 2.38, 2.41, 2.49, 2.51], 1, 4), atol=1.0)
    @test isapprox(eu, repeat([2.3, 2.37, 2.42, 2.51, 2.51], 1, 4), atol=1.0)


## -- bin_bsr_ratios

    x = 0:100; num = 0:100; denom=reverse(num)
    xmin = 0; xmax = 100; nbins = 5
    step = (xmax-xmin)/nbins
    (c,m,el,eu) = bin_bsr_ratios(x, num, denom, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.33, 8.99], rtol=0.1)
    @test isapprox(el, [0.03, 0.05, 0.09, 0.26, 2.11], rtol=0.4)
    @test isapprox(eu, [0.03, 0.05, 0.1, 0.29, 3.03], rtol=0.4)

    (c,m,el,eu) = bin_bsr_ratios(x, num, denom, xmin:step:xmax, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.33, 8.99], rtol=0.1)
    @test isapprox(el, [0.03, 0.05, 0.09, 0.26, 2.11], rtol=0.4)
    @test isapprox(eu, [0.03, 0.05, 0.1, 0.29, 3.03], rtol=0.4)

    # With weights
    (c,m,el,eu) = bin_bsr_ratios(x, num, denom, xmin, xmax, nbins, ones(101), x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.33, 8.99], rtol=0.1)
    @test isapprox(el, [0.03, 0.05, 0.09, 0.26, 2.11], rtol=0.4)
    @test isapprox(eu, [0.03, 0.05, 0.1, 0.29, 3.03], rtol=0.4)

    # Medians
    (c,m,el,eu) = bin_bsr_ratio_medians(x, num, denom, xmin, xmax, nbins, x_sigma=ones(101))
    @test c == 10.0:20.0:90.0
    @test isapprox(m, [0.11, 0.43, 1.0, 2.34, 9.25], rtol=0.1)
    @test isapprox(el, [0.05, 0.08, 0.15, 0.4, 3.1], rtol=0.5)
    @test isapprox(eu, [0.05, 0.09, 0.17, 0.51, 6.42], rtol=0.5)

## --- Constant silica reference model

    N = 10000
    SiO2 = rand(N) * 40 .+ 40
    Age = rand(N) * 3800
    Y = 1 .+ rand(N) .* (SiO2 .- 40)./10 .* (4000 .- Age)./1000
    ds = (; Y=Y, Age=Age, SiO2=SiO2, Y_sigma=Y*0.05, Age_sigma=Age*0.05)

    c,m,el,eu = constproportion(bin_bsr_means, ds, "Age", "Y", 0, 4000, 8)
    @test c == 250:500:3750
    @test round.(m, sigdigits=2) ≈  [4.9, 4.4, 3.8, 3.4, 2.8, 2.3, 1.8, 1.4] rtol=0.2
    @test round.(el, sigdigits=2) ≈ [0.12, 0.1, 0.091, 0.08, 0.063, 0.043, 0.029, 0.021] rtol=0.2
    @test round.(eu, sigdigits=2) ≈ [0.12, 0.1, 0.092, 0.074, 0.062, 0.043, 0.03, 0.023] rtol=0.2

    Num = 10 .+ rand(N) .* (SiO2 .- 40)./10 .* (4000 .- Age)./1000
    Denom = 100 .- rand(N) .* (SiO2 .- 40)./10 .* Age./1000
    ds = (; ds..., Num=Num, Denom=Denom, Num_sigma=0.05Num, Denom_sigma=0.05Denom)

    c,m,el,eu = constproportion(bin_bsr_ratio_medians, ds, "Age", "Num", "Denom", 0, 4000, 8)
    @test c == 250:500:3750
    @test round.(m, sigdigits=2) ≈ [0.14, 0.13, 0.13, 0.12, 0.12, 0.12, 0.11, 0.11] rtol=0.2
    @test round.(el, sigdigits=2) ≈ [0.0019, 0.0018, 0.0016, 0.0014, 0.0011, 0.00087, 0.00073, 0.00084] rtol=0.2
    @test round.(eu, sigdigits=2) ≈ [0.002, 0.0018, 0.0017, 0.0013, 0.0011, 0.00089, 0.00076, 0.00083] rtol=0.2


## --- Monte Carlo interpolation/fitting

    (c,m) = mcfit(0:11, ones(12), 0:11, ones(12), 1, 10, 10)
    @test c == 1:10
    @test isapprox(m, [1.15, 2.02, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.98, 9.85], atol = 0.25)

## --- Downsampling

    @test downsample(1:100, 10) == 10:10:100
    A = reshape(1:100,10,10)
    @test downsample(A, 2) == [12 32 52 72 92; 14 34 54 74 94; 16 36 56 76 96; 18 38 58 78 98; 20 40 60 80 100]
    @test downsample(collect(A), 2) == [12 32 52 72 92; 14 34 54 74 94; 16 36 56 76 96; 18 38 58 78 98; 20 40 60 80 100]

## ---
