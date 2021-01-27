## --- Resampling.jl

    index = Array{Int64}(undef,1000)
    resampled = Array{Int64}(undef,1000)
    @test bsr!(resampled,index,1:10,0,0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,0,fill(0.5,1000)) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),0.5) |> unique ⊆ 1:10
    @test bsr!(resampled,index,1:10,fill(0,10),fill(0.5,1000)) |> unique ⊆ 1:10

    resampled = Array{Int64}(undef,1000)
    bsr!(resampled, index, ones(10), 0.5, 1000, 0.5)
    @test isapprox(mean(reasampled), 1, atol=0.15)
    bsr!(resampled, index, ones(10), 0.5, 1000, fill(0.5,1000))
    @test isapprox(mean(reasampled), 1, atol=0.15)
    bsr!(resampled, index, ones(10), fill(0.5,10), 1000, 0.5)
    @test isapprox(mean(reasampled), 1, atol=0.15)
    bsr!(resampled, index, ones(10), fill(0.5,10), 1000, fill(0.5,1000))
    @test isapprox(mean(reasampled), 1, atol=0.15)

    @test bsresample(1:10,fill(0.5,10),1000,0.5)::Array{Float64} |> length == 1000
    @test bsresample(1:10,fill(0.5,10),1000,fill(0.5,10))::Array{Float64} |> length == 1000

    @test bincounts(1:100,0,100,10) == (5:10:95, fill(10,10))
    @test binmeans(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(0.9574271077563381,10))
    @test binmedians(1:100,1:100,0,100,10) == (5:10:95, 5.5:10:95.5, fill(1.1720982147414096,10))

    @test randsample(1:10,1000)::Array{Int64} |> length == 1000
    @test unique(randsample(1:10,1000)) ⊆ 1:10

## --- Invweight

    @test invweight(0:10, 0:10, 0:10) ≈ [13.092772378121769, 13.759663290331229, 14.079390874013654, 14.244556812410089, 14.327747696132253, 14.354508911206949, 14.331218676773712, 14.251150311763046, 14.088257739618454, 13.76917452212827, 13.101581593462868]
    @test invweight_location(0:10, 0:10) ≈ [2.3478642777118957, 2.950680392065609, 3.2200556525889006, 3.348975353894235, 3.4103053158720016, 3.4297605864726126, 3.413776296513463, 3.3555688532471923, 3.2289225181937002, 2.9601916238626513, 2.3566734930529947]
    @test invweight_age(0:10) ≈ [10.744908100409873, 10.808982898265619, 10.859335221424754, 10.895581458515855, 10.91744238026025, 10.924748324734335, 10.91744238026025, 10.895581458515855, 10.859335221424754, 10.808982898265619, 10.744908100409873]

## --- Monte Carlo interpolation/fitting

    (c,m) = mcfit(0:11, ones(12), 0:11, ones(12), 1, 10, 10)
    @test c == 1:10
    @test isapprox(m, [1.15, 2.02, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.98, 9.85], atol = 0.15)

## --- Downsampling

    @test downsample(1:100,10) == 10:10:100
    @test downsample(reshape(1:100,10,10),2) == [12 32 52 72 92; 14 34 54 74 94; 16 36 56 76 96; 18 38 58 78 98; 20 40 60 80 100]

## ---
