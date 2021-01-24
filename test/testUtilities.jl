
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

    # Invweight...

## ---
