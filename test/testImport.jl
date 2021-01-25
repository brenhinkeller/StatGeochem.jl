## --- String parsing functions

    @test parsedlm("1,2,3\n4,5,6\n7,8,9\n", ',', Float64) == reshape(1:9,3,3)'
    @test parsedlm("1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", ',', Int64) == reshape(1:16,4,4)'

## --- Elementify/unelementify functions

    elements = string.(permutedims(unique(rand("abcdefghijklmnopqrstuvwxyz",11))))
    data = vcat(elements, hcat(randn(1000, length(elements)-1), string.(rand("abcdefghijklmnopqrstuvwxyz",1000))))
    datatuple = elementify(data,importas=:Tuple)::NamedTuple
    datadict = elementify(data,importas=:Dict)::Dict

    @test isa(datatuple, NamedTuple)
    @test unelementify(datatuple) == data
    @test isa(datadict, Dict)
    @test unelementify(datadict) == data

## --- Import / export functions

    @test exportdataset(datatuple, "tupledataset.csv", ',') == nothing
    @test importdataset("tupledataset.csv", ',', importas=:Tuple) == datatuple

    @test exportdataset(datadict, datadict["elements"], "dictdataset.csv", ',') == nothing
    @test importdataset("dictdataset.csv", ',', importas=:Dict) == datadict


## --  Normalization functions

    elements = string.(permutedims(unique(rand("abcdefghijklmnopqrstuvwxyz",11))))
    dataarray = rand(1000,length(elements))
    data = vcat(elements, dataarray)
    datadict = elementify(data,importas=:Dict)::Dict
    datatuple = elementify(data,importas=:Tuple)::NamedTuple

    # Renormalization functions on arrays
    dataarray[rand(1:length(dataarray),100)] .= NaN
    renormalize!(dataarray, total=100)
    @test nansum(dataarray) ≈ 100
    renormalize!(dataarray, dim=1, total=100)
    @test all(nansum(dataarray, dims=1) .≈ 100)
    renormalize!(dataarray, dim=2, total=100)
    @test all(nansum(dataarray, dims=2) .≈ 100)

    # Renormalization functions on NamedTuple-based dataset
    renormalize!(datatuple, total=100)
    @test all(sum(unelementify(datatuple, floatout=true),dims=2) .≈ 100)

    # Renormalization functions on Dict-based dataset
    renormalize!(datadict, datadict["elements"], total=100)
    @test all(sum(unelementify(datadict, floatout=true),dims=2) .≈ 100)

## --- Concatenating and merging datasets

    d2 = concatenatedatasets(datadict, datadict)
    @test isa(d2, Dict)

    d2array = unelementify(d2, floatout=true)
    @test isa(d2array, Array{Float64,2})
    @test size(d2array) == (2000, length(elements))

## ---
