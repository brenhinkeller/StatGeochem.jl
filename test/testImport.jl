## --- String parsing functions

    @test parsedlm("1,2,3\n4,5,6\n7,8,9\n", ',', Float64) == reshape(1:9,3,3)'
    @test parsedlm("1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", ',', Int64) == reshape(1:16,4,4)'

    A = delim_string_function(x -> delim_string_parse(x, ',', Float32), "1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", '\n', Array{Float32,1})
    @test isa(A, Array{Array{Float32,1},1})
    @test all([A[i][j] == (i-1)*4 + j for i=1:4, j=1:4])

    A = delim_string_function(x -> delim_string_parse(x, ',', Int64, merge=true, undefval=0),
        "\n1,2,3,,4\n5,6,,7,8\n9,10,,,,11,12\n\n\n13,14,15,16", '\n', Array{Int64,1}, merge=true)
    @test all([A[i][j] == (i-1)*4 + j for i=1:4, j=1:4])

## --- Elementify/unelementify functions

    elements = string.(permutedims(unique(rand("abcdefghijklmnopqrstuvwxyz",11))))
    data = vcat(elements, hcat(randn(1000, length(elements)-1), string.(rand("abcdefghijklmnopqrstuvwxyz0123456789",1000))))
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
    @test importdataset("dictdataset.csv", ',', importas=:Dict, mindefinedcolumns=2) == datadict


## --  Normalization functions

    # Renormalization functions on arrays
    dataarray = unelementify(datadict, findnumeric=true, floatout=true)
    dataarray[rand(1:length(dataarray),100)] .= NaN
    renormalize!(dataarray, total=100)
    @test nansum(dataarray) ≈ 100
    renormalize!(dataarray, dim=1, total=100)
    @test all(nansum(dataarray, dims=1) .≈ 100)
    renormalize!(dataarray, dim=2, total=100)
    @test all(nansum(dataarray, dims=2) .≈ 100)

    # Renormalization functions on NamedTuple-based dataset
    datatuple = elementify(unelementify(datatuple, findnumeric=true), importas=:Tuple)
    renormalize!(datatuple, total=100)
    @test all(sum(unelementify(datatuple, floatout=true),dims=2) .≈ 100)

    # Renormalization functions on Dict-based dataset
    datadict = elementify(unelementify(datadict, findnumeric=true), importas=:Dict)
    renormalize!(datadict, datadict["elements"], total=100)
    @test all(sum(unelementify(datadict, floatout=true),dims=2) .≈ 100)

    # Internal standardization functions
    @test isnan(StatGeochem.floatify("asdf"))
    @test StatGeochem.floatify("12345", Float64) === 12345.0
    @test StatGeochem.floatify("12345", Float32) === 12345f0
    @test StatGeochem.floatify(12345, Float64) === 12345.0
    @test StatGeochem.floatify(12345, Float32) === 12345f0

    @test isa(StatGeochem._columnformat(["asdf","qwer","zxcv"], false), Array{String,1})
    @test isa(StatGeochem._columnformat([1f0, 2f0, 3f0], false), Array{Float32,1})
    @test isa(StatGeochem._columnformat([1., 2., 3.], false), Array{Float64,1})
    @test isa(StatGeochem._columnformat([0x01,0x02,0x03], false), Array{UInt8,1})
    @test isa(StatGeochem._columnformat([1,2,3], false), Array{Int64,1})
    @test all(StatGeochem._columnformat([0x01,2,"3"], false) .=== [0x01,2,"3"])
    @test StatGeochem._columnformat([0x01,2,"3"], true) == [1,2,3]

## --- Concatenating and merging datasets

    d2 = concatenatedatasets(datadict, datadict)
    @test isa(d2, Dict)

    d2array = unelementify(d2, floatout=true)
    @test isa(d2array, Array{Float64,2})
    @test size(d2array) == (2000, length(datadict["elements"]))

## ---
