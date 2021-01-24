## --- Elementify/unelementify functions
    elements = string.(permutedims(unique(rand("abcdefghijklmnopqrstuvwxyz",10))))
    data = vcat(elements, randn(1000, length(elements)))
    datatuple = elementify(data,importas=:Tuple)
    datadict = elementify(data,importas=:Dict)

    @test isa(datatuple, NamedTuple)
    @test unelementify(datatuple) == data
    @test isa(datadict, Dict)
    @test unelementify(datadict) == data

## --  Normalization functions

    # Renormalization functions on arrays
    dataarray = unelementify(datatuple, floatout=true)
    renormalize!(dataarray, total=100)
    @test sum(dataarray) ≈ 100
    renormalize!(dataarray, dim=1, total=100)
    @test all(sum(dataarray, dims=1) .≈ 100)
    renormalize!(dataarray, dim=2, total=100)
    @test all(sum(dataarray, dims=2) .≈ 100)

    # Renormalization functions on NamedTuple-based dataset
    renormalize!(datatuple, total=100)
    @test all(sum(unelementify(datatuple, floatout=true),dims=2) .≈ 100)

    # Renormalization functions on Dict-based dataset
    renormalize!(datadict, datadict["elements"], total=100)
    @test all(sum(unelementify(datadict, floatout=true),dims=2) .≈ 100)

## ---
