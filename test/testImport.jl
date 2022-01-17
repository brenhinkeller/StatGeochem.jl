## --- String parsing functions

    @test parsedlm("1,2,3\n4,5,6\n7,8,9\n", ',', Float64) == reshape(1:9,3,3)'
    @test parsedlm("1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", ',', Int64) == reshape(1:16,4,4)'

    A = delim_string_function(x -> delim_string_parse(x, ',', Float32),
        "1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", '\n', Array{Float32,1})
    @test isa(A, Array{Array{Float32,1},1})
    @test all([A[i][j] == (i-1)*4 + j for i=1:4, j=1:4])

    A = delim_string_function(x -> delim_string_parse(x, ',', Int64, merge=true, undefval=0),
        "\n1,2,3,,4\n5,6,,7,8\n9,10,,,,11,12\n\n\n13,14,15,16", '\n', Array{Int64,1}, merge=true)
    @test all([A[i][j] == (i-1)*4 + j for i=1:4, j=1:4])

## --- Elementify/unelementify functions

    elements = ["U" "Lv" "Te" "O" "Er" "W" "Re" "j" "asdf" "Zr" "Al" "S" "K" "V" "N" "Ga" "I"]
    data = vcat(elements, hcat(rand(1000, length(elements)-1), string.(rand("abcdefghijklmnopqrstuvwxyz0123456789",1000))))
    datatuple = elementify(data,importas=:Tuple)::NamedTuple
    datadict = elementify(data,importas=:Dict)::Dict

    @test isa(datatuple, NamedTuple)
    @test unelementify(datatuple) == data
    @test isa(datadict, Dict)
    @test unelementify(datadict) == data

    # Test adding or averaging option for numeric elements
    addtest = ["a" "b" "a";1 2 3]
    avg = elementify(addtest, importas=:Dict)
    add = elementify(addtest, importas=:Dict, sumduplicates=true)
    @test avg["elements"] == avg["elements"]
    @test avg["a"] == 2
    @test add["a"] == 4

## --- Import / export functions

    @test exportdataset(datatuple, "tupledataset.csv", ',') == nothing
    @test importdataset("tupledataset.csv", ',', importas=:Tuple) == datatuple

    @test exportdataset(datatuple, "tupledataset.csv", ',', digits=6) == nothing
    @test importdataset("tupledataset.csv", ',', importas=:Tuple).Lv == round.(datatuple.Lv, digits=6)

    @test exportdataset(datatuple, "tupledataset.csv", ',', sigdigits=5) == nothing
    @test importdataset("tupledataset.csv", ',', importas=:Tuple).Lv == round.(datatuple.Lv, sigdigits=5)

    @test exportdataset(datadict, datadict["elements"], "dictdataset.csv", ',') == nothing
    @test importdataset("dictdataset.csv", ',', importas=:Dict, mindefinedcolumns=2) == datadict


## --  Normalization functions

    dataarray = rand(1000, length(elements))
    data = vcat(elements, dataarray)
    datatuple = elementify(data,importas=:Tuple)::NamedTuple
    datadict = elementify(data,importas=:Dict)::Dict

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
    renormalize!(datadict, datadict["elements"], total=100.)
    @test all(sum(unelementify(datadict, floatout=true),dims=2) .≈ 100)

    # Internal standardization functions
    @test isnan(StatGeochem.floatify("asdf"))
    @test StatGeochem.floatify("12345", Float64) === 12345.0
    @test StatGeochem.floatify("12345", Float32) === 12345f0
    @test StatGeochem.floatify(12345, Float64) === 12345.0
    @test StatGeochem.floatify(12345, Float32) === 12345f0

    @test isa(StatGeochem.columnformat(["asdf","qwer","zxcv"], false), Array{String,1})
    @test isa(StatGeochem.columnformat([1f0, 2f0, 3f0], false), Array{Float32,1})
    @test isa(StatGeochem.columnformat([1., 2., 3.], false), Array{Float64,1})
    @test isa(StatGeochem.columnformat([0x01,0x02,0x03], false), Array{UInt8,1})
    @test isa(StatGeochem.columnformat([1,2,3], false), Array{Int64,1})
    @test all(StatGeochem.columnformat([0x01,2,"3"], false) .=== [0x01,2,"3"])
    @test StatGeochem.columnformat([0x01,2,"3"], true) == [1,2,3]
    @test StatGeochem.columnformat(["asdf","qwer","zxcv"], true) == ["asdf","qwer","zxcv"]
    @test StatGeochem.columnformat(["","","zxcv"], true) == ["","","zxcv"]
    @test isequal(StatGeochem.columnformat(["","","5"], true), [NaN, NaN, 5.0])

    @test StatGeochem.isnumeric(missing) == false
    @test StatGeochem.nonnumeric(missing) == false
    @test StatGeochem.isnumeric("") == false
    @test StatGeochem.nonnumeric("") == false
    @test StatGeochem.isnumeric("5") == true
    @test StatGeochem.nonnumeric("5") == false
    @test StatGeochem.isnumeric('x') == false
    @test StatGeochem.nonnumeric('x') == true
    @test StatGeochem.isnumeric(NaN) == true
    @test StatGeochem.nonnumeric(NaN) == false

    @test isequal(StatGeochem.emptys(Any,3), [missing, missing, missing])
    @test isequal(StatGeochem.emptys(String,3), ["", "", ""])
    @test all(StatGeochem.emptys(Float16,3) .=== Float16[NaN, NaN, NaN])
    @test all(StatGeochem.emptys(Float64,3) .=== [NaN, NaN, NaN])
    @test all(StatGeochem.emptys(Int64,3) .=== [NaN, NaN, NaN])


## --- Concatenating and merging datasets

    d2 = concatenatedatasets(datadict, datadict)
    @test isa(d2, Dict)

    d2array = unelementify(d2, floatout=true)
    @test isa(d2array, Array{Float64,2})
    @test size(d2array) == (2000, length(datadict["elements"]))

    A = ["La" "Ce" "Pr" "ID"; 1.5 1.1 1.0 "x9"; 3.7 2.9 2.5 "SJ21-12"]
    B = ["La" "Yb"; 1.5 1.1; 1.0 3.7; 2.9 2.5]
    a = elementify(A, importas=:Tuple)
    b = elementify(B, importas=:Tuple)
    d = concatenatedatasets(a,b)
    @test isa(d, NamedTuple)
    @test isequal(d.La, [1.5, 3.7, 1.5, 1.0, 2.9])
    @test isequal(d.Yb, [NaN, NaN, 1.1, 3.7, 2.5])
    @test isequal(d.ID, ["x9", "SJ21-12", "", "", ""])

    darray = unelementify(d, floatout=true)
    @test isa(darray, Array{Float64,2})
    @test size(darray) == (5, length(keys(d)))

## ---
