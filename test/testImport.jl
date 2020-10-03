elements = string.(permutedims(unique(rand("abcdefghijklmnopqrstuvwxyz",10))))
data = vcat(elements, randn(1000, length(elements)))
datatuple = elementify(data,importas=:Tuple)
datadict = elementify(data,importas=:Dict)

@test isa(datatuple, NamedTuple)
@test unelementify(datatuple) == data
@test isa(datadict, Dict)
@test unelementify(datadict) == data
