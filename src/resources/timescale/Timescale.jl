gts = importdataset(joinpath(moduleresourcepath,"timescale","gts_intervals.tsv"), '\t', importas=:Tuple)

e = (:Age_min, :Age_min_sigma, :Age_max, :Age_max_sigma)

const timescale = NamedTuple{e}(((Dict{String,Float64}() for _ in e)...,))
for i in eachindex(gts.Name)
    timescale.Age_max[gts.Name[i]] = gts.Age_max[i]
    timescale.Age_max_sigma[gts.Name[i]] = gts.Age_max_sigma[i]
    timescale.Age_min[gts.Name[i]] = gts.Age_min[i]
    timescale.Age_min_sigma[gts.Name[i]] = gts.Age_min_sigma[i]
end
export timescale 

