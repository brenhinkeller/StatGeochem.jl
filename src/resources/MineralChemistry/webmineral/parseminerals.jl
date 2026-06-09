## --- Clean up CSV files
    
    cd(@__DIR__)
    # # Only need to run once, should already be set
    # for file in filter(d->contains(d,".csv"), readdir(@__DIR__))
    #     run(`sed -i.bak -r -e 's/([0-9])%[A-Z][a-z]?/\1/' -e 's/[[:space:]!?*]*,/,/g' -e 's/[[:space:]]*"/"/g' $file`)
    # end

## --- Organize as per-mineral dictionaries
    using StatGeochemBase

    minerals = Dict{String,Dict{Symbol,Float64}}()
    for file in filter(d->contains(d,".csv"), readdir(@__DIR__))
        elem = replace(file, ".csv"=>"")
        data = importdataset(file)
        names = data["Mineral Name"]
        concentrations = data["%"*elem]
        if length(concentrations) > 1
            for i in eachindex(names, concentrations)
                if !haskey(minerals, names[i])
                    minerals[names[i]] = Dict{Symbol,Float64}()
                end
                minerals[names[i]][Symbol(elem)] = concentrations[i]
            end
        else
            if !haskey(minerals, names)
                minerals[names] = Dict{Symbol,Float64}()
            end
            minerals[names][Symbol(elem)] = concentrations
        end
end

## -- Fix minerals with missing REE

    minerals["Ashanite"] = merge(minerals["Ashanite"], Dict(:Y=>2.5)) # By differnce
    minerals["Arsenoflorencite-(La)"] = merge(minerals["Arsenoflorencite-(La)"], Dict(:Nd=>2.55, :Th=>2.55)) # Assign remainder of bulk REE to Nd and Th
    minerals["Arsenoflorencite-(Ce)"] = merge(minerals["Arsenoflorencite-(Ce)"], Dict(:Nd=>2.45, :Th=>2.45)) # Assign remainder of bulk REE to Nd and Th
    minerals["Calciobetafite"] = merge(minerals["Calciobetafite"], Dict(:Ce=>3.95, :Y=>3.35)) # Rescaled from yttrobetafite
    minerals["Cerotungstite-(Ce)"] = merge(minerals["Cerotungstite-(Ce)"], Dict(:Th=>1.07)) # By difference
    minerals["Cheralite-(Ce)"] = merge(minerals["Cheralite-(Ce)"], Dict(:Pr=>1.15, :Sm=>1.85)) # Rescaled to Chondritic abundances for the remainig LREE
    minerals["Dissakisite-(Ce)"] = merge(minerals["Dissakisite-(Ce)"], Dict(:Y=>3.03)) # By difference; the only REE omitted
    minerals["Dualite"] = merge(minerals["Dualite"], Dict(:La=>0.77, :Ce=>1.61, :Nd=>0.69)) # From https://handbookofmineralogy.org/pdfs/Dualite.pdf
    minerals["Fergusonite-(Nd)"] = merge(minerals["Fergusonite-(Nd)"], Dict(:Pr=>5.99, :Sm=>6.39, :Y=>3.78)) # Rescaled from Fergusonite-beta-(Nd)
    minerals["Fergusonite-beta-(Ce)"] = merge(minerals["Fergusonite-beta-(Ce)"], Dict(:Th=>2.44, :Sm=>1.6, :Y=>0.96)) # Th from REE scaling, Sm & Y rescaled from Fergusonite-beta-(Nd)
    minerals["Ferrokentbrooksite"] = merge(minerals["Ferrokentbrooksite"], Dict(:Tb => 0.21, :Dy => 1.23, :Yb => 0.13, :Ho => 0.17, :Tm => 0.12, :Sm => 0.47, :Gd => 0.95, :Pr => 0.05, :Er => 0.41, :Ce => 0.18, :Nd => 0.48)) # Rescaled to match bulk "REE", using pattern from Samarskite-(Y)
    minerals["Gadolinite-(Ce)"] = merge(minerals["Gadolinite-(Ce)"], Dict(:Th=>2.28, :Pr=>2.07, :Nd=>10.19, :Sm=>3.34)) # Rescaled to Chondritic abundances for the remainig LREE
    minerals["Hainite"] = merge(minerals["Hainite"], Dict(:Dy => 0.16, :La => 0.13, :Yb => 0.15, :Pr => 0.04, :Er => 0.12, :Lu => 0.03, :Y => 0.88, :Ce => 0.32, :Nd => 0.15)) # Rescaled from RRUFF
    minerals["Hellandite-(Ce)"] = merge(minerals["Hellandite-(Ce)"], Dict(:Th=>2.53)) # By difference
    minerals["Hellandite-(Y)"] = merge(minerals["Hellandite-(Y)"], Dict(:Th=>2.7)) # By difference
    minerals["Khristovite-(Ce)"] = merge(minerals["Khristovite-(Ce)"], Dict(:Th=>1.4, :Y=>1.85)) # Th from REE scaling, Y by difference
    minerals["Murataite"] = merge(minerals["Murataite"], Dict(:Tb => 0.41, :Dy => 2.41, :Yb => 0.25, :Ho => 0.34, :Tm => 0.24, :Sm => 0.93, :Gd => 1.86, :Pr => 0.1, :Er => 0.8, :Ce => 0.36, :Nd => 0.94)) # Rescaled to match bulk "REE", using pattern from Samarskite-(Y)
    minerals["Okanoganite-(Y)"] = merge(minerals["Okanoganite-(Y)"], Dict(:Tb => 0.4729, :Dy => 3.2394, :Yb => 2.1044, :Ho => 0.733, :Tm => 0.33103, :Sm => 1.9625, :Gd => 2.601, :Pr => 1.2295, :Er => 2.1044, :Ce => 16.102, :Nd => 6.0295)) # Rescaled to match bulk "REE", using pattern from Hellandite-(Y)
    minerals["Plumbobetafite"] = merge(minerals["Plumbobetafite"], Dict(:Ce=>1.41, :Y=>1.19, :Th=>0.78)) # Rescaled from yttrobetafite
    minerals["Perrierite-(Ce)"] = merge(minerals["Perrierite-(Ce)"], Dict(:Pr=>0.64, :Nd=>3.13, :Sm=>1.02)) # Remaining LREE scaled to chondritic
    minerals["Samarskite-(Y)"] = merge(minerals["Samarskite-(Y)"], Dict(:Tb => 0.68, :Dy => 4.02, :Yb => 0.41, :Ho => 0.56, :Tm => 0.4, :Sm => 1.55, :Gd => 3.1, :Pr => 0.16, :Er => 1.34, :Ce => 0.6, :Nd => 1.58)) # Rescaled to match bulk "REE", using pattern from RRUFF
    minerals["Saryarkite-(Y)"] = merge(minerals["Saryarkite-(Y)"], Dict(:Tb => 0.29, :Dy => 1.74, :Yb => 0.18, :Ho => 0.24, :Tm => 0.17, :Sm => 0.67, :Gd => 1.34, :Pr => 0.07, :Er => 0.58, :Ce => 0.26, :Nd => 0.68)) # Rescaled to match bulk "REE", using pattern from Samarskite-(Y)
    minerals["Semenovite"] = merge(minerals["Semenovite"], Dict(:Nd=>2.54, :Sm=>0.83)) # Extrapolating from La and Ce, assuming a chondritic ratio
    minerals["Titanite"] = merge(minerals["Titanite"], Dict(:Tb => 0.02, :Dy => 0.11, :Yb => 0.04, :Ho => 0.02, :Tm => 0.01, :Sm => 0.14, :Y => 0.53, :Gd => 0.13, :La => 0.46, :Pr => 0.15, :Er => 0.06, :Ce => 1.31, :Eu => 0.03, :Nd => 0.63)) # Rescaled to match bulk "REE", using garnet-absent average from https://doi.org/10.1016/S0009-2541(01)00355-2
    minerals["Thomasclarkite-(Y)"] = merge(minerals["Thomasclarkite-(Y)"], Dict(:Th=>5.25)) # By difference
    minerals["Tombarthite-(Y)"] = merge(minerals["Tombarthite-(Y)"], Dict(:Tb => 0.91, :Dy => 5.4, :Yb => 0.55, :Ho => 0.76, :Tm => 0.54, :Sm => 2.08, :Gd => 4.17, :Pr => 0.22, :Er => 1.8, :Ce => 0.8, :Nd => 2.12)) # Rescaled to match bulk "REE", using pattern from Samarskite-(Y)
    minerals["Vyuntspakhkite-(Y)"] = merge(minerals["Vyuntspakhkite-(Y)"], Dict(:Th => 3.09, :Dy => 2.27, :Er => 4.03, :Ho => 0.87, :Tm => 0.83, :Lu => 1.54)) # Rescaled to fit bulk "REE", using pattern from zircon
    minerals["Yttrotungstite-(Y)"] = merge(minerals["Yttrotungstite-(Y)"], Dict(:Tb => 0.43, :Dy => 2.56, :Yb => 0.26, :Ho => 0.36, :Tm => 0.26, :Sm => 0.98, :Gd => 1.97, :Pr => 0.1, :Er => 0.85, :Ce => 0.38, :Nd => 1.0)) # Rescaled to match bulk "REE", using pattern from Samarskite-(Y)
    minerals["Zircon"] = merge(minerals["Zircon"], Dict(:Th => 0.18, :Dy => 0.13, :Yb => 0.46, :U => 1.25, :Er => 0.23, :Ho => 0.05, :Tm => 0.05, :Lu => 0.09, :Y => 1.34)) # Rescaled to match bulk "REE", using pattern from https://doi.org/10.1016/j.chemgeo.2019.01.006, sample EAF60


## --- Check for approximate normalization

    for name in keys(minerals)
        mineral = minerals[name]
        total = sum(k->mineral[k], keys(mineral))
        if !(99 < total < 103)
            @info name, total
            @info mineral
            println()
        end
    end

## --- Renormalize everything to 100%

    for name in keys(minerals)
        mineral = minerals[name]
        normconst = sum(k->mineral[k], keys(mineral))
        for k in keys(mineral)
            mineral[k] = round(mineral[k]*100/normconst, digits=2)
        end
    end

## --- Save results

    f = open(joinpath("..","mineralcompositions.jl"), "a")
    print(f, "\nmineralcompositions = $minerals\nexport mineralcompositions\n")
    close(f)

## --- End of File