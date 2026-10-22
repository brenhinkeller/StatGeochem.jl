cd(@__DIR__)
date = "15-OCT-2026"
for mineral in ("Albite", "Anorthite", "Amphibole", "Biotite", "Clinopyroxene", "Ilmenite", "Magnetite", "Orthopyroxene", "Olivine", "Garnet", "Spinel")
    ds = importdataset("$mineral.tsv", importas=:Tuple)

    parsed = (;
        Rock_Type = String[],
        Mineral	= String[],
        Z = Int[],
        Elem = String[],
        Value = Float64[],
        SD = Float64[],
        Low = Float64[],
        High = Float64[],
        Kd_Type	= String[],
        Reference = String[],
        Definition = String[],
        Notes = String[],
        Inserted = String[],
        SiO2 = Float64[],
    )

    for i in eachindex(ds.SiO2)
        if 0 < ds.SiO2[i] < 100
            for e in keys(atomicnumber)
                if haskey(ds, e) && (0 < ds[e][i])
                    descriptor = if ds.Kd_Type[i] == "Experimental"
                        "Synthetic ($(strip(ds.Sample[i]))), $(ds.P_GPa[i]) GPa, $(ds.T_C[i]) C" * (isempty(ds.fO2_condition[i]) ? "" : ", $(ds.fO2_condition[i])")
                    else
                        strip(ds.Sample[i])
                    end
                    push!(parsed.Rock_Type, descriptor)
                    detailedmineral = (haskey(ds, :Mineral) && !isempty(ds.Mineral[i])) ? ds.Mineral[i] : mineral
                    if (detailedmineral == "Albite") || (detailedmineral == "Anorthite")
                        push!(parsed.Mineral, "Plagioclase (An$(round(Int, ds._An[i])))")
                    else
                        push!(parsed.Mineral, detailedmineral)
                    end
                    push!(parsed.Z, atomicnumber[e])
                    push!(parsed.Elem, String(e))
                    push!(parsed.Value, ds[e][i])
                    e_sigma = Symbol("$(e)_sigma")
                    push!(parsed.SD, (haskey(ds, e_sigma) ? ds[e_sigma][i] : NaN))
                    push!(parsed.Low, NaN)
                    push!(parsed.High, NaN)
                    push!(parsed.Kd_Type, ds.Kd_Type[i])
                    push!(parsed.Reference, ds.Citation[i])
                    push!(parsed.Definition, "Solid-Melt")
                    push!(parsed.Notes, ds.Compilation[i]) 
                    push!(parsed.Inserted, date)
                    push!(parsed.SiO2, ds.SiO2[i])
                end
            end
        end
    end

    exportdataset(parsed, "$(mineral)_parsed.tsv")
end