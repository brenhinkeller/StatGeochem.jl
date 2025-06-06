## --- Load required packages
    using StatGeochem
    using Plots
    using LsqFit: curve_fit

## --- Load partition coefficients from GERM into a dict

    # path = joinpath(moduleresourcepath,"PartitionCoefficients")
    path = joinpath("src/resources/","PartitionCoefficients")

    # Obtain the names of downloaded mineral files
    filenames = readdir(path)

    # Load every `.tsv` file into a struct
    pr = Dict{String,Any}()
    pr["minerals"] = String[]
    for n in filenames
        m = match(r"(.*).tsv", n)
        if ~ isnothing(m)
            name = m[1]
            push!(pr["minerals"], name)
            pr[name] = importdataset(joinpath(path, name*".tsv"), '\t', mindefinedcolumns=2, importas=:Dict)
        end
    end

    # Approximate content of various rock types, for when actual melt SiO2 not reported
    str = ["basalt","andesit","dacit","rhyolit","hawaiite","alkali basalt","basanit","benmoreite","camptonite","eclogit","garnet pyroxenit","kimberlit","komatiit","lampro","latite","leucitit","lherzolit","mugearite","pantellerite","trachyt","peridotit","phonolit","picrit","syen","morb","tholeiit","tonalit","granodiorit","granit","shoshonit","carbonatit","aplite","monzonit","synthetic"]
    val = [49,57.5,67.5,73,48.5,47,44,57,42,49,49,40,45,46.5,62.5,47,42,52,73,61,42,59,45.5,61,49,50,65,67.5,72,55,2,73,67.4,60]
    err = [5,6.5,5.5,5,3.5,4,3,4,3,5,5,4,5,11.5,11.5,5,5,5,5,9,5,5,5.5,9,5,4,8,7.5,5,5,5,5,5,20]

    # Determine the unique reference / rock type pairs
    samples = Dict{String,Vector{Tuple{String,String,String,Float64,Float64}}}()
    allsamples = Tuple{String,String,String,Float64,Float64}[]
    allelements = String[]
    for m in pr["minerals"]
        if ~haskey(pr[m],"SiO2")
            pr[m]["SiO2"] = fill(NaN, size(pr[m]["Reference"]))
        end
        if ~haskey(pr[m],"SiO2_sigma")
            pr[m]["SiO2_sigma"] = fill(NaN, size(pr[m]["Reference"]))
        end
        # Estimate SiO2 from rock type
        needssi = isnan.(pr[m]["SiO2"])
        for i ∈ eachindex(str)
            t = containsi.(pr[m]["Rock Type"], str[i]) .& needssi
            # Average together, such that "basaltic andesite" is halfway between "basalt" and "andesite"
            pr[m]["SiO2"][t] = nanmean([pr[m]["SiO2"][t] val[i]*ones(count(t))], dim=2)
            pr[m]["SiO2_sigma"][t] = nanmean([pr[m]["SiO2_sigma"][t] err[i]*ones(count(t))], dim=2)
        end

        # Save unique samples
        println(m)
        flush(stdout)
        samples[m] = collect(zip(string.(pr[m]["Reference"]), string.(pr[m]["Rock Type"]), string.(pr[m]["Mineral"]), pr[m]["SiO2"], pr[m]["SiO2_sigma"]))
        append!(allsamples, unique(samples[m]))
        any(isequal("Eu3"), pr[m]["Elem"]) && @info m

        # Save unique elements
        append!(allelements, unique(pr[m]["Elem"]))
    end
    append!(allelements, ["Eu3", "Eu2", "Ce3", "Ce4", "U4",])

    # Deduplicate
    allsamples = sort(unique(allsamples))
    allelements = sort(unique(allelements))


    # Second dict: sort the data by sample ID
    pd = Dict{String, Union{Vector, Dict{String, Union{Vector{Float64}, Vector{String}}}}}()
    pd["samples"] = allsamples
    pd["Reference"] = allsamples .|> x -> x[1]
    pd["Rock Type"] = allsamples .|> x -> x[2]
    pd["SiO2"] = allsamples .|> x -> x[4]
    pd["SiO2_sigma"] = allsamples .|> x -> x[5]
    pd["minerals"] = pr["minerals"]
    for m in pd["minerals"]
        pd[m] = Dict{String, Union{Vector{Float64}, Vector{String}}}()
        pd[m]["elements"] = allelements
        for e in allelements
            pd[m][e] = fill(NaN, length(allsamples))
            pd[m][e*"_sigma"] = fill(NaN, length(allsamples))

            t = pr[m]["Elem"] .== e
            kd = map(x->(x>0 ? log10(x) : NaN), pr[m]["Value"][t])
            kd_sigma = log10.(pr[m]["Value"][t] .+ pr[m]["SD"][t]) .- kd
            kdl = map(x->(x>0 ? log10(x) : NaN), pr[m]["Low"][t])
            kdh = map(x->(x>0 ? log10(x) : NaN), pr[m]["High"][t])

            if any(t)
                rows = findmatches(samples[m][t], allsamples)
                pd[m][e][rows] = nanmean([kd nanmean([kdl kdh], dim=2)], dim=2)
                pd[m][e*"_sigma"][rows] = nanmean([kd_sigma (kdh .- kdl)/4], dim=2)
            end
        end
    end

## --- Fit invariant-charge elements as a function of activation energy and bulk modulus

    # Minerals with well-behaved 1+ onuma diagrams
    m1p = ["Albite", "Allanite", "Amphibole", "Apatite", "Biotite", "Bridgmanite", "Olivine", "Orthopyroxene",]

    # Minerals with well-behaved 2+ onuma diagrams
    m2p = ["Albite", "Apatite", "Cordierite", "Perovskite",]

    # Minerals with well-behaved ree 3+ onuma diagrams
    # mree3p = ["Albite", "Allanite", "Amphibole", "Anorthite", "Apatite", "Baddeleyite", "Biotite", "Bridgmanite", "Chevkinite", "Clinopyroxene", "Cordierite", "Garnet", "Ilmenite", "Leucite", "Magnetite", "Melilite", "Monazite", "Muscovite", "Nepheline", "Olivine", "Orthoclase", "Orthopyroxene", "Perovskite", "Phlogopite", "Rutile", "Sphene", "Whitlockite", "Xenotime", "Zircon", "Zoisite",]
    mree3p = pr["minerals"] # All minerals

    # Minerals with well-behaved 3+ onuma diagrams
    m3p = ["Bridgmanite", "Clinopyroxene", "Garnet", "Orthopyroxene", "Perovskite", "Rutile",]

    m4p = ["Allanite", "Baddeleyite", "Bridgmanite", "Ilmenite", "Orthopyroxene", "Perovskite", "Phlogopite", "Rutile", "Sphene", "Zircon",]


    # Need more data:
    # Spinel (different kinds?) check existing REE data

    e1 = (:Li, :Na, :Ag, :K, :Rb, :Cs)
    e1r = [ionicradius[e] for e in  e1]
    e1_bonus = ()
    e1r_bonus = Float64[]

    e2 = (:Be, :Fe2, :Co, :Ni, :Mg, :Zn, :Cd, :Ca, :Eu2, :Sr, :Pb, :Ba)
    e2r = [ionicradius[e] for e in  e2]
    e2_bonus = (:Eu,)
    e2r_bonus = [ionicradius[e] for e in  (:Eu2,)]

    ree3 = (:La, :Ce3, :Pr, :Nd, :Sm, :Eu3, :Gd, :Tb, :Dy, :Ho, :Y, :Er, :Tm, :Yb, :Lu,)
    ree3r = [ionicradius[e] for e in  ree3]
    ree3_bonus = (:Ce, :Eu,)
    ree3r_bonus = [ionicradius[e] for e in (:Ce3,:Eu3,)]

    e3 = (:Cr, :Fe3, :Sc, :V3, :Lu, :Y, ree3...,)
    e3r = [ionicradius[e] for e in (:Cr3, :Fe3, :Sc, :V3, :Lu, :Y,  ree3...,)]
    e3_bonus = ()
    e3r_bonus = Float64[]

    e4 = (:Ti, :Sn, :Hf, :Zr, :Ce4, :U4, :Th)
    e4r = [ionicradius[e] for e in  e4]
    e4_bonus = (:Ce, :U,)
    e4r_bonus = [ionicradius[e] for e in (:Ce4, :U4,)]

    # Equation we're fitting to (from Blundy and Wood, 1994):
    # log10D0 + 4π a * (r0/2*(x-r0)^2 + 1/3*(x-r0)^3)/log(10)
    # where a = E Na / RT
    function blundy_wood(x,param)
        logD0, a, r0 = param
        return @. logD0 + 4π*a * (r0/2*(x-r0)^2 + 1/3*(x-r0)^3)/log(10)
    end

    #      log10(d0)   Ea/RT   r0 [pm]
    lower = [-10.0,    -1.0,      0]
    upper = [ 10.0,     0.0,    200]

    for j in 1:5
        charge = ("1+", "2+", "3+ree", "3+", "4+")[j]
        e = String.((e1, e2, ree3, e3, e4)[j])
        r = (e1r, e2r, ree3r, e3r, e4r)[j]
        e_bonus = String.((e1_bonus, e2_bonus, ree3_bonus, e3_bonus, e4_bonus)[j])
        r_bonus = (e1r_bonus, e2r_bonus, ree3r_bonus, e3r_bonus, e4r_bonus)[j]

        rplot = range(extrema(r)..., length=100)

        for m in (m1p, m2p, mree3p, m3p, m4p)[j]

            h = plot(title=m, 
                legend=:bottomleft, 
                ylabel="log10(Kd)",
                xlabel="Ionic radius [pm]",
                xticks=(r,e), 
                xrotation=-45, 
                framestyle=:box,)
            for i ∈ eachindex(pd["samples"])
                kD = fill(NaN,length(e))
                for j ∈ eachindex(e)
                    eⱼ = e[j]
                    kD[j] = haskey(pd[m], eⱼ) ? pd[m][eⱼ][i] : NaN
                end
                if count(!isnan, kD) > 0
                    plot!(h, r, kD, label="", seriestype=:scatter, color=lines[mod(i,length(lines))+1], msalpha=0)
                end

                kD_bonus = fill(NaN,length(e_bonus))
                for j ∈ eachindex(e_bonus)
                    eⱼ = e_bonus[j]
                    kD_bonus[j] = haskey(pd[m], eⱼ) ? pd[m][eⱼ][i] : NaN
                end
                if count(!isnan, kD_bonus) > 0
                    plot!(h, r_bonus, kD_bonus, label="", seriestype=:scatter, color=lines[mod(i,length(lines))+1], alpha=0.35, msalpha=0)
                end

                if (count(.~isnan.(kD)) > 2) && ((nanrange(r[.~isnan.(kD)]) > 0.55*nanrange(r)) || ((nanrange(r[.~isnan.(kD)]) > 0.45*nanrange(r)) && charge=="3+ree" ))
                    # Fit to Blundy and Wood curve
                    t = .~( isnan.(kD) .| isinf.(kD) )
                    param = [maximum(kD[t]), -1e-4, nanmean(r)] # initial guess
                    f = curve_fit(blundy_wood, r[t], kD[t], param; lower, upper) # Fit
                    @info "$m : $(f.param)"
                    if (f.param[2] < 0) || charge=="3+ree" # If there is no curvature, we probably haven't fit very well
                        plot!(h, rplot, blundy_wood(rplot,f.param), label="", color=lines[mod(i,length(lines))+1])
            
                        # Replace stored partition coefficients with new fits
                        for j in eachindex(e)
                            eⱼ = e[j]
                            if haskey(pd[m], eⱼ) && (isnan(pd[m][eⱼ][i]) || charge=="3+ree")
                                pd[m][eⱼ][i] = blundy_wood(r[j],f.param)
                            end
                        end
                    end
                end
            end
            savefig(h, "Calc_$(m)_$(charge).pdf")
        end
    end

    # Interpolate Eu partition coefficients where missing, assuming
    # 60% Eu as Eu2+ (c.f. Sr, Pb) and 40% as Eu3+ (c.f. Sm, Gd)
    for m in pd["minerals"]
        for i ∈ eachindex(pd["samples"])
            isnan(pd[m]["Eu2"][i]) && (pd[m]["Eu2"][i] = nanmean([pd[m]["Sr"][i], pd[m]["Pb"][i]]))
            isnan(pd[m]["Eu3"][i]) && (pd[m]["Eu3"][i] = nanmean([pd[m]["Sm"][i], pd[m]["Gd"][i]]))
            if isnan(pd[m]["Eu"][i])
                Eu_est = nanadd(0.6*exp10(pd[m]["Eu2"][i]), 0.4*exp10(pd[m]["Eu3"][i]))
                Eu_est > 0 && (pd[m]["Eu"][i] = log10(Eu_est))
            end
        end
    end

## --- Convert data to average as a function of SiO2

    # Convert from row-based to Si-based Dict
    kd = Dict{String, Union{Vector{String}, Vector{Float64}, Dict{String,Union{Float64, Vector{String}, Vector{Float64}}}}}()
    kd["minerals"] = pd["minerals"]
    for m in kd["minerals"]
        kd[m] = Dict{String,Union{Float64, Vector{String}, Vector{Float64}}}()
        kd[m]["elements"] = allelements
        for e in allelements
            t = .!isnan.(pd[m][e])
            if (count(t) > 2) && (nanrange(pd["SiO2"][t]) > 8)
                kd[m][e] = round.(mcfit(pd["SiO2"], pd["SiO2_sigma"], pd[m][e], pd[m][e*"_sigma"], 40, 80, 41, binwidth=5)[2], sigdigits=7)
            else
                kd[m][e] = round.(ones(41) * nanmean(pd[m][e]), sigdigits=7)
            end
            kd[m][e*"_sigma"] = round(nanstd(pd[m][e]),sigdigits=7)
        end
    end
    kd["SiO2"] = collect(40:80.)


## --- Save results

    # using MAT
    # matwrite(joinpath(path,"partitioncoeffs.mat"),p)
    f = open(joinpath(path,"PartitionCoefficients.jl"), "a")
    print(f, "\ngerm_kd = $kd\nexport germ_kd\n")
    close(f)

## --- End of File
