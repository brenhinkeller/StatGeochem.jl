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
            pr[name] = importdataset(joinpath(path, name*".tsv"), '\t', mindefinedcolumns=2)
        end
    end

    # Approximate content of various rock types, for when actual melt SiO2 not reported
    str = ["basalt","andesit","dacit","rhyolit","hawaiite","alkali basalt","basanit","benmoreite","camptonite","eclogit","garnet pyroxenit","kimberlit","komatiit","lampro","latite","leucitit","lherzolit","mugearite","pantellerite","trachyt","peridotit","phonolit","picrit","syen","morb","tholeiit","tonalit","granodiorit","granit","shoshonit","carbonatit","aplite","monzonit","synthetic"]
    val = [49,57.5,67.5,73,48.5,47,44,57,42,49,49,40,45,46.5,62.5,47,42,52,73,61,42,59,45.5,61,49,50,65,67.5,72,55,2,73,67.4,60]
    err = [5,6.5,5.5,5,3.5,4,3,4,3,5,5,4,5,11.5,11.5,5,5,5,5,9,5,5,5.5,9,5,4,8,7.5,5,5,5,5,5,20]

    # Determine the unique reference / rock type pairs
    samples = Dict{String,Vector{Tuple{String,String,Float64,Float64}}}()
    allsamples = Tuple{String,String,Float64,Float64}[]
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
        for i=1:length(str)
            t = containsi.(pr[m]["Rock Type"], str[i]) .& needssi
            # Average together, such that "basaltic andesite" is halfway between "basalt" and "andesite"
            pr[m]["SiO2"][t] = nanmean([pr[m]["SiO2"][t] val[i]*ones(count(t))], dim=2)
            pr[m]["SiO2_sigma"][t] = nanmean([pr[m]["SiO2_sigma"][t] err[i]*ones(count(t))], dim=2)
        end

        # Save unique samples
        println(m)
        flush(stdout)
        samples[m] = collect(zip(string.(pr[m]["Reference"]), string.(pr[m]["Rock Type"]), pr[m]["SiO2"], pr[m]["SiO2_sigma"]))
        append!(allsamples, unique(samples[m]))

        # Save unique elements
        append!(allelements, unique(pr[m]["Elem"]))
    end
    allsamples = sort(unique(allsamples))
    allelements = sort(unique(allelements))


    # Second dict: sort the data by sample ID
    pd = Dict{String, Union{Vector, Dict{String, Union{Vector{Float64}, Vector{String}}}}}()
    pd["samples"] = allsamples
    pd["Reference"] = allsamples .|> x -> x[1]
    pd["Rock Type"] = allsamples .|> x -> x[2]
    pd["SiO2"] = allsamples .|> x -> x[3]
    pd["SiO2_sigma"] = allsamples .|> x -> x[4]
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

## --- Fit REEs as a function of activation energy and bulk modulus

    ree3 = ["La","Pr","Nd","Sm","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",] # 3+ rare earth elements
    r = [0.106,0.101,0.100,0.096,0.094,0.092,0.091,0.089,0.088,0.087,0.086,0.085,] # Atomic Radii
    rp = 0.085:0.001:0.107 # Radius range to plot over

    # Equation we're fitting to (from Blundy and Wood, 1994):
    # lnD0 + a * (r0/2*(x-r0)^2 + 1/3*(x-r0)^3)
    # where a = 4π E Na / RT
    @. blundy_wood(x,param) = param[1] + param[2] * (param[3]/2*(x-param[3])^2 + 1/3*(x-param[3])^3)

    # for mineral={'Apatite','Amphibole','Clinopyroxene','Orthopyroxene','Garnet','Sphene','Allanite','Baddeleyite'}
    for m in pd["minerals"]
        h = plot(title=m, legend=:bottomleft)
        for i = 1:length(pd["samples"])
            kD = fill(NaN,length(ree3))
            for j = 1:length(ree3)
                kD[j] = pd[m][ree3[j]][i]
            end

            if count(.~isnan.(kD)) > 1
                # plot!(h, r, kD, label="$i", seriestype=:scatter, color=lines[mod(i,length(lines))+1], msalpha=0)
                plot!(h, r, kD, label="", seriestype=:scatter, color=lines[mod(i,length(lines))+1], msalpha=0)
            end

            if (count(.~isnan.(kD)) > 3) && (nanrange(r[.~isnan.(kD)]) > 0.013)
                # Fit to Blundy and Wood curve
                t = .~( isnan.(kD) .| isinf.(kD) )
                param = [maximum(kD[t]), -10000, 0.095] # initial guess
                f = curve_fit(blundy_wood, r[t], kD[t], param) # Fit
                plot!(h, r, blundy_wood(r,f.param), label="", color=lines[mod(i,length(lines))+1]) # Plot

                # Replace stored partition coefficients with new fits
                for j = 1:length(ree3)
                    pd[m][ree3[j]][i] = blundy_wood(r[j],f.param)
                end
            end
        end
        savefig(h, "Calc_$(m)_REE.pdf")
    end

    # Interpolate Eu partition coefficients where missing, assuming
    # 60% Eu as Eu2+ (c.f. Ba, Sr, Ca) and 40% as Eu3+ (c.f. Sm, Gd)
    # for m in = ["Albite","Anorthite","Orthoclase","Apatite"]
    for m in pd["minerals"]
        for i = 1:length(pd["samples"])
            if isnan(pd[m]["Eu"][i])
                pd[m]["Eu"][i] = log10(0.6*10^nanmean([pd[m]["Ba"][i], pd[m]["Sr"][i]]) + 0.4*10^nanmean([pd[m]["Sm"][i], pd[m]["Gd"][i]]))
            end
        end
    end
    for m in ("Monazite", "Xenotime", "Allanite")
        for i = 1:length(pd["samples"])
            if isnan(pd[m]["Eu"][i])
                pd[m]["Eu"][i] = log10(0.6*0 + 0.4*10^nanmean([pd[m]["Sm"][i], pd[m]["Gd"][i]]))
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
                kd[m][e] = mcfit(pd["SiO2"], pd["SiO2_sigma"], pd[m][e], pd[m][e*"_sigma"], 40, 80, 41, binwidth=5)[2]
            else
                kd[m][e] = ones(41) * nanmean(pd[m][e])
            end
            kd[m][e*"_sigma"] = nanstd(pd[m][e])
        end
    end
    kd["SiO2"] = collect(40:80.)

    # Set Albite partiton coefficients
    for e in kd["Albite"]["elements"]
        kd["Albite"][e] = nanmean([kd["Albite"][e] kd["Orthoclase"][e] kd["Anorthite"][e]], dim=2)
    end
    kd["note"] = ["kd for Albite is nanmean of AlkaliFeldspar, Orthoclase, and Anorthite",]

## --- Save results

    # using MAT
    # matwrite(joinpath(path,"partitioncoeffs.mat"),p)
    f = open(joinpath(path,"PartitionCoefficients.jl"), "a")
    print(f, "\ngerm_kd = $kd\nexport germ_kd\n")
    close(f)

## --- End of File
