## --- Load (and install if neccesary) the StatGeochem package which has the resampling functions we'll want

    using StatGeochem
    using Plots; gr();

    if VERSION>=v"0.7"
        using Statistics, DelimitedFiles, SpecialFunctions
    end

## --- Download and unzip Keller and Schoene (2012) dataset

    if ~isfile("ign.h5") # Unless it already exists
        download("https://storage.googleapis.com/statgeochem/ign.h5.gz","./ign.h5.gz")
        download("https://storage.googleapis.com/statgeochem/err2srel.csv","./err2srel.csv")
        run(`gunzip -f ign.h5.gz`) # Unzip file
    end

    # Read HDF5 file
    using HDF5
    ign = h5read("ign.h5","vars")


## --- Compute proximity coefficients (inverse weights)

    # # Compute inverse weights
    # k = invweight(ign["Latitude"], ign["Longitude"], ign["Age"])

    # Since this is pretty computatually intensive, let's load a precomputed version instead
    k = ign["k"]

    # Probability of keeping a given data point when sampling
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the data in each resampling

    # Calculate age uncertainty
    ign["Age_sigma"] = (ign["Age_Max"]-ign["Age_Min"])/2;
    t = (ign["Age_sigma"] .< 50) .| isnan.(ign["Age_sigma"]) # Find points with < 50 Ma absolute uncertainty
    ign["Age_sigma"][t] .= 50 # Set 50 Ma minimum age uncertainty (1-sigma)

    # Set absolute uncertainties for each element where possible, using errors defined inerr2srel.csv
    ign["err2srel"] = importdataset("err2srel.csv", ',')
    for e in ign["elements"]
        # If there's an err2srel for this variable, create a "_sigma" if possible
        if haskey(ign["err2srel"], e) && !haskey(ign, e*"_sigma")
            ign[e*"_sigma"] = ign[e] .* (ign["err2srel"][e] / 2);
        end
    end

## --- Resample a single variable

    xmin = 0 # Minimum Age
    xmax = 1200 # Maximum Age
    nbins = 39
    elem = "K2O" # Element to plot

    # Look only at samples from a specific silica range
    t = 43 .< ign["SiO2"] .< 51 # Mafic
    # t = 51 .< ign["SiO2"] .< 62 # Intermediate
    # t = 62 .< ign["SiO2"] .< 74 # Felsic
    # t = trues(size(ign[elem])) # Everything

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_means(ign["Age"][t],ign[elem][t],xmin,xmax,nbins, p=p[t], x_sigma=ign["Age_sigma"][t])

    # Plot results
    plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkblue,markerstrokecolor=:auto,label="")
    plot!(xlabel="Age (Ma)", ylabel="$elem (wt. %)",xlims=(xmin,xmax),framestyle=:box,grid=:off,xflip=true) # Format plot


## ---  Resample a ratio

    tmin = 0 # Minimum age
    tmax = 3900 # Maximum age
    nbins = 39
    num = "La" # Numerator
    denom = "Yb" # Denominator

    # Look only at samples from a specific silica range
    t = 43 .< ign["SiO2"] .< 51 # Mafic
    # t = 51 .< ign["SiO2"] .< 62 # Intermediate
    # t = 62 .< ign["SiO2"] .< 74 # Felsic

    # Exclude outliers
    t = t .& inpctile(ign[num], 99) .& inpctile(ign[denom], 99)

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_ratios(ign["Age"][t],ign[num][t],ign[denom][t],tmin,tmax,nbins, p=p[t],
                    x_sigma=ign["Age_sigma"][t])

    # Plot results
    h = plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkred,markerstrokecolor=:auto,label="")
    plot!(h, xlabel="Age (Ma)", ylabel="$(num) / $(denom)",xlims=(tmin,tmax),framestyle=:box,grid=:off,xflip=true) # Format plot
    display(h)

    savefig(h,"$(num)$(denom)_$(tmax)-$(tmin)Ma.pdf")


## --- Single element differentiation example

    xelem = "SiO2"
    xmin = 45
    xmax = 75
    nbins = 8
    elem = "K2O"

    h = plot(xlabel=xelem, ylabel="$(elem)",xlims=(xmin,xmax),framestyle=:box,grid=:off,fg_color_legend=:white) # Format plot

    rt = [0,1,2,3,4] # Time range (Ga)
    colors = reverse(resize_colormap(viridis[1:end-20],length(rt)-1))
    for i=1:length(rt)-1
        t = rt[i]*1000 .< ign["Age"] .< rt[i+1]*1000

        # Resample, returning binned means and uncertainties
        # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
        (c,m,el,eu) = bin_bsr_means(ign[xelem][t],ign[elem][t],xmin,xmax,nbins, p=p[t],
                        x_sigma=ign[xelem*"_sigma"][t], y_sigma=ign[elem*"_sigma"][t])

        # Plot results
        plot!(h, c,m,yerror=(el,eu),color=colors[i],seriestype=:scatter,markerstrokecolor=:auto,label="$(rt[i])-$(rt[i+1]) Ga")
        plot!(h, c,m,style=:dot,color=colors[i],markerstrokecolor=:auto,label="")
    end
    # savefig(h,"$(xelem)_$(num)$(denom).pdf")
    display(h)

## --- Ratio differentiation example

    xelem = "SiO2"
    xmin = 45
    xmax = 75
    nbins = 8
    num = "Sc" # Numerator
    denom = "Yb" # Denominator

    h = plot(xlabel=xelem, ylabel="$(num) / $(denom)",xlims=(xmin,xmax),framestyle=:box,grid=:off,legend=:topleft,fg_color_legend=:white) # Format plot

    rt = [0,1,2,3,4]
    colors = reverse(resize_colormap(viridis[1:end-20],length(rt)-1))
    for i=1:length(rt)-1
        t = rt[i]*1000 .< ign["Age"] .< rt[i+1]*1000

        # Resample, returning binned means and uncertainties
        # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
        (c,m,el,eu) = bin_bsr_ratio_medians(ign[xelem][t],ign[num][t],ign[denom][t],xmin,xmax,nbins, p=p[t],
                        x_sigma=ign[xelem*"_sigma"][t], num_sigma=ign[num*"_sigma"][t], denom_sigma=ign[denom*"_sigma"][t])

        # Plot results
        plot!(h, c,m,yerror=(el,eu),color=colors[i],seriestype=:scatter,markerstrokecolor=:auto,label="$(rt[i])-$(rt[i+1]) Ga")
        plot!(h, c,m,style=:dot,color=colors[i],markerstrokecolor=:auto,label="")
    end
    display(h)

## --- Ratio differentiation

    xelem = "SiO2"
    xmin = 40 # Minimum age
    xmax = 80 # Maximum age
    nbins = 20
    num = "Sc" # Numerator
    denom = "Yb" # Denominator

    # Exclude outliers
    t = inpctile(ign[num], 99) .& inpctile(ign[denom], 99)

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_ratios(ign[xelem][t],ign[num][t],ign[denom][t],xmin,xmax,nbins, p=p[t],
                    x_sigma=ign[xelem][t]*0.01, num_sigma=ign[num][t]*0.05, denom_sigma=ign[denom][t]*0.05)

    # Plot results
    h = plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkred,markerstrokecolor=:auto,label="")
    plot!(h, xlabel=xelem, ylabel="$(num) / $(denom)",xlims=(xmin,xmax),framestyle=:box,grid=:off) # Format plot
    display(h)



## --- Export differentiation trends

    xelem = "SiO2"
    xmin = 45
    xmax = 75
    nbins = 10

    rt = [0,1,2,3,4] # Time range (Ga)
    data = Dict()
    for elem in ["Al2O3", "MgO", "Na2O", "Fe2O3T", "K2O", "CaO"]
        for i=1:length(rt)-1
            t = rt[i]*1000 .< ign["Age"] .< rt[i+1]*1000

            # Resample, returning binned means and uncertainties
            (c,m,e) = bin_bsr(ign[xelem][t],ign[elem][t],xmin,xmax,nbins, p=p[t],
                            x_sigma=ign[xelem*"_sigma"][t], y_sigma=ign[elem*"_sigma"][t])
            data[xelem] = c
            data[elem*"_$(rt[i])-$(rt[i+1])Ga"] = m
            data[elem*"_$(rt[i])-$(rt[i+1])Ga_sigma"] = e
        end
    end
    exportdataset(data,"MajorDifferentiation.csv",',')

## --- End of File
