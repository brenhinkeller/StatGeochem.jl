## --- Load (and install if neccesary) the StatGeochem package which has the resampling functions we'll want

    using StatGeochem
    using Plots; gr();

    if VERSION>=v"0.7"
        using Statistics, DelimitedFiles, SpecialFunctions
    end

## --- Download and unzip Keller and Schoene (2012) dataset

    if ~isfile("ign.h5") # Unless it already exists
        download("https://storage.googleapis.com/statgeochem/ign.h5.gz","./ign.h5.gz")
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


## --- Resample a single variable

    nresamplings=1000
    xmin = 0 # Minimum Age
    xmax = 1200 # Maximum Age
    nbins = 39
    elem = "K2O" # Element to plot

    # Look only at samples from a specific silica range
    t = (ign["SiO2"].>43) .& (ign["SiO2"].<51) # Mafic
    # t = (ign["SiO2"].>51) .& (ign["SiO2"].<62) # Intermediate
    # t = (ign["SiO2"].>62) .& (ign["SiO2"].<74) # Felsic
    # t = trues(size(ign[elem])) # Everything

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_means(ign["Age"][t],ign[elem][t],xmin,xmax,nbins,ign["Age_sigma"][t],nresamplings,p[t])

    # Plot results
    plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkblue,markerstrokecolor=:auto,label="")
    plot!(xlabel="Age (Ma)", ylabel="$elem (wt. %)",xlims=(xmin,xmax),framestyle=:box,grid=:off,xflip=true) # Format plot


## ---  Resample a ratio

    nresamplings=1000
    tmin = 0 # Minimum age
    tmax = 3900 # Maximum age
    nbins = 39
    num = "La" # Numerator
    denom = "Yb" # Denominator

    # Look only at samples from a specific silica range
    t = (ign["SiO2"].>43) .& (ign["SiO2"].<51) # Mafic
    # t = (ign["SiO2"].>51) .& (ign["SiO2"].<62) # Intermediate
    # t = (ign["SiO2"].>62) .& (ign["SiO2"].<74) # Felsic
    # t = trues(size(ign[elem])) # Everything

   # Exclude outliers
    t = t .& (ign[num] .> pctile(ign[num],0.5)) .& (ign[num] .< pctile(ign[num],99.5))
    t = t .& (ign[denom] .> pctile(ign[denom],0.5)) .& (ign[denom] .< pctile(ign[denom],99.5))

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_ratios(ign["Age"][t],ign[num][t],ign[denom][t],tmin,tmax,nbins,ign["Age_sigma"][t],ign[num][t]*0.05,ign[denom][t]*0.05,nresamplings,p[t])

    # Plot results
    h = plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkred,markerstrokecolor=:auto,label="")
    plot!(h, xlabel="Age (Ma)", ylabel="$(num) / $(denom)",xlims=(tmin,tmax),framestyle=:box,grid=:off,xflip=true) # Format plot
    display(h)

    savefig(h,"$(num)$(denom)_$(tmax)-$(tmin)Ma.pdf")

## --- Ratio differentiation

    nresamplings=1000
    xelem = "SiO2"
    xmin = 40 # Minimum age
    xmax = 80 # Maximum age
    nbins = 20
    num = "Sc" # Numerator
    denom = "Yb" # Denominator

    # Look only at samples from a specific silica range
    t = (ign["SiO2"].>40) .& (ign["SiO2"].<80) # Mafic

   # Exclude outliers
    t = t .& (ign[num] .> pctile(ign[num],0.5)) .& (ign[num] .< pctile(ign[num],99.5))
    t = t .& (ign[denom] .> pctile(ign[denom],0.5)) .& (ign[denom] .< pctile(ign[denom],99.5))

    # Resample, returning binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    (c,m,el,eu) = bin_bsr_ratios(ign[xelem][t],ign[num][t],ign[denom][t],xmin,xmax,nbins,ign[xelem][t]*0.01,ign[num][t]*0.05,ign[denom][t]*0.05,nresamplings,p[t])

    # Plot results
    h = plot(c,m,yerror=(el,eu),seriestype=:scatter,color=:darkred,markerstrokecolor=:auto,label="")
    plot!(h, xlabel=xelem, ylabel="$(num) / $(denom)",xlims=(xmin,xmax),framestyle=:box,grid=:off) # Format plot
    display(h)

    savefig(h,"$(xelem)_$(num)$(denom).pdf")

## --- End of File
