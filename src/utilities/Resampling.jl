## --- Bootstrap resampling

    # Bootstrap resample (with uncertainty) a variable up to size nrows.
    # Optionally provide weights in p
    function bsresample(data::Array{<:Number}, sigma, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1
        pm = Progress(nrows, 1, "Resampling: ")
        while i <= nrows
            # Update progress
            update!(pm, i)

            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma,1) > 1
                    serr = sigma[t,:]
                else
                    serr = ones(size(sdata)) .* sigma
                end
            else # If only one sample
                sdata = data
                serr = sigma
            end

            # Randomize data over uncertainty interval
            sdata += randn(Float64, size(sdata)) .* serr

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end

        # Complete progress meter
        update!(pm, nrows)

        return resampled
    end

    # Second method for bsresample that takes a dictionary as input. Yay multiple dispatch!
    function bsresample(in::Dict, nrows, elements=in["elements"], p=min(0.2,nrows/length(in[elements[1]])))
        # 2d array of nominal values
        data = unelementify(in, elements, floatout=true)

        # 2d array of absolute 1-sigma uncertainties
        if haskey(in,"err") && isa(in["err"], Dict)
            sigma = unelementify(in["err"], elements, floatout=true)
        else
            sigma = unelementify(in, elements.*"_sigma", floatout=true)
        end

        # Resample
        sdata = bsresample(data, sigma, nrows, p)
        return elementify(sdata, elements, skipstart=0)
    end
    export bsresample

    # As bsresample, but with a uniform distribution stretching from data-sigma to data+sigma
    function bsresample_unif(data::Array{<:Number}, sigma, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1
        pm = Progress(nrows, 1, "Resampling: ")
        while i <= nrows
            # Update progress
            update!(pm, i)

            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma,1) > 1
                    serr = sigma[t,:]
                else
                    serr = ones(size(sdata)) .* sigma
                end
            else # If only one sample
                sdata = data
                serr = sigma
            end

            # Randomize data over uncertainty interval
            sdata += (2 .* rand(Float64, size(sdata)) .* serr) .- serr

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end

        # Complete progress meter
        update!(pm, nrows)

        return resampled
    end

    # Second method for bsresample_unif that takes a dictionary as input
    function bsresample_unif(in::Dict, nrows, elements=in["elements"], p=min(0.2,nrows/length(in[elements[1]])))
        # 2d array of nominal values
        data = unelementify(in, elements, floatout=true)

        # 2d array of absolute uncertainties
        if haskey(in,"err") && isa(in["err"], Dict)
            sigma = unelementify(in["err"], elements, floatout=true)
        else
            sigma = unelementify(in, elements.*"_sigma", floatout=true)
        end

        # Resample
        sdata = bsresample_unif(data, sigma, nrows, p)
        return elementify(sdata, elements, skipstart=0)
    end
    export bsresample_unif

    # As bsresample, but with a uniform distribution stretching from data-sigma to data+sigma, AND a gaussian component
    function bsresample_unif_norm(data::Array{<:Number}, sigma_unif, sigma_norm, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma_unif,1) > 1
                    serr_unif = sigma_unif[t,:]
                    serr_norm = sigma_norm[t,:]
                else
                    serr_unif = ones(size(sdata)) .* sigma_unif
                    serr_norm = ones(size(sdata)) .* sigma_norm
                end
            else # If only one sample
                sdata = data
                serr_unif = sigma_unif
                serr_norm = sigma_norm
            end

            # Randomize data over uncertainty interval
            sdata += (2 .* rand(Float64, size(sdata)) .* serr_unif) .- serr_unif # Uniform component
            sdata += randn(Float64, size(sdata)) .* serr_norm # Gaussian component

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end
        return resampled
    end
    export bsresample_unif_norm

    # As bsresample, but also return an index of the rows included from data
    function bsresample_index(data::Array{<:Number}, sigma, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))
        index = Array{Int}(undef,nrows)

        # Resample
        i = 1
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]
                sindex = findall(t)

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma,1) > 1
                    serr = sigma[t,:]
                else
                    serr = ones(size(sdata)) .* sigma
                end
            else # If only one sample
                sdata = data
                serr = sigma
                sindex = 1
            end

            # Randomize data over uncertainty interval
            sdata += randn(Float64, size(sdata)) .* serr

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
                index[i:i+size(sdata,1)-1] = sindex
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
                index[i:end] = sindex[1:nrows-i+1]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end
        return (resampled, index)
    end
    export bsresample_index

    # As bsresample_unif, but also return an index of the rows included from data
    function bsresample_unif_index(data::Array{<:Number}, sigma, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))
        index = Array{Int}(undef,nrows)

        # Resample
        i = 1
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]
                sindex = findall(t)

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma,1) > 1
                    serr = sigma[t,:]
                else
                    serr = ones(size(sdata)) .* sigma
                end
            else # If only one sample
                sdata = data
                serr = sigma
                sindex = 1
            end

            # Randomize data over uncertainty interval (uniform distribution)
            sdata += (2 .* rand(Float64, size(sdata)) .* serr) .- serr

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
                index[i:i+size(sdata,1)-1] = sindex
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
                index[i:end] = sindex[1:nrows-i+1]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end
        return (resampled, index)
    end
    export bsresample_unif_index

    # As bsresample, but with a uniform distribution stretching from age-sigma to age+sigma, AND a gaussian component
    function bsresample_unif_norm_index(data::Array{<:Number}, sigma_unif, sigma_norm, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))
        index = Array{Int}(undef,nrows)

        # Resample
        i = 1
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]
                sindex = findall(t)

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma_unif,1) > 1
                    serr_unif = sigma_unif[t,:]
                    serr_norm = sigma_norm[t,:]
                else
                    serr_unif = ones(size(sdata)) .* sigma_unif
                    serr_norm = ones(size(sdata)) .* sigma_norm
                end
            else # If only one sample
                sdata = data
                serr_unif = sigma_unif
                serr_norm = sigma_norm
                sindex = 1
            end

            # Randomize data over uncertainty interval
            sdata += (2 .* rand(Float64, size(sdata)) .* serr_unif) .- serr_unif # Uniform component
            sdata += randn(Float64, size(sdata)) .* serr_norm # Gaussian component

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
                index[i:i+size(sdata,1)-1] = sindex
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
                index[i:end] = sindex[1:nrows-i+1]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end
        return (resampled, index)
    end
    export bsresample_unif_norm_index


    # Bootstrap resample (without uncertainty) a variable to size nrows.
    # Optionally provide weights in p
    function randsample(data::Array{<:Number}, nrows::Number, p = min(0.2,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(Float64, size(data,1)) .< p
                sdata = data[t,:]
            else # If only one sample
                sdata = data
            end

            # Figure out how much of our resampled data to output
            if (i+size(sdata,1)-1) <= nrows
                resampled[i:i+size(sdata,1)-1,:] = sdata
            else
                resampled[i:end,:] = sdata[1:nrows-i+1,:]
            end

            # Keep track of current filled rows
            i += size(sdata,1)
        end
        return resampled
    end
    # Second method for randsample that takes a dictionary as input
    function randsample(in::Dict, nrows, elements=in["elements"], p=min(0.2,nrows/length(in[elements[1]])))
        data = unelementify(in, elements, floatout=true)
        sdata = randsample(data, nrows, p)
        return elementify(sdata, elements, skipstart=0)
    end
    export randsample


## --- Bin a dataset by a given independent variable

    function binmeans(x, y, min::Number, max::Number, nbins::Integer; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        means = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
            means[i] = mean(y[t])
            errors[i] = std(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
        end

        return (bincenters, means, errors)
    end
    function binmeans(x, y, min::Number, max::Number, nbins::Integer, tally::Bool; resamplingratio::Number=1)
        # As binmeans, but tally the number of samples (either resampled or corrected/original) that fall into bin
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        means = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        if tally
            N = Array{Int64}(undef,nbins)
            for i = 1:nbins
                t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
                means[i] = mean(y[t])
                errors[i] = std(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
                N[i] = count(t)
            end
        else
            N = Array{Float64}(undef,nbins)
            for i = 1:nbins
                t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
                means[i] = mean(y[t])
                errors[i] = std(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
                N[i] = count(t) / resamplingratio
            end
        end
        return (bincenters, means, errors, N)
    end
    export binmeans


    function binmedians(x, y, min::Number, max::Number, nbins::Integer; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        medians = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
            medians[i] = median(y[t])
            errors[i] = 1.4826 * nanmad(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
        end

        return (bincenters, medians, errors)
    end
    function binmedians(x, y, min::Number, max::Number, nbins::Ingeger, tally::Bool; resamplingratio::Number=1)
        # As binmedians, but tally the number of samples (either resampled or corrected/original) that fall into bin
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        medians = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        if tally
            N = Array{Int64}(undef,nbins)
            for i = 1:nbins
                t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
                medians[i] = median(y[t])
                errors[i] = 1.4826 * nanmad(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
                N[i] = count(t)
            end
        else
            N = Array{Float64}(undef,nbins)
            for i = 1:nbins
                t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
                medians[i] = median(y[t])
                errors[i] = 1.4826 * nanmad(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
                N[i] = count(t) / resamplingratio
            end
        end
        return (bincenters, medians, errors, N)
    end
    export binmedians

## --- Bin bootstrap resampled data

    function bin_bsr(x, y, min::Number, max::Number, nbins::Integer,x_sigma, nresamples::Integer, p::Number=0.2)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(x_sigma)))

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x),p)
            (c,m,s) = binmeans(dbs[:,1], dbs[:,2], min, max, nbins)
            means[:,i] = m
        end

        m = nanmean(means,dim=2)
        e = nanstd(means,dim=2)

        return (c, m, e)
    end
    export bin_bsr

    function bin_bsr_means(x, y, min::Number, max::Number, nbins::Integer, x_sigma, nresamples::Integer, p::Number=0.2)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(x_sigma)))

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x),p)
            (c,m,s) = binmeans(dbs[:,1], dbs[:,2], min, max, nbins)
            means[:,i] = m
        end

        m = nanmean(means,dim=2)
        el = m .- pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) .- m

        return (c, m, el, eu)
    end
    export bin_bsr_means

    function bin_bsr_medians(x, y, min::Number, max::Number, nbins::Integer, x_sigma, nresamples::Integer, p::Number=0.2)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(x_sigma)))

        medians = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x),p)
            (c,m,s) = binmedians(dbs[:,1], dbs[:,2], min, max, nbins)
            medians[:,i] = m
        end

        m = nanmedian(medians,dim=2)
        el = m .- pctile(medians,2.5,dim=2)
        eu = pctile(medians,97.5,dim=2) .- m

        return (c, m, el, eu)
    end
    export bin_bsr_medians

    function bin_bsr_ratios(x, num, denom, min::Number, max::Number, nbins::Integer, x_sigma, num_sigma, denom_sigma, nresamples::Integer, p::Number=0.2)
        data = hcat(x, num, denom)
        sigma = hcat(x_sigma, num_sigma, denom_sigma)

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x),p)
            (c,m,s) = binmeans(dbs[:,1], dbs[:,2] ./ (dbs[:,2] .+ dbs[:,3]), min, max, nbins)
            means[:,i] = m ./ (1 .- m)
        end

        m = nanmean(means,dim=2)
        el = m .- pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) .- m

        return (c, m, el, eu)
    end
    export bin_bsr_ratios


## --- Downsample an image / array

    function downsample(matrix::Array,factor::Int,jfactor=factor::Int)
        if ndims(matrix)==2
            rows = floor(Int,size(matrix,1)/factor)
            cols = floor(Int,size(matrix,2)/jfactor)

            downsampled = typeof(matrix)(undef,rows,cols)
            for i=1:rows
                for j=1:cols
                    downsampled[i,j]=matrix[i*factor,j*jfactor]
                end
            end
        else
            downsampled = matrix[factor:factor:end]
        end
        return downsampled
    end
    export downsample

## --- Spatiotemporal sample weighting

    # Produce a weighting coefficient for each row of data corresponding
    # to the input lat, lon, and age that is inversely proportional to the
    # spatiotemporal data concentration
    function invweight(lat::Array{<:Number}, lon::Array{<:Number}, age::Array{<:Number};
        lp::Number=2, spatialscale::Number=1.8, agescale::Number=38.0)

        # Check if there is lat, lon, and age data
        nodata = isnan.(lat) .| isnan.(lon) .| isnan.(age)

        k = Array{Float64}(undef,length(lat))
        @showprogress 1 "Calculating weights: " for i=1:length(lat)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( 1.0 ./ ( (arcdistance(lat[i],lon[i],lat,lon)./spatialscale).^lp .+ 1.0) .+ 1.0./((abs.(age[i] .- age)./agescale).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight

    # Produce a weighting coefficient for each row of data corresponding
    # to the input lat, lon, and age that is inversely proportional to the
    # spatial data concentration
    function invweight_location(lat::Array{<:Number}, lon::Array{<:Number};
        lp::Number=2, spatialscale::Number=1.8)

        # Check if there is lat, lon, and age data
        nodata = isnan.(lat) .| isnan.(lon)

        k = Array{Float64}(undef,length(lat))
        @showprogress 1 "Calculating weights: " for i=1:length(lat)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( 1.0 ./ ( (arcdistance(lat[i],lon[i],lat,lon)/spatialscale).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight_location

    # Produce a weighting coefficient for each row of data corresponding
    # to the input age that is inversely proportional to the
    # temporal data concentration
    function invweight_age(age::Array{<:Number}; lp::Number=2, agescale::Number=38.0)

        # Check if there is lat, lon, and age data
        nodata = isnan.(age)

        k = Array{Float64}(undef,length(age))
        @showprogress 1 "Calculating weights: " for i=1:length(age)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum(1.0 ./ ( (abs.(age[i] .- age)./agescale).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight_age

## --- End of file
