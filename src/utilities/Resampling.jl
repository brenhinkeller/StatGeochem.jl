## --- Bootstrap resampling

    # Bootstrap resample (with uncertainty) a variable up to size nrows.
    # Optionally provide weights in p
    function bsresample(data::Array{<:Number}, sigma, nrows::Number, p = min(0.5,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1;
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(size(data,1)) .< p
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
            sdata += randn(size(sdata)) .* serr

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

    # Second method for bsresample that takes a dictionary as input
    function bsresample(in::Dict, nrows, elements=in["elements"], p=min(0.5,nrows/length(in[elements[1]])))
        data = unelementify(in, elements, floatout=true)
        sigma = unelementify(in, elements.*"_sigma", floatout=true)
        sdata = bsresample(data, sigma, nrows, p)
        return elementify(sdata, elements)
    end
    export bsresample

    # Bootstrap resample (without uncertainty) a variable to size nrows.
    # Optionally provide weights in p
    function randsample(data::Array{<:Number}, nrows::Number, p = min(0.5,nrows/size(data,1)))
        # Allocate output array
        resampled = Array{Float64}(undef,nrows,size(data,2))

        # Resample
        i = 1;
        while i <= nrows
            # If we have more than one sample
            if size(data,1) > 1
                # Select weighted sample of data
                t = rand(size(data,1)) .< p
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
    function randsample(in::Dict, nrows, elements=in["elements"], p=min(0.5,nrows/length(in[elements[1]])))
        data = unelementify(in, elements, floatout=true)
        sdata = randsample(data, nrows, p)
        return elementify(sdata, elements)
    end
    export randsample
## --- Bin a dataset by a given independent variable

    function binmeans(x,y,min,max,nbins; resamplingratio=1)
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        averages = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
            averages[i] = mean(y[t])
            errors[i] = std(y[t]) ./ sqrt(sum(t)) .* sqrt(resamplingratio)
        end

        return(bincenters, averages, errors)
    end
    export binmeans

## --- Bin bootstrap resampled data

    function bin_bsr_means(x,y,min,max,nbins,x_sigma,nresamples,p=0.5)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(x_sigma)))

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x),p)
            (c,m,s) = binmeans(dbs[:,1], dbs[:,2], min, max, nbins)
            means[:,i] = m;
        end

        m = nanmean(means,dim=2)
        el = m .- pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) .- m

        return (c, m, el, eu)
    end
    export bin_bsr_means

## --- Downsample an image / array

    function downsample(matrix::Array,factor::Int,jfactor=factor::Int)
        if ndims(matrix)==2
            rows = floor(Int,size(matrix,1)/factor)
            cols = floor(Int,size(matrix,2)/jfactor)

            downsampled = typeof(matrix)(rows,cols)
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
    function invweight(lat::Array{<:Number}, lon::Array{<:Number}, age::Array{<:Number}; lp=2)

        # Check if there is lat, lon, and age data
        nodata = isnan.(lat) .| isnan.(lon) .| isnan.(age)

        k = Array{Float64}(undef,length(lat))
        for i=1:length(lat)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( 1.0 ./ ( (arcdistance(lat[i],lon[i],lat,lon)./1.8).^lp .+ 1.0) .+ 1.0./((abs.(age[i] .- age)./38.0).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight

    # Produce a weighting coefficient for each row of data corresponding
    # to the input lat, lon, and age that is inversely proportional to the
    # spatial data concentration
    function invweight_location(lat::Array{<:Number}, lon::Array{<:Number}; lp=2)

        # Check if there is lat, lon, and age data
        nodata = isnan.(lat) .| isnan.(lon)

        k = Array{Float64}(undef,length(lat))
        for i=1:length(lat)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( 1.0 ./ ( (arcdistance(lat[i],lon[i],lat,lon)/1.8).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight_location

    # Produce a weighting coefficient for each row of data corresponding
    # to the input age that is inversely proportional to the
    # temporal data concentration
    function invweight_age(age::Array{<:Number}; lp=2)

        # Check if there is lat, lon, and age data
        nodata = isnan.(age)

        k = Array{Float64}(undef,length(lat))
        for i=1:length(lat)
            if nodata[i] # If there is no data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum(1.0 ./ ( (abs.(age[i] .- age)./38.0).^lp .+ 1.0) )
            end
        end
        return k
    end
    export invweight_age

## --- End of file
