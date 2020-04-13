## --- Bootstrap resampling

    """
    ```julia
    resampled = bsr(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
    ```

    Return data boostrap resampled from of a (sample-per-row / element-per-column)
    dataset `data` with uncertainties `sigma` and variable resampling probabilities `p`
    """
    function bsr(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
        # Allocate output array
        nrows_initial = size(data,1)
        ncolumns = size(data,2)
        r = randn(Float64,nrows,ncolumns) # All the normal RVs we'll need
        resampled = Array{Float64}(undef,nrows,ncolumns) # Allocate output array

        # Resample
        n = 1
        accepted = Array{Int}(undef, nrows_initial)
        while n <= nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:nrows_initial
                if rand(Float64) < p[i]
                    nrows_accepted += 1
                    accepted[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n + 1)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                k = n
                # Optimized inner loop
                @inbounds @simd for i = 1:nrows_new
                    a = accepted[i]
                    resampled[k,j] = data[a,j] + r[k,j] * sigma[a,j]
                    k += 1
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    """
    ```julia
    resampled = bsr(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
    ```

    Return data boostrap resampled from of a (sample-per-row / element-per-column)
    dataset `data` with uncertainties `sigma` and uniform resampling probability `p`
    """
    function bsr(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::Number = min(0.36787944,nrows/size(data,1)))
        # Allocate output array
        nrows_initial = size(data,1)
        ncolumns = size(data,2)
        r = randn(Float64,nrows,ncolumns) # All the normal RVs we'll need
        resampled = Array{Float64}(undef,nrows,ncolumns) # Allocate output array

        # Resample
        n = 1
        accepted = Array{Int}(undef, nrows_initial)
        while n <= nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:nrows_initial
                if rand(Float64) < p
                    nrows_accepted += 1
                    accepted[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n + 1)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                k = n
                # Optimized inner loop
                @inbounds @simd for i = 1:nrows_new
                    a = accepted[i]
                    resampled[k,j] = data[a,j] + r[k,j] * sigma[a,j]
                    k += 1
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    export bsr


    """
    ```julia
    bsr!(resampled::AbstractArray{<:Number}, data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
    ```

    Fill `resampled` with data boostrap resampled from a (sample-per-row / element-per-column)
    dataset `data` with uncertainties `sigma` and variable resampling probabilities `p`
    """
    function bsr!(resampled::AbstractArray{<:Number}, data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
        # Prepare
        nrows_initial = size(data,1)
        ncolumns = size(data,2)
        r = randn(Float64,nrows,ncolumns) # All the normal RVs we'll need
        accepted = Array{Int}(undef, nrows_initial) # Array to hold list of accepted samples

        # Resample
        n = 1
        while n <= nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:nrows_initial
                if rand(Float64) < p[i]
                    nrows_accepted += 1
                    accepted[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n + 1)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                k = n
                # Optimized inner loop
                @inbounds @simd for i = 1:nrows_new
                    a = accepted[i]
                    resampled[k,j] = data[a,j] + r[k,j] * sigma[a,j]
                    k += 1
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    """
    ```julia
    bsr!(resampled::AbstractArray{<:Number}, data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::AbstractArray{<:Number})
    ```

    Fill `resampled` with data boostrap resampled from a (sample-per-row / element-per-column)
    dataset `data` with uncertainties `sigma` and uniform resampling probability `p`
    """
    function bsr!(resampled::AbstractArray{<:Number}, data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}, nrows::Integer, p::Number = min(0.36787944,nrows/size(data,1)))
        # Prepare
        nrows_initial = size(data,1)
        ncolumns = size(data,2)
        r = randn(Float64,nrows,ncolumns) # All the normal RVs we'll need
        accepted = Array{Int}(undef, nrows_initial) # Array to hold list of accepted samples

        # Resample
        n = 1
        while n <= nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:nrows_initial
                if rand(Float64) < p
                    nrows_accepted += 1
                    accepted[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n + 1)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                k = n
                # Optimized inner loop
                @inbounds @simd for i = 1:nrows_new
                    a = accepted[i]
                    resampled[k,j] = data[a,j] + r[k,j] * sigma[a,j]
                    k += 1
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    export bsr!

    # Bootstrap resample (with uncertainty) a variable up to size nrows.
    # Optionally provide weights in p
    function bsresample(data::Array{<:Number}, sigma::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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

        return resampled
    end
    # Second method for bsresample that takes a dictionary as input. Yay multiple dispatch!
    function bsresample(in::Dict, nrows::Integer, elements=in["elements"],
        p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/length(in[elements[1]])))

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
    function bsresample_unif(data::Array{<:Number}, sigma::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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

        return resampled
    end
    # Second method for bsresample_unif that takes a dictionary as input
    function bsresample_unif(in::Dict, nrows::Integer, elements=in["elements"],
        p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/length(in[elements[1]])))

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
    function bsresample_unif_norm(data::Array{<:Number},
        sigma_unif::Union{Number,Array{<:Number}}, sigma_norm::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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
    function bsresample_index(data::Array{<:Number}, sigma::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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
    function bsresample_unif_index(data::Array{<:Number}, sigma::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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
    function bsresample_unif_norm_index(data::Array{<:Number},
        sigma_unif::Union{Number,Array{<:Number}}, sigma_norm::Union{Number,Array{<:Number}},
        nrows::Integer, p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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
    function randsample(data::Array{<:Number}, nrows::Integer,
        p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/size(data,1)))

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
    function randsample(in::Dict, nrows::Integer, elements=in["elements"],
        p::Union{Number,AbstractVector{<:Number}} = min(0.2,nrows/length(in[elements[1]])))

        data = unelementify(in, elements, floatout=true)
        sdata = randsample(data, nrows, p)

        return elementify(sdata, elements, skipstart=0)
    end
    export randsample


## --- Bin a dataset by a given independent variable

    """
    ```julia
    (bincenters, N) = bincounts(x::AbstractArray{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
    ```

    Tally the number of samples that fall into each of `nbins` equally spaced
    bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function bincounts(x::AbstractArray{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Tally the number of samples (either resampled or corrected/original) that fall into each bin
        binwidth = (xmax-xmin)/nbins
        bincenters = (xmin+binwidth/2):binwidth:(xmax-binwidth/2)

        # Calculate index from x value
        scalefactor = nbins / (xmax - xmin)
        index = (x .- xmin) .* scalefactor

        # Add up the results
        N = fill(0,nbins)
        @inbounds for i in index
            if 0 < i < nbins
                N[ceil(Int, i)] += 1
            end
        end
        return bincenters, N
    end
    export bincounts

    # The nanmean of y binned by x, returning bincenters, means, and standard error of the mean
    function binmeans(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}, xmin::Number, xmax::Number, nbins::Integer; resamplingratio::Number=1)
        binwidth = (xmax-xmin)/nbins
        bincenters = (xmin+binwidth/2):binwidth:(xmax-binwidth/2)

        # Calculate index from x value
        scalefactor = nbins / (xmax - xmin)
        index_float = (x .- xmin) .* scalefactor

        # Calculate the nanmeans and nansems
        N = fill(0,nbins)
        mu = fill(0.0,nbins)
        s2 = fill(0.0,nbins)
        for i = 1:length(x)
            if (0 < index_float[i] < nbins) && !isnan(y[i])
                index = ceil(Int, index_float[i])
                N[index] += 1
                delta = y[i] - mu[index]
                mu[index] += delta / N[index]
                s2[index] += delta * (y[i] - mu[index])
            end
        end

        # Calculate standard errors
        se = Array{Float64}(undef, nbins)
        for i=1:nbins
            if N[i] > 0
                s2[i] /= N[i] - 1 # Variance
                se[i] = sqrt(s2[i] * resamplingratio / N[i])
            else
                mu[i] = NaN
                se[i] = NaN
            end
        end
        return bincenters, mu, se
    end
    function binmeans(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}, min::Number, max::Number, nbins::Integer, weight::AbstractArray{<:Number}; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        means = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<=binedges[i+1]) .& (.~isnan.(y))
            w = ProbabilityWeights(weight[t])
            means[i] = mean(y[t], w)
            errors[i] = std(y[t], w, corrected=true) * sqrt(resamplingratio) / sqrt(count(t))
        end

        return bincenters, means, errors
    end
    export binmeans

    # The nanmedian of y binned by x, retunring bincenters, medians, and equivalent standard error of the mean (1.4828 * median abolute deviation)
    function binmedians(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}, min::Number, max::Number, nbins::Integer; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = linsp(min,max,nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        medians = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<=binedges[i+1]) .& (.~isnan.(y))
            medians[i] = median(y[t])
            errors[i] = 1.4826 * nanmad(y[t]) * sqrt(resamplingratio / count(t))
        end

        return bincenters, medians, errors
    end
    export binmedians

## --- Bin bootstrap resampled data

    """
    ```julia
    (c, m, e) = bin_bsr(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer,
        \tx_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)
    ```

    Returns the bincenters `c`, means `m`, and 1σ standard errors of the mean `e` for
    a variable `y` binned by variable `x` into `nbins` equal bins between `xmin` and `xmax`,
    after `nresamples` boostrap resamplings with acceptance probability `p`.
    """
    function bin_bsr(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(y)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                bsr!(dbs,data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            nanmean!(view(means,:,i), view(dbs,ns,1), view(dbs,ns,2), xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        e = nanstd(means,dim=2) # Standard deviation of means (sem)

        return c, m, e
    end
    """
    ```julia
    (c, m, e) = bin_bsr(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer,
        \t x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}, w::AbstractVector{<:Number})
    ```

    Returns the bincenters `c`, means `m`, and 1σ standard errors of the mean `e` for
    a variable `y` binned by variable `x` into `nbins` equal bins between `xmin` and `xmax`,
    after `nresamples` boostrap resamplings with acceptance probability `p` and weight `w`.
    """
    function bin_bsr(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}, w::AbstractVector{<:Number})
        data = hcat(x, y, w)
        sigma = hcat(x_sigma, zeros(size(y)), zeros(size(w)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                bsr!(dbs,data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            nanmean!(view(means,:,i), view(dbs,ns,1), view(dbs,ns,2), view(dbs,ns,3), xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        e = nanstd(means,dim=2) # Standard deviation of means (sem)

        return c, m, e
    end
    """
    ```julia
    (c, M, E) = bin_bsr(x::AbstractVector{<:Number}, Y::AbstractMatrix{<:Number}, xmin::Number, xmax::Number, nbins::Integer,
        \tx_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2;
        \ty_sigma::AbstractMatrix{<:Number}=zeros(size(Y)))
    ```

    Returns the bincenters `c`, means `M`, and 1σ standard errors of the mean `E` for
    every (column-wise) variable in a 2-d array `Y`, binned by variable `x` into `nbins` equal bins between `xmin` and `xmax`,
    after `nresamples` boostrap resamplings with acceptance probability `p`.
    """
    function bin_bsr(x::AbstractVector{<:Number}, Y::AbstractMatrix{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2; Y_sigma::AbstractMatrix{<:Number}=zeros(size(Y)))
        data = hcat(x, Y)
        sigma = hcat(x_sigma, Y_sigma)
        dtype = float(eltype(data))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{dtype}(undef, nrows*chunk, ncols)
        means = Array{dtype}(undef, nbins, nresamples, size(Y,2))
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                bsr!(dbs,data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            nanmean!(view(means,:,i,:), view(dbs,ns,1), view(dbs,ns,2:1+size(Y,2)), xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = Array{dtype}(undef, nbins, size(Y,2))
        e = Array{dtype}(undef, nbins, size(Y,2))
        for j = 1:size(Y,2)
            m[:,j] .= nanmean(view(means,:,:,j),dim=2) # Mean-of-means
            e[:,j] .= nanstd(view(means,:,:,j),dim=2) # Standard deviation of means (sem)
        end
        return c, m, e
    end
    export bin_bsr


    function bin_bsr_means(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)

        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(y)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                bsr!(dbs,data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            nanmean!(view(means,:,i), view(dbs,ns,1), view(dbs,ns,2), xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    function bin_bsr_means(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}, w::AbstractVector{<:Number})

        data = hcat(x, y, w)
        sigma = hcat(x_sigma, zeros(size(y)), zeros(size(w)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                bsr!(dbs,data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            nanmean!(view(means,:,i), view(dbs,ns,1), view(dbs,ns,2), view(dbs,ns,3), xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    export bin_bsr_means

    function bin_bsr_medians(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer, x_sigma::AbstractVector{<:Number}, nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)

        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(y)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        medians = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                dbs .= bsr(data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            medians[:,i] = nanmedian(dbs[ns,1], dbs[ns,2], xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmedian(medians,dim=2) # Median-of-medians
        el = m .- pctile(medians,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(medians,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    export bin_bsr_medians

    function bin_bsr_ratios(x::AbstractVector{<:Number}, num::AbstractVector{<:Number}, denom::AbstractVector{<:Number},
        xmin::Number, xmax::Number, nbins::Integer,
        x_sigma::AbstractVector{<:Number}, num_sigma::AbstractVector{<:Number}, denom_sigma::AbstractVector{<:Number},
        nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)

        data = hcat(x, num, denom)
        sigma = hcat(x_sigma, num_sigma, denom_sigma)
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                dbs .= bsr(data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            means[:,i] = nanmean(dbs[ns,1], dbs[ns,2], dbs[ns,3], xmin, xmax, nbins)
            m = nanmean(dbs[ns,1], dbs[ns,2] ./ (dbs[ns,2] .+ dbs[ns,3]), xmin, xmax, nbins)
            means[:,i] = m ./ (1 .- m)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    function bin_bsr_ratios(x::AbstractVector{<:Number}, num::AbstractVector{<:Number}, denom::AbstractVector{<:Number},
        xmin::Number, xmax::Number, nbins::Integer,
        x_sigma::AbstractVector{<:Number}, num_sigma::AbstractVector{<:Number}, denom_sigma::AbstractVector{<:Number},
        nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}, w::AbstractVector{<:Number})

        data = hcat(x, num, denom, w)
        sigma = hcat(x_sigma, num_sigma, denom_sigma, zeros(size(w)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        means = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                dbs .= bsr(data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            means[:,i] = nanmean(dbs[ns,1], dbs[ns,2], dbs[ns,3], xmin, xmax, nbins)
            m = nanmean(dbs[ns,1], dbs[ns,2] ./ (dbs[ns,2] .+ dbs[ns,3]), dbs[ns,4], xmin, xmax, nbins)
            means[:,i] = m ./ (1 .- m)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    export bin_bsr_ratios


    function bin_bsr_ratio_medians(x::AbstractVector{<:Number}, num::AbstractVector{<:Number}, denom::AbstractVector{<:Number},
        xmin::Number, xmax::Number, nbins::Integer,
        x_sigma::AbstractVector{<:Number}, num_sigma::AbstractVector{<:Number}, denom_sigma::AbstractVector{<:Number},
        nresamples::Integer, p::Union{Number,AbstractVector{<:Number}}=0.2)

        data = hcat(x, num, denom)
        sigma = hcat(x_sigma, num_sigma, denom_sigma)
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        chunk = ceil(Int,10^6/nrows)
        dbs = Array{Float64}(undef, nrows*chunk, ncols)
        medians = Array{Float64}(undef, nbins, nresamples)
        # Resample
        for i=1:nresamples
            k = mod(i-1,chunk)
            if k == 0
                dbs .= bsr(data,sigma,nrows*chunk,p) # Boostrap Resampling
            end
            ns = (k*nrows+1):((k+1)*nrows)
            medians[:,i] = nanmedian(dbs[ns,1], dbs[ns,2] ./ dbs[ns,3], xmin, xmax, nbins)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmedian(medians,dim=2) # Median-of-medians
        el = m .- pctile(medians,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(medians,97.5,dim=2) .- m # Upper bound of central 95% CI

        return (c, m, el, eu)
    end
    export bin_bsr_ratio_medians


## --- Quick Monte Carlo binning/interpolation functions

    function mcfit(x::AbstractVector, sx::AbstractVector, y::AbstractVector, sy::AbstractVector,
        xmin::Number, xmax::Number, nbins::Integer=10; binwidth::Number=(xmax-xmin)/(nbins-1), minrows::Number=100000)
        # (c,m)=monteCarloFit(x,sx,y,sy,xmin,xmax,nbins,binwidth)
        # Run a simplified Monte Carlo fit with nbins of witdth binwidth between xmin and xmax

        # Fill in variances where not provided explicitly
        sx[isnan.(sx) .& .!isnan.(x)] .= nanstd(x)
        sy[isnan.(sy) .& .!isnan.(y)] .= nanstd(y)

        # Remove missing data
        hasdata = .!(isnan.(x) .| isnan.(y) .| isnan.(sx) .| isnan.(sy))
        x = x[hasdata]; y = y[hasdata];
        sx = sx[hasdata]; sy = sy[hasdata];

        # Increase x uncertainty if x sampling is sparse
        xsorted = sort(x)
        minerr = maximum(xsorted[2:end] - xsorted[1:end-1]) / 2
        sx[sx .< minerr] .= minerr

        # Run the Monte Carlo
        halfwidth = binwidth / 2
        c = collect(xmin:(xmax-xmin)/(nbins-1):xmax)

        nsims = ceil(Int, minrows/length(x))
        xm = repeat(x,nsims) + randn(length(x)*nsims) .* repeat(sx,nsims)
        ym = repeat(y,nsims) + randn(length(y)*nsims) .* repeat(sy,nsims)

        m = fill(NaN, nbins)
        for i = 1:nbins
            inbin = ym[(xm .> c[i] .- halfwidth) .& (xm .< c[i] .+ halfwidth)]
            m[i] = nanmean(inbin)
        end

        return (c, m)
    end
    export mcfit

## --- Downsample an image / array

    function downsample(matrix::Array, factor::Integer, jfactor=factor::Integer)
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

    """
    ```julia
    k = invweight(lat::AbstractVector{<:Number}, lon::AbstractVector{<:Number}, age::AbstractVector{<:Number};
        \tlp::Number=2, spatialscale::Union{Number,AbstractVector{<:Number}}=1.8, agescale::Union{Number,AbstractVector{<:Number}}=38.0)
    ```

    Find the inverse weights `k` (proportional to spatiotemporal sample density) for
    a set of geological samples with specified latitude (`lat`), logitude (`lon`),
    and `age` (of crystallization, deposition, etc.).

    The default `spatialscale` and `agescale` are taken from Keller and Schoene 2012.
    However, alternative scalings can be supplied. If an array is supplied for either
    `spatialscale`, `agescale`, or both, a 3-d matrix of `k` values will be returned,
    with dimensions length(`spatialscale`)*length(`agescale`)*nrows.
    """
    function invweight(lat::AbstractVector{<:Number}, lon::AbstractVector{<:Number}, age::AbstractVector{<:Number};
        lp::Number=2, spatialscale::Union{Number,AbstractVector{<:Number}}=1.8, agescale::Union{Number,AbstractVector{<:Number}}=38.0)

        # Check if there is lat, lon, and age data
        nodata = isnan.(lat) .| isnan.(lon) .| isnan.(age)

        # Convert lat and lon to radians
        latr = lat/180*pi
        lonr = lon/180*pi
        spatialscalr= spatialscale/180*pi

        # Precalculate some sines and cosines
        latsin = sin.(latr)
        latcos = cos.(latr)

        # Allocate and fill ks
        if isa(spatialscale, Number) && isa(agescale, Number)
            k = Array{Float64}(undef,length(lat))
            @showprogress 1 "Calculating weights: " for i=1:length(lat)
                if nodata[i] # If there is missing data, set k=inf for weight=0
                    k[i] = Inf
                else # Otherwise, calculate weight
                    k[i] = nansum( @avx @. 1.0 / ( (acos(min( latsin[i] * latsin + latcos[i] * latcos * cos(lonr[i] - lonr), 1.0 ))/spatialscalr)^lp + 1.0) + 1.0 / ((abs(age[i] - age)/agescale)^lp + 1.0) )
                end
            end
        else
            k = Array{Float64}(undef,length(spatialscale),length(agescale),length(lat))
            spatialdistr = similar(lat)
            @showprogress 1 "Calculating weights: " for i=1:length(lat)
                if nodata[i] # If there is missing data, set k=inf for weight=0
                    k[:,:,i] .= Inf
                else # Otherwise, calculate weight
                    @avx @. spatialdistr = acos(min( latsin[i] * latsin + latcos[i] * latcos * cos(lonr[i] - lonr), 1.0 ))
                    Threads.@threads for g = 1:length(spatialscale)
                        for h = 1:length(agescale)
                            k[g,h,i] = nansum( @avx @. 1.0 / ( (spatialdistr/spatialscalr[g])^lp + 1.0) + 1.0 / ((abs(age[i] - age)/agescale[h])^lp + 1.0) )
                        end
                    end
                end
            end
        end
        return k
    end
    export invweight


    """
    ```julia
    k = invweight_location(lat::Array{<:Number}, lon::Array{<:Number};
        \tlp::Number=2, spatialscale::Number=1.8)
    ```

    Find the inverse weights `k` (proportional to spatial sample density) for
    a set of geological samples with specified latitude (`lat`), and logitude (`lon`).
    """
    function invweight_location(lat::Array{<:Number}, lon::Array{<:Number};
        lp::Number=2, spatialscale::Number=1.8)

        # Check if there is lat, lon data
        nodata = isnan.(lat) .| isnan.(lon)

        # Convert lat and lon to radians
        latr = lat/180*pi
        lonr = lon/180*pi
        spatialscalr= spatialscale/180*pi

        # Precalculate some sines and cosines
        latsin = sin.(latr)
        latcos = cos.(latr)

        # Allocate and fill ks
        k = Array{Float64}(undef,length(lat))
        @showprogress 1 "Calculating weights: " for i=1:length(lat)
            if nodata[i] # If there is missing data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( @avx @. 1.0 / ( (acos(min( latsin[i] * latsin + latcos[i] * latcos * cos(lonr[i] - lonr), 1.0 ))/spatialscalr)^lp + 1.0) )
            end
        end
        return k
    end
    export invweight_location


    """
    ```julia
    k = invweight(nums::Array{<:Number}, scale::Number; lp=2)
    ```

    Find the inverse weights for a single array `nums` for a given `scale`, and
    exponent `lp` (default lp = 2).

    Returns an array k where k[i] is the "inverse weight" for element i of the
    input array.
    """
    function invweight(nums::Array{<:Number}, scale::Number; lp=2)
        # Check if there is data
        nodata = isnan.(nums)

        k = Array{Float64}(undef,length(nums))
        @showprogress 1 "Calculating weights: " for i=1:length(nums)
            if nodata[i] # If there is missing data, set k=inf for weight=0
                k[i] = Inf
            else # Otherwise, calculate weight
                k[i] = nansum( @avx @. 1.0 / ( (abs(nums[i] - nums)/scale)^lp + 1.0) )
            end
        end
        return k
    end
    export invweight


    """
    ```julia
    k = invweight_age(age::Array{<:Number}; lp::Number=2, agescale::Number=38.0)
    ```

    Find the inverse weights `k` (proportional to temporal sample density) for
    a set of geological samples with specified `age` (of crystallization, deposition, etc.).
    """
    function invweight_age(age::Array{<:Number}; lp::Number=2, agescale::Number=38.0)
        return invweight(age, agescale, lp=lp)
    end
    export invweight_age

## --- End of file
