## --- Bootstrap resampling

    # Kernel functions for bsr!
    @inline function uniform(rng::AbstractRNG, mu, halfwidth)
        mu + halfwidth * (2*rand(rng)-1)
    end
    export uniform
    @inline function triangular(rng::AbstractRNG, mu, halfwidth)
        mu + halfwidth * (rand(rng) - rand(rng))
    end
    export triangular
    @inline function gaussian(rng::AbstractRNG, mu, sigma)
        mu + sigma * randn(rng)
    end
    export gaussian

    """
    ```julia
    bsr!([f::Function], resampled::Array, index::Vector{Int}, data, sigma, p;
        \trng::AbstractRNG=MersenneTwister()
    )
    ```

    Fill `resampled` with data boostrap resampled from a (sample-per-row / element-per-column)
    dataset `data` with uncertainties `sigma` and resampling probabilities `p`, optionally using
    random numbers generated by `f` where `f` is a function of the form `f(rng, data[i], sigma[i])`
    """
    function bsr!(f::Function, resampled::Array, index::Vector{Int}, data::AbstractArray, sigma::Number, p::Number; rng::AbstractRNG=MersenneTwister())
        # Prepare
        ndata = size(data,1)
        nrows = size(resampled,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        nrows_accepted = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            @inbounds for i=1:ndata
                if rand(rng) < p
                    nrows_accepted += 1
                    index[nrows_accepted] = i
                    if nrows_accepted == nrows
                        break
                    end
                end
            end
            nrows_new = min(nrows_accepted - n, nrows - n)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                # Optimized inner loop
                @inbounds for i = 1:nrows_new
                    row = index[n+i]
                    resampled[n+i,j] = f(rng, data[row,j], sigma)
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    function bsr!(f::Function, resampled::Array, index::Vector{Int}, data::AbstractArray, sigma::Number, p::AbstractVector; rng::AbstractRNG=MersenneTwister())
        # Prepare
        ndata = size(data,1)
        nrows = size(resampled,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        nrows_accepted = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            @inbounds for i=1:ndata
                if rand(rng) < p[i]
                    nrows_accepted += 1
                    index[nrows_accepted] = i
                    if nrows_accepted == nrows
                        break
                    end
                end
            end
            nrows_new = min(nrows_accepted - n, nrows - n)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                # Optimized inner loop
                @inbounds for i = 1:nrows_new
                    row = index[n+i]
                    resampled[n+i,j] = f(rng, data[row,j], sigma)
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    function bsr!(f::Function, resampled::Array, index::Vector{Int}, data::AbstractArray, sigma::AbstractArray, p::Number; rng::AbstractRNG=MersenneTwister())
        # Prepare
        ndata = size(data,1)
        nrows = size(resampled,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        nrows_accepted = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            @inbounds for i=1:ndata
                if rand(rng) < p
                    nrows_accepted += 1
                    index[nrows_accepted] = i
                    if nrows_accepted == nrows
                        break
                    end
                end
            end
            nrows_new = min(nrows_accepted - n, nrows - n)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                # Optimized inner loop
                @inbounds for i = 1:nrows_new
                    row = index[n+i]
                    resampled[n+i,j] = f(rng, data[row,j], sigma[row,j])
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    function bsr!(f::Function, resampled::Array, index::Vector{Int}, data::AbstractArray, sigma::AbstractArray, p::AbstractVector; rng::AbstractRNG=MersenneTwister())
        # Prepare
        ndata = size(data,1)
        nrows = size(resampled,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        nrows_accepted = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            @inbounds for i=1:ndata
                if rand(rng) < p[i]
                    nrows_accepted += 1
                    index[nrows_accepted] = i
                    if nrows_accepted == nrows
                        break
                    end
                end
            end
            nrows_new = min(nrows_accepted - n, nrows - n)

            # Columns go in outer loop because of column major indexing
            for j=1:ncolumns
                # Optimized inner loop
                @inbounds for i = 1:nrows_new
                    row = index[n+i]
                    resampled[n+i,j] = f(rng, data[row,j], sigma[row,j])
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    # default method if f not specified
    bsr!(resampled::Array, index::Vector{Int}, data::AbstractArray, sigma, p; rng=MersenneTwister()) = bsr!(gaussian, resampled, index, data, sigma, p, rng=rng)
    export bsr!

    """
    ```julia
    resampled = bsresample(data::AbstractArray, sigma, nrows, [p];
        \t kernel = gaussian,
        \t rng = MersenneTwister(),
        \t return_index = false
    )
    ```
    Bootstrap resample a (sample-per-row / element-per-column) array of `data`
    with uncertainties `sigma` and resampling probabilities `p`
    """
    function bsresample(data::AbstractArray, sigma, nrows, p=min(0.2,nrows/size(data,1));
            kernel = gaussian,
            rng = MersenneTwister(),
            return_index = false
        )
        index = Array{Int}(undef, nrows)
        resampled = Array{float(eltype(data))}(undef, nrows, size(data,2))
        bsr!(kernel, resampled, index, data, sigma, p, rng=rng)
        if return_index
            return resampled, index
        else
            return resampled
        end
    end
    """
    ```julia
    resampled = bsresample(dataset::Dict, nrows, [elements], [p];
        \t kernel = gaussian,
        \t rng = MersenneTwister()
    )
    ```
    Bootstrap resample a dictionary-based `dataset` with uncertainties stored either
    in `dataset["err"]` or `dataset["[variable]_sigma"]`
    """
    function bsresample(dataset::Dict, nrows, elements=dataset["elements"], p=min(0.2,nrows/length(in[elements[1]]));
            kernel = gaussian,
            rng = MersenneTwister()
        )
        # 2d array of nominal values
        data = unelementify(dataset, elements, floatout=true)

        # 2d array of absolute 1-sigma uncertainties
        if haskey(dataset, "err") && isa(dataset["err"], Dict)
            sigma = unelementify(dataset["err"], elements, floatout=true)
        else
            sigma = unelementify(dataset, elements.*"_sigma", floatout=true)
        end

        # Resample
        sdata = bsresample(data, sigma, nrows, p, kernel=kernel, rng=rng)
        return elementify(sdata, elements, skipstart=0)
    end
    export bsresample


    function randsample!(resampled::AbstractArray, data::AbstractArray, nrows::Integer, p::Number, rng::AbstractRNG=MersenneTwister(), buffer::Vector{Int}=Array{Int}(undef,size(data,1)))
        # Prepare
        ndata = size(data,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:ndata
                if rand(rng) < p
                    nrows_accepted += 1
                    buffer[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n)

            # Columns go in outer loop because of column major indexing
            @inbounds @simd for j=1:ncolumns
                # Optimized inner loop
                for i = 1:nrows_new
                    resampled[n+i,j] = data[buffer[i],j]
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end
    function randsample!(resampled::AbstractArray, data::AbstractArray, nrows::Integer, p::AbstractVector, rng::AbstractRNG=MersenneTwister(), buffer::Vector{Int}=Array{Int}(undef,size(data,1)))
        # Prepare
        ndata = size(data,1)
        ncolumns = size(resampled,2)

        # Resample
        n = 0
        while n < nrows

            # Compare acceptance probability p against Unif(0,1)
            nrows_accepted = 0
            @inbounds for i=1:ndata
                if rand(rng) < p[i]
                    nrows_accepted += 1
                    buffer[nrows_accepted] = i
                end
            end
            nrows_new = min(nrows_accepted, nrows - n)

            # Columns go in outer loop because of column major indexing
            @inbounds @simd for j=1:ncolumns
                # Optimized inner loop
                for i = 1:nrows_new
                    resampled[n+i,j] = data[buffer[i],j]
                end
            end

            # Keep track of current filled rows
            n += nrows_new
        end

        return resampled
    end

    """
    ```julia
    randsample(data::Array, nrows, [p])
    ```
    Bootstrap resample (without uncertainty) a `data` array to length `nrows`.
    Optionally provide weights `p` either as a vector (one-weight-per-sample) or scalar.
    """
    function randsample(data::AbstractArray, nrows::Integer, p=min(0.2,nrows/size(data,1));
            rng::AbstractRNG=MersenneTwister(),
            buffer::Vector{Int}=Array{Int}(undef,size(data,1))
        )
        resampled = Array{eltype(data)}(undef,nrows,size(data,2))
        return randsample!(resampled, data, nrows, p, rng, buffer)
    end
    # Second method for randsample that takes a dictionary as input
    """
    ```julia
    randsample(dataset::Dict, nrows, [elements], [p])
    ```
    Bootstrap resample (without uncertainty) a `dataset` dict to length `nrows`.
    Optionally provide weights `p` either as a vector (one-weight-per-sample) or scalar.
    """
    function randsample(dataset::Dict, nrows::Integer, elements=in["elements"], p=min(0.2,nrows/length(in[elements[1]])))
        data = unelementify(dataset, elements, floatout=true)
        sdata = randsample(data, nrows, p)
        return elementify(sdata, elements, skipstart=0)
    end
    export randsample


## --- Bin a dataset by a given independent variable

    """
    ```julia
    (bincenters, N) = bincounts(x::AbstractArray, xmin::Number, xmax::Number, nbins::Integer)
    ```
    Tally the number of samples that fall into each of `nbins` equally spaced `x`
    bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function bincounts(x::AbstractArray, xmin::Number, xmax::Number, nbins::Integer)
        # Tally the number of samples (either resampled or corrected/original) that fall into each bin
        binwidth = (xmax-xmin)/nbins
        bincenters = (xmin+binwidth/2):binwidth:(xmax-binwidth/2)

        # Calculate index from x value
        scalefactor = nbins / (xmax - xmin)
        index = (x .- xmin) .* scalefactor

        # Add up the results
        N = fill(0,nbins)
        @inbounds for i in index
            if 0 < i <= nbins
                N[ceil(Int, i)] += 1
            end
        end
        return bincenters, N
    end
    export bincounts

    """
    ```julia
    (c,m,e) = binmeans(x, y, xmin, xmax, nbins, [weight]; resamplingratio::Number=1)
    ```
    The means (ignoring NaNs) of `y` values binned by `x`, into each of `nbins`
    equally spaced `x` bins between `xmin` and `xmax`, returning bincenters,
    means, and standard errors of the mean.

    To calculate binned medians only (without uncertainties), see `nanmean`
    """
    function binmeans(x::AbstractArray, y::AbstractArray, xmin::Number, xmax::Number, nbins::Integer; resamplingratio::Number=1)
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
            if (0 < index_float[i] <= nbins) && !isnan(y[i])
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
    function binmeans(x::AbstractArray, y::AbstractArray, min::Number, max::Number, nbins::Integer, weight::AbstractArray; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = range(min,max,length=nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        means = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (binedges[i] .< x .<= binedges[i+1]) .& (y.==y)
            w = ProbabilityWeights(weight[t])
            means[i] = mean(y[t], w)
            errors[i] = std(y[t], w, corrected=true) * sqrt(resamplingratio) / sqrt(count(t))
        end

        return bincenters, means, errors
    end
    export binmeans

    """
    ```julia
    (c,m,e) = binmedians(x, y, xmin, xmax, nbins; resamplingratio::Number=1)
    ```

    The medians (ignoring NaNs) of `y` values binned by `x`, into each of `nbins`
    equally spaced `x` bins between `xmin` and `xmax`, returning bincenters, medians,
    and equivalent standard errors of the mean (1.4828 * median abolute deviation)

    To calculate binned medians only (without uncertainties), see `nanmedian`
    """
    function binmedians(x::AbstractArray, y::AbstractArray, min::Number, max::Number, nbins::Integer; resamplingratio::Number=1)
        binwidth = (max-min)/nbins
        binedges = range(min,max,length=nbins+1)
        bincenters = (min+binwidth/2):binwidth:(max-binwidth/2)

        medians = Array{Float64}(undef,nbins)
        errors = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (binedges[i] .< x .<= binedges[i+1]) .& (y.==y)
            medians[i] = median(y[t])
            errors[i] = 1.4826 * nanmad(y[t]) * sqrt(resamplingratio / count(t))
        end

        return bincenters, medians, errors
    end
    export binmedians

## --- Bin bootstrap resampled data

    """
    ```julia
    bin_bsr([f!::Function], x::Vector, y::VecOrMat, xmin, xmax, nbins, [w];
        \tx_sigma = zeros(size(x)),
        \ty_sigma = zeros(size(y)),
        \tnresamplings = 1000,
        \tsem = :sigma,
        \tp = 0.2
    )
    ```
    Returns the bincenters `c`, means or medians `m`, and uncertainties of the
    mean or median for a variable `y` binned by independent variable `x` into
    `nbins` equal bins between `xmin` and `xmax`, after `nresamplings` boostrap
    resamplings with acceptance probability `p`.

    If a 2-d array (matrix) of `y` values is provided, each column will be treated
    as a separate variable, means and uncertainties will be returned column-wise.

    Optional keyword arguments and defaults:

        x_sigma = zeros(size(x))

    A vector representing the uncertainty (standard deviation) of each x value

        y_sigma = zeros(size(y))

    A vector representing the uncertainty (standard deviation) of each y value

        nresamplings = 1000

    The number of resamplings to conduct

        sem = :sigma

    Format of the uncertainty estimate of the distribution of the mean.
    If `:sigma` is chosen, a tuple of three vectors `(c, m, e)` will be returned,
    where `e` is the standard error of the mean.
    If `:CI` or `:pctile` is chosen, a tuple of four vectors `(c, m, el, eu)`
    will be returned, where `el` and `eu` are the lower and upper bounds of the 95%
    confidence interval.

        p = 0.2

    Resampling probabilities, either as a scalar or a vector of the same length as `x`


    # Examples:
    ```julia
    (c,m,e) = bin_bsr(nanmedian!, x, y, 0, 4000, 40, x_sigma=0.05x, p=probability, sem=:sigma)
    ```
    ```julia
    (c,m,el,eu) = bin_bsr(nanmean!, x, y, 0, 4000, 40, x_sigma=0.05x, p=probability, sem=:pctile)
    ```
    """
    function bin_bsr(f!::Function, x::AbstractVector, y::AbstractVector, xmin, xmax, nbins::Integer;
            x_sigma = zeros(size(x)),
            y_sigma = zeros(size(y)),
            nresamplings = 1000,
            sem = :credibleinterval,
            p = 0.2
        )

        data = hcat(x, y)
        sigma = hcat(x_sigma, y_sigma)
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        dbs = Array{Float64}(undef, nrows, ncols)
        index = Array{Int}(undef, nrows) # Must be preallocated even if we don't want it later
        means = Array{Float64}(undef, nbins, nresamplings)
        rng = MersenneTwister()
        N = Array{Int}(undef, nbins)
        # Resample
        for i=1:nresamplings
            bsr!(dbs, index, data, sigma, p, rng=rng) # Boostrap Resampling
            f!(view(means,:,i), N, view(dbs,:,1), view(dbs,:,2), xmin, xmax, nbins)
        end

        # Return summary of results
        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        if sem == :sigma
            # Standard deviation of means (sem)
            e = nanstd(means,dim=2)
            return c, m, e
        elseif sem == :credibleinterval || sem == :CI || sem == :pctile
            # Lower bound of central 95% CI of means
            el = m .- pctile(means,2.5,dim=2)
            # Upper bound of central 95% CI of means
            eu = pctile(means,97.5,dim=2) .- m
            return c, m, el, eu
        else
            return c, means
        end
    end
    function bin_bsr(f!::Function, x::AbstractVector, y::AbstractMatrix, xmin, xmax, nbins::Integer;
            x_sigma = zeros(size(x)),
            y_sigma = zeros(size(y)),
            nresamplings = 1000,
            sem = :credibleinterval,
            p = 0.2
        )

        data = hcat(x, y)
        sigma = hcat(x_sigma, y_sigma)
        dtype = float(eltype(data))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        dbs = Array{dtype}(undef, nrows, ncols)
        means = Array{dtype}(undef, nbins, nresamplings, size(y,2))
        index = Array{Int}(undef, nrows) # Must be preallocated even if we don't want it later
        rng = MersenneTwister()
        N = Array{Int}(undef, nbins, size(y,2))
        # Resample
        for i=1:nresamplings
            bsr!(dbs, index, data, sigma, p, rng=rng) # Boostrap Resampling
            f!(view(means,:,i,:), N, view(dbs,:,1), view(dbs,:,2:1+size(y,2)), xmin, xmax, nbins)
        end

        # Return summary of results
        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        if sem == :sigma
            m = Array{dtype}(undef, nbins, size(y,2))
            e = Array{dtype}(undef, nbins, size(y,2))
            for j = 1:size(y,2)
                m[:,j] .= nanmean(view(means,:,:,j),dim=2) # Mean-of-means
                e[:,j] .= nanstd(view(means,:,:,j),dim=2) # Standard deviation of means (sem)
            end
            return c, m, e
        elseif sem == :credibleinterval || sem == :CI || sem == :pctile
            m = Array{dtype}(undef, nbins, size(y,2))
            el = Array{dtype}(undef, nbins, size(y,2))
            eu = Array{dtype}(undef, nbins, size(y,2))
            for j = 1:size(y,2)
                m[:,j] .= nanmean(view(means,:,:,j),dim=2) # Mean-of-means
                el[:,j] .= m[:,j] .- pctile(view(means,:,:,j),dim=2)
                eu[:,j] .= pctile(view(means,:,:,j),dim=2) .- m[:,j]
            end
            return c, m, el, eu
        else
            return c, means
        end
    end
    function bin_bsr(f!::Function, x::AbstractVector, y::AbstractVector, xmin, xmax, nbins::Integer, w::AbstractVector;
            x_sigma = zeros(size(x)),
            y_sigma = zeros(size(x)),
            nresamplings = 1000,
            sem = :credibleinterval,
            p = 0.2
        )

        data = hcat(x, y, w)
        sigma = hcat(x_sigma, y_sigma, zeros(size(w)));

        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        dbs = Array{Float64}(undef, nrows, ncols)
        means = Array{Float64}(undef, nbins, nresamplings)
        index = Array{Int}(undef, nrows) # Must be preallocated even if we don't want it later
        rng = MersenneTwister()
        N = Array{Int}(undef, nbins)
        # Resample
        for i=1:nresamplings
            bsr!(dbs, index, data, sigma, p, rng=rng) # Boostrap Resampling
            f!(view(means,:,i), N, view(dbs,:,1), view(dbs,:,2), view(dbs,:,3), xmin, xmax, nbins)
        end

        # Return summary of results
        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        if sem == :sigma
            # Standard deviation of means (sem)
            e = nanstd(means,dim=2)
            return c, m, e
        elseif sem == :credibleinterval || sem == :CI || sem == :pctile
            # Lower bound of central 95% CI of means
            el = m .- pctile(means,2.5,dim=2)
            # Upper bound of central 95% CI of means
            eu = pctile(means,97.5,dim=2) .- m
            return c, m, el, eu
        else
            return c, means
        end
    end
    bin_bsr(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:sigma, p=0.2) =
        bin_bsr(nanmean!,x,y,xmin,xmax,nbins,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    bin_bsr(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer, w::AbstractVector; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:sigma, p=0.2) =
        bin_bsr(nanmean!,x,y,xmin,xmax,nbins,w,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    export bin_bsr

    bin_bsr_means(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr(nanmean!,x,y,xmin,xmax,nbins,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    bin_bsr_means(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer, w::AbstractVector; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr(nanmean!,x,y,xmin,xmax,nbins,w,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    export bin_bsr_means

    bin_bsr_medians(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr(nanmedian!,x,y,xmin,xmax,nbins,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    bin_bsr_medians(x::AbstractVector, y::AbstractVecOrMat, xmin, xmax, nbins::Integer, w::AbstractVector; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr(nanmedian!,x,y,xmin,xmax,nbins,w,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,sem=sem,p=p)
    export bin_bsr_medians

    """
    ```julia
    (c, m, el, eu) = bin_bsr_ratios([f!::Function], x::Vector, num::Vector, denom::Vector, xmin, xmax, nbins, [w];
        \tx_sigma = zeros(size(x)),
        \tnum_sigma = zeros(size(num)),
        \tdenom_sigma = zeros(size(denom)),
        \tnresamplings = 1000,
        \tp::Union{Number,Vector} = 0.2
    )
    ```

    Returns the bincenters `c`, means `m`, as well as upper (`el`) and lower (`eu`) 95% CIs of the mean
    for a ratio `num`/`den` binned by `x` into `nbins` equal bins between `xmin` and `xmax`,
    after `nresamplings` boostrap resamplings with acceptance probability `p`.
    """
    function bin_bsr_ratios(f!::Function, x::AbstractVector, num::AbstractVector, denom::AbstractVector, xmin, xmax, nbins::Integer;
            x_sigma::AbstractVector=zeros(size(x)),
            num_sigma::AbstractVector=zeros(size(num)),
            denom_sigma::AbstractVector=zeros(size(denom)),
            nresamplings=1000,
            p::Union{Number,AbstractVector}=0.2
        )

        data = hcat(x, num, denom)
        sigma = hcat(x_sigma, num_sigma, denom_sigma)
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        dbs = Array{Float64}(undef, nrows, ncols)
        index = Array{Int}(undef, nrows) # Must be preallocated even if we don't want it later
        means = Array{Float64}(undef, nbins, nresamplings)
        fractions = Array{Float64}(undef, nrows)
        fraction_means = Array{Float64}(undef, nbins)
        rng = MersenneTwister()
        N = Array{Int}(undef, nbins) # Array of bin counts -- Not used but preallocated for speed
        # Resample
        for i=1:nresamplings
            bsr!(dbs, index, data, sigma, p, rng=rng) # Boostrap Resampling
            @views @avx @. fractions = dbs[:,2] / (dbs[:,2] + dbs[:,3])
            f!(fraction_means, N, view(dbs,:,1), fractions, xmin, xmax, nbins)
            @. means[:,i] = fraction_means / (1 - fraction_means)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    function bin_bsr_ratios(f!::Function, x::AbstractVector, num::AbstractVector, denom::AbstractVector, xmin, xmax, nbins::Integer, w::AbstractVector;
            x_sigma::AbstractVector=zeros(size(x)),
            num_sigma::AbstractVector=zeros(size(num)),
            denom_sigma::AbstractVector=zeros(size(denom)),
            nresamplings=1000,
            p::Union{Number,AbstractVector}=0.2
        )

        data = hcat(x, num, denom, w)
        sigma = hcat(x_sigma, num_sigma, denom_sigma, zeros(size(w)))
        binwidth = (xmax-xmin)/nbins
        nrows = size(data,1)
        ncols = size(data,2)

        # Preallocate
        dbs = Array{Float64}(undef, nrows, ncols)
        index = Array{Int}(undef, nrows) # Must be preallocated even if we don't want it later
        means = Array{Float64}(undef, nbins, nresamplings)
        fractions = Array{Float64}(undef, nrows)
        fraction_means = Array{Float64}(undef, nbins)
        rng = MersenneTwister()
        W = Array{Float64}(undef, nbins) # Array of bin weights -- Not used but preallocated for speed
        # Resample
        for i=1:nresamplings
            bsr!(dbs, index, data, sigma, p, rng=rng) # Boostrap Resampling
            @views @avx @. fractions = dbs[:,2] / (dbs[:,2] + dbs[:,3])
            f!(fraction_means, W, view(dbs,:,1), fractions, view(dbs,:,4), xmin, xmax, nbins)
            @. means[:,i] = fraction_means / (1 - fraction_means)
        end

        c = (xmin+binwidth/2):binwidth:(xmax-binwidth/2) # Bin centers
        m = nanmean(means,dim=2) # Mean-of-means
        el = m .- pctile(means,2.5,dim=2) # Lower bound of central 95% CI
        eu = pctile(means,97.5,dim=2) .- m # Upper bound of central 95% CI

        return c, m, el, eu
    end
    bin_bsr_ratios(x::AbstractVector, num::AbstractVector, denom::AbstractVector, xmin, xmax, nbins::Integer; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr_ratios(nanmean!,x,num,denom,xmin,xmax,nbins,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,p=p)
    bin_bsr_ratios(x::AbstractVector, num::AbstractVector, denom::AbstractVector, xmin, xmax, nbins::Integer, w::AbstractVector; x_sigma=zeros(size(x)), y_sigma=zeros(size(y)), nresamplings=1000, sem=:pctile, p=0.2) =
        bin_bsr_ratios(nanmean!,x,num,denom,xmin,xmax,nbins,w,x_sigma=x_sigma,y_sigma=y_sigma,nresamplings=nresamplings,p=p)
    export bin_bsr_ratios


## --- Quick Monte Carlo binning/interpolation functions

    function mcfit(x::AbstractVector, σx::AbstractVector, y::AbstractVector, σy::AbstractVector,
            xmin::Number, xmax::Number, nbins::Integer=10;
            binwidth::Number=(xmax-xmin)/(nbins-1),
            minrows::Number=100000
        )
        # Run a simplified Monte Carlo fit with nbins of witdth binwidth between xmin and xmax

        # Remove missing data
        hasdata = .!(isnan.(x) .| isnan.(y))
        x = x[hasdata]
        y = y[hasdata]
        σx = σx[hasdata]
        σy = σy[hasdata]

        # Fill in variances where not provided explicitly
        σx[isnan.(σx)] .= nanstd(x)
        σy[isnan.(σy)] .= nanstd(y)

        # Increase x uncertainty if x sampling is sparse
        xsorted = sort(x)
        minerr = maximum(xsorted[2:end] - xsorted[1:end-1]) / 2
        σx[σx .< minerr] .= minerr

        # Bin centers
        c = xmin:(xmax-xmin)/(nbins-1):xmax
        halfwidth = binwidth / 2

        # Run the Monte Carlo
        N = fill(0, nbins)
        m = Array{float(eltype(y))}(undef, nbins)
        xresampled = similar(x, float(eltype(x)))
        yresampled = similar(y, float(eltype(y)))
        @inbounds for n = 1:ceil(Int, minrows/length(x))
            xresampled .= x .+ σx .* randn!(xresampled)
            yresampled .= y .+ σy .* randn!(yresampled)
            for j = 1:nbins
                l = (c[j] - halfwidth)
                u = (c[j] + halfwidth)
                for i ∈ eachindex(xresampled)
                    if l < xresampled[i] < u
                        m[j] += yresampled[i]
                        N[j] += 1
                    end
                end
            end
        end
        m ./= N

        return c, m
    end
    export mcfit

## --- Downsample an image / array

    function downsample(A::AbstractArray, factor::Integer, jfactor=factor::Integer)
        if ndims(A)==2
            rows = size(A,1) ÷ factor
            cols = size(A,2) ÷ jfactor

            result = Array{eltype(A)}(undef,rows,cols)
            @avx for i=1:rows
                for j=1:cols
                    result[i,j]=A[i*factor,j*jfactor]
                end
            end
        else
            result = A[factor:factor:end]
        end
        return result
    end
    export downsample

## --- Spatiotemporal sample weighting

    """
    ```julia
    k = invweight(lat::AbstractArray, lon::AbstractArray, age::AbstractArray;
        \tlp::Number=2,
        \tspatialscale=1.8,
        \tagescale=38.0
    )
    ```

    Find the inverse weights `k` (proportional to spatiotemporal sample density) for
    a set of geological samples with specified latitude (`lat`), logitude (`lon`),
    and `age` (of crystallization, deposition, etc.).

    The default `spatialscale` and `agescale` are taken from Keller and Schoene 2012.
    However, alternative scalings can be supplied. If an array is supplied for either
    `spatialscale`, `agescale`, or both, a 3-d matrix of `k` values will be returned,
    with dimensions length(`spatialscale`)*length(`agescale`)*nrows.
    """
    function invweight(lat::AbstractArray, lon::AbstractArray, age::AbstractArray;
        lp::Number=2, spatialscale=1.8, agescale=38.0)


        # Check if there is lat, lon, and age data
        nodata = vec(isnan.(lat) .| isnan.(lon) .| isnan.(age))

        # Convert lat and lon to radians
        latr = vec(lat/180*pi)
        lonr = vec(lon/180*pi)
        spatialscalr = spatialscale/180*pi

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
    k = invweight_location(lat::AbstractArray, lon::AbstractArray;
        \tlp::Number=2,
        \tspatialscale::Number=1.8
    )
    ```

    Find the inverse weights `k` (proportional to spatial sample density) for
    a set of geological samples with specified latitude (`lat`), and logitude (`lon`).
    """
    function invweight_location(lat::AbstractArray, lon::AbstractArray;
        lp::Number=2, spatialscale::Number=1.8)

        # Check if there is lat, lon data
        nodata = vec(isnan.(lat) .| isnan.(lon))

        # Convert lat and lon to radians
        latr = vec(lat/180*pi)
        lonr = vec(lon/180*pi)
        spatialscalr = spatialscale/180*pi

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
    k = invweight(nums::AbstractArray, scale::Number; lp=2)
    ```

    Find the inverse weights for a single array `nums` for a given `scale`, and
    exponent `lp` (default lp = 2).

    Returns an array k where k[i] is the "inverse weight" for element i of the
    input array.
    """
    function invweight(nums::AbstractArray, scale::Number; lp=2)
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
    k = invweight_age(age::AbstractArray; lp::Number=2, agescale::Number=38.0)
    ```

    Find the inverse weights `k` (proportional to temporal sample density) for
    a set of geological samples with specified `age` (of crystallization, deposition, etc.).
    """
    function invweight_age(age::AbstractArray; lp::Number=2, agescale::Number=38.0)
        return invweight(age, agescale, lp=lp)
    end
    export invweight_age

## --- End of file
