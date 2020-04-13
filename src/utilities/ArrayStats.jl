## --- Weighted mean of an array

    """
    ```julia
    (wx, wσ, mswd) = awmean(x, σ)
    ```
    Weighted mean, absent the MSWD correction to uncertainty.
    """
    function awmean(x, σ)
        n = length(x)

        if n == 1
            wx = x[1]
            mswd = NaN
            wσ = σ[1]
        else
            sum_of_values = sum_of_weights = χ2 = 0.0
            for i=1:n
                sum_of_values += x[i] / (σ[i]*σ[i])
                sum_of_weights += 1 / (σ[i]*σ[i])
            end
            wx = sum_of_values / sum_of_weights

            for i=1:n
                χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
            end
            mswd = χ2 / (n-1)
            wσ = sqrt(1.0 / sum_of_weights)
        end
        return wx, wσ, mswd
    end
    export awmean

    """
    ```julia
    (wx, wσ, mswd) = gwmean(x, σ)
    ```
    Geochronologist's weighted mean, with "MSWD correction" to uncertainty,
    i.e., wσ is increased by a factor of sqrt(mswd)
    """
    function gwmean(x, σ)
        n = length(x)

        if n == 1
            wx = x[1]
            mswd = NaN
            wσ = σ[1]
        else
            sum_of_values = sum_of_weights = χ2 = 0.0
            for i=1:n
                sum_of_values += x[i] / (σ[i]*σ[i])
                sum_of_weights += 1 / (σ[i]*σ[i])
            end
            wx = sum_of_values / sum_of_weights

            for i=1:n
                χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
            end
            mswd = χ2 / (n-1)
            wσ = sqrt(mswd / sum_of_weights)
        end
        return wx, wσ, mswd
    end
    export gwmean


    """
    ```julia
    MSWD(x, σ)
    ```
    Return the Mean Square of Weighted Deviates (AKA the reduced chi-squared
    statistic) of a dataset with values `x` and one-sigma uncertainties `σ`
    """
    function MSWD(x, σ)
        sum_of_values = sum_of_weights = χ2 = 0.0
        n = length(x)

        for i=1:n
            w = 1 / (σ[i]*σ[i])
            sum_of_values += w * x[i]
            sum_of_weights += w
        end
        wx = sum_of_values / sum_of_weights

        for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end

        return χ2 / (n-1)
    end
    export MSWD

## --- Percentile statistics, excluding NaNs

    """
    ```julia
    pctile(A, p; dim=0)
    ```
    Find the `p`th percentile of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dim`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    function pctile(A, p; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? percentile(A[t],p) : eltype(A)(NaN)
        end
        return result
    end
    export pctile

    """
    ```julia
    inpctile(A, p::Number; dim=0)
    ```
    Return a boolean array that identifies which values of the iterable
    collection `A` fall within the central `p`th percentile, optionally along a
    dimension specified by `dim`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    function inpctile(A, p::Number; dim=0)
        offset = (100 - p) / 2
        return (A .> pctile(A, offset, dim=dim)) .& (A .< pctile(A, 100-offset, dim=dim))
    end
    export inpctile


## --- Summary statistics of arrays with NaNs

    """
    ```julia
    nansum(A; dim=0)
    ```
    Calculate the sum of an indexable collection `A`, ignoring NaNs, optionally
    along a dimension specified by `dim`.
    """
    function nansum(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? sum(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? sum(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? sum(A[t]) : NaN
        end
        return result
    end
    export nansum


    """
    ```julia
    nanminimum(A; dim=0)
    ```
    Find the smallest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dim`.
    """
    function nanminimum(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? minimum(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? minimum(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? minimum(A[t]) : NaN
        end
        return result
    end
    export nanminimum


    """
    ```julia
    nanmaximum(A; dim=0)
    ```
    Find the largest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dim`.
    """
    function nanmaximum(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? maximum(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? maximum(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? maximum(A[t]) : NaN
        end
        return result
    end
    export nanmaximum


    """
    ```julia
    nanextrema(A)
    ```
    Find the extrema (max & min) of an indexable collection `A`, ignoring NaNs.
    """
    function nanextrema(A)
        t = .~ isnan.(A)
        return extrema(A[t])
    end
    export nanextrema


    """
    ```julia
    nanrange(A; dim=0)
    ```
    Calculate the range (max-min) of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dim`.
    """
    function nanrange(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                if any(t)
                    extr = extrema(A[i,t])
                    result[i] = extr[2] - extr[1]
                else
                    result[i] = 0
                end
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                if any(t)
                    extr = extrema(A[t,i])
                    result[i] = extr[2] - extr[1]
                else
                    result[i] = 0
                end
            end
        else
            t = .~ isnan.(A)
            if any(t)
                extr = extrema(A[t])
                result = extr[2] - extr[1]
            else
                result = 0
            end
        end
        return result
    end
    export nanrange


    """
    ```julia
    nanmean(A; dim=0)
    ```
    Calculate the mean, ignoring NaNs, of an indexable collection `A`, optionally
    along a dimension specified by `dim`.
    """
    function nanmean(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? mean(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? mean(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? mean(A[t]) : NaN
        end
        return result
    end
    """
    ```julia
    nanmean(A, W; dim=0)
    ```
    Calculate the weighted mean, ignoring NaNs, of an indexable collection `A`
    with weights `W`, optionally along a dimension specified by `dim`.
    """
    function nanmean(A, W; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? mean(A[i,t], ProbabilityWeights(W[i,t])) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? mean(A[t,i], ProbabilityWeights(W[t,i])) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? mean(A[t], ProbabilityWeights(W[t])) : NaN
        end
        return result
    end


    """
    ```julia
    mu = nanmean(x::AbstractVector{<:Number}, y::AbstractVector{<:Number},
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Return the mean, ignoring NaNs, of `y` values that fall into each of `nbins`
    equally spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function nanmean(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0,nbins)
        mu = fill(float(eltype(y))(0),nbins)
        for i = 1:length(x)
            if (0 < index_float[i] < nbins) && !isnan(y[i])
                index = ceil(Int, index_float[i])
                N[index] += 1
                mu[index] += y[i]
            end
        end
        mu ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    """
    ```julia
    mu = nanmean(x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Where `y` is a (sample-per-row, element-per-column) 2-d array, return a
    (bin-per-row, element-per-column) 2-d array `mu` of the means,
    ignoring NaNs, of `y` values that fall into each of `nbins` equally spaced
    `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function nanmean(x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        bin_index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0, nbins, size(y,2))
        mu = fill(float(eltype(y))(0), nbins, size(y,2))
        for i = 1:length(x)
            if (0 < bin_index_float[i] < nbins)
                bin_index = ceil(Int, bin_index_float[i])
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        mu[bin_index,j] += y[i,j]
                        N[bin_index,j] += 1
                    end
                end
            end
        end
        mu ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    """
    ```julia
    mu = nanmean(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, w::AbstractVector{<:Number},
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Return the weighted mean, ignoring NaNs, of `y` values that fall into each of
    `nbins` equally spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`.
    """
    function nanmean(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, w::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0.0,nbins)
        mu = fill(0.0,nbins)
        for i = 1:length(x)
            if (0 < index_float[i] < nbins) && !isnan(y[i])
                index = ceil(Int, index_float[i])
                N[index] += w[i]
                mu[index] += y[i]*w[i]
            end
        end
        mu ./= N # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    export nanmean



    """
    ```julia
    nanmean!(mu::AbstractVector{<:Number}, x::AbstractVector{<:Number}, y::AbstractVector{<:Number},
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Fill the 1-d array `mu`, with the means, ignoring NaNs, of `y` values that fall into each of
    `nbins` equally spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function nanmean!(mu::AbstractVector{<:Number}, x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        bin_index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0, nbins)
        fill!(mu, 0) # Fill the output array with zeros to start
        for i = 1:length(x)
            if (0 < bin_index_float[i] < nbins) && !isnan(y[i])
                bin_index = ceil(Int, bin_index_float[i])
                N[bin_index] += 1
                mu[bin_index] += y[i]
            end
        end
        mu ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    """
    ```julia
    nanmean!(mu::AbstractMatrix{<:Number}, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Where `y` is a (sample-per-row, element-per-column) 2-d array, fill
    the (bin-per-row, element-per-column) 2-d array `mu`, with the means,
    ignoring NaNs, of `y` values that fall into each of `nbins` equally spaced
    bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function nanmean!(mu::AbstractMatrix{<:Number}, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        bin_index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0, nbins, size(y,2))
        fill!(mu, 0) # Fill the output array with zeros to start
        for i = 1:length(x)
            if (0 < bin_index_float[i] < nbins)
                bin_index = ceil(Int, bin_index_float[i])
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        mu[bin_index,j] += y[i,j]
                        N[bin_index,j] += 1
                    end
                end
            end
        end
        mu ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    """
    ```julia
    nanmean!(mu::AbstractVector{<:Number},x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, w::AbstractVector{<:Number},
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Fill the 1-d array `mu`, with the weighted mean, ignoring NaNs, of `y` values
    that fall into each of `nbins` equally spaced `x` bins between `xmin` and `xmax`,
    aligned with bin edges as `xmin:(xmax-xmin)/nbins:xmax`.
    """
    function nanmean!(mu::AbstractVector{<:Number}, x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, w::AbstractVector{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)
        index_float = (x .- xmin) .* scalefactor

        # Calculate the means for each bin, ignoring NaNs
        N = fill(0.0,nbins)
        fill!(mu, 0) # Fill the output array with zeros to start
        for i = 1:length(x)
            if (0 < index_float[i] < nbins) && !isnan(y[i])
                index = ceil(Int, index_float[i])
                N[index] += w[i]
                mu[index] += y[i]*w[i]
            end
        end
        mu ./= N # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return mu
    end
    export nanmean!


    """
    ```julia
    nanstd(A; dim=0)
    ```
    Calculate the standard deviation, ignoring NaNs, of an indexable collection `A`,
    optionally along a dimension specified by `dim`.
    """
    function nanstd(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? std(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? std(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? std(A[t]) : NaN
        end
        return result
    end
    """
    ```julia
    nanstd(A, W; dim=0)
    ```
    Calculate the weighted standard deviation, ignoring NaNs, of an indexable
    collection `A` with weights `W`, optionally along a dimension specified by `dim`.
    """
    function nanstd(A, W; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? std(A[i,t], ProbabilityWeights(W[i,t]), corrected=false) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? std(A[t,i], ProbabilityWeights(W[t,i]), corrected=false) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? std(A[t], ProbabilityWeights(W[t]), corrected=false) : NaN
        end
        return result
    end
    export nanstd


    """
    ```julia
    nanmedian(A; dim=0)
    ```
    Calculate the median, ignoring NaNs, of an indexable collection `A`,
    optionally along a dimension specified by `dim`.
    """
    function nanmedian(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef, s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? median(A[i,t]) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? median(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? median(A[t]) : NaN
        end
        return result
    end

    """
    ```julia
    nanmedian(x::AbstractArray{<:Number}, y::AbstractArray{<:Number},
        \txmin::Number, xmax::Number, nbins::Integer)
    ```
    Calculate the median, ignoring NaNs, of y values that fall into each of
    `nbins` equally spaced bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin:(xmax-xmin)/nbins:xmax`
    """
    function nanmedian(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}, xmin::Number, xmax::Number, nbins::Integer)
        binedges = linsp(xmin,xmax,nbins+1)
        medians = Array{Float64}(undef,nbins)
        for i = 1:nbins
            t = (x.>binedges[i]) .& (x.<=binedges[i+1]) .& (.~isnan.(y))
            medians[i] = any(t) ? median(y[t]) : NaN
        end

        return medians
    end
    export nanmedian


    """
    ```julia
    nanmad(A; dim=0)
    ```
    Median absolute deviation from the median, ignoring NaNs, of an indexable
    collection `A`, optionally along a dimension specified by `dim`.
    Note that for a Normal distribution, sigma = 1.4826 * MAD
    """
    function nanmad(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef, s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? median(abs.( A[i,t] .- median(A[i,t]) )) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? median(abs.( A[t,i] .- median(A[t,i]) )) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? median(abs.( A[t] .- median(A[t]) )) : NaN
        end
        return result
    end
    export nanmad


    """
    ```julia
    nanaad(A; dim=0)
    ```
    Mean (average) absolute deviation from the mean, ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dim`.
    Note that for a Normal distribution, sigma = 1.253 * AAD
    """
    function nanaad(A; dim=0)
        s = size(A)
        if dim == 2
            result = Array{eltype(A)}(undef, s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? mean(abs.( A[i,t] .- mean(A[i,t]) )) : NaN
            end
        elseif dim == 1
            result = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? mean(abs.( A[t,i] .- mean(A[t,i]) )) : NaN
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? mean(abs.( A[t] .- mean(A[t]) )) : NaN
        end
        return result
    end
    export nanaad


## --- Array construction

    # Construct linearly spaced array with n points between l and u
    # (linspace replacement )
    if VERSION>=v"0.7"
        function linsp(l::Number,u::Number,n::Number)
            return range(l,stop=u,length=n)
        end
    else
        function linsp(l::Number,u::Number,n::Number)
            return linspace(l,u,n)
        end
    end
    export linsp

## --- Interpolating arrays

    """
    ```julia
    cntr(edges::AbstractArray{<:Number})
    ```
    Given an array of bin edges, return a corresponding vector of bin centers
    """
    function cntr(edges::AbstractArray{<:Number})
        centers = (edges[1:end-1] + edges[2:end]) ./ 2
        return centers
    end
    export cntr

    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, outboundsval=NaN)
        if i > 1 && i < length(y)
            i_below = floor(Int, i)
            i_above = i_below + 1
            f = i - i_below
            return @inbounds Float64(f*y[i_above] + (1-f)*y[i_below])
        else
            return Float64(outboundsval)
        end
    end
    export linterp_at_index

    # Interpolate y-value at xq
    # Linear interpolation, sorting inputs
    if VERSION>v"0.7"
        function linterp1(x,y,xq)
            itp = LinearInterpolation(x,y, extrapolation_bc = Line())
            yq = itp(xq) # Interpolate value of y at queried x values
            return yq
        end
    else
        function linterp1(x,y,xq)
            itp = interpolate((x,),y, Gridded(Linear()))
            yq = itp[xq] # Interpolate value of y at queried x values
            return yq
        end
    end
    export linterp1

    # Sort x and interpolate y-value at xq
    if VERSION>v"0.7"
        function linterp1s(x,y,xq)
            sI = sortperm(x) # indices to construct sorted array
            itp = LinearInterpolation(x[sI], y[sI], extrapolation_bc = Line())
            yq = itp(xq) # Interpolate value of y at queried x values
            return yq
        end
    else
        function linterp1s(x,y,xq)
            sI = sortperm(x) # indices to construct sorted array
            itp = interpolate((x[sI],), y[sI], Gridded(Linear()))
            yq = itp[xq] # Interpolate value of y at queried x values
            return yq
        end
    end
    export linterp1s

    """
    ```julia
    movmean(x::AbstractArray, n::Number)
    ```
    Simple moving average of `x` in 1 or 2 dimensions, spanning `n` bins (or n*n in 2D)
    """
    function movmean(x::AbstractArray, n::Number)
        halfspan = ceil((n-1)/2)
        m = Array{Float64}(undef,size(x))

        # 2-D case
        if length(size(x)) == 2
            iind = repmat(1:size(x,1), 1, size(x,2))
            jind = repmat((1:size(x,2))', size(x,1), 1)
            for k = 1:length(x)
                i = iind[k]
                j = jind[k]
                t = (iind .>= (i-halfspan)) .& (iind .<= (i+halfspan)) .& (jind .>= (j-halfspan)) .& (jind .<= (j+halfspan))
                m[i,j] = mean(x[t])
            end
        # Treat all others as 1-D
        else
            ind = 1:length(x)
            for i in ind
                t = (ind .>= ceil(i-halfspan)) .& (ind .<= ceil(i+halfspan))
                m[i] = mean(x[t])
            end
        end
        return m
    end
    export movmean

## --- Searching arrays

    """
    ```julia
    findmatches(source, target)
    ```
    Return the index of the first value in `target` (if any) that is equal to
    a given value in `source` for each value in `source`; else 0.
    """
    function findmatches(source, target)
        # Allocate output array, initializing with zeros
        index = fill(0, size(source))
        # Allocate match-test array
        t = Array{Bool}(undef,length(target))
        # Loop through source and find first match for each (if any)
        @inbounds for i = 1:length(index)
            for j = 1:length(target)
                t[j] = isequal(source[i], target[j])
            end
            if any(t)
                index[i] = findfirst(t)
            end
        end
        return index
    end
    export findmatches

    """
    ```julia
    findclosest(source, target)
    ```
    Return the index of the numerically closest value in the indexable collection
    `target` for each value in `source`.
    If muliple values are equally close, the first one is used
    """
    function findclosest(source, target)
        # Allocate index and difference arrays
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        # Find closest (numerical) match in target for each value in source
        @inbounds for i = 1:length(source)
            for j = 1:length(diff)
                diff[j] = abs(target[j] - source[i])
            end
            index[i] = argmin(diff)
        end
        return index
    end
    export findclosest

    """
    ```julia
    findclosestbelow(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is less than (i.e., "below") each value in `source`.
    If no such target values exist in `target`, returns an index of 0.
    """
    function findclosestbelow(source, target)
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        @inbounds for i = 1:length(source)
            j = 0
            closestbelow = 0
            while j < length(diff)
                j += 1
                if target[j] < source[i]
                    diff[j] = source[i] - target[j]
                    closestbelow = j
                    break
                end
            end
            while j < length(diff)
                j += 1
                if target[j] < source[i]
                    diff[j] = source[i] - target[j]
                    if diff[j] < diff[closestbelow]
                        closestbelow = j
                    end
                end
            end
            index[i] = closestbelow
        end
        return index
    end
    export findclosestbelow

    """
    ```julia
    findclosestabove(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is greater than (i.e., "above") each value in `source`.
    If no such values exist in `target`, returns an index of 0.
    """
    function findclosestabove(source, target)
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        @inbounds for i = 1:length(source)
            j = 0
            closestabove = 0
            while j < length(diff)
                j += 1
                if target[j] > source[i]
                    diff[j] = target[j] - source[i]
                    closestabove = j
                    break
                end
            end
            while j < length(diff)
                j += 1
                if target[j] > source[i]
                    diff[j] = target[j] - source[i]
                    if diff[j] < diff[closestabove]
                        closestabove = j
                    end
                end
            end
            index[i] = closestabove
        end
        return index
    end
    export findclosestabove

    """
    ```julia
    findnth(t::AbstractArray{Bool}, n::Integer)
    ```
    Return the index of the `n`th true value of `t`, else length(`t`)
    """
    function findnth(t::AbstractArray{Bool}, n::Integer)
        N = 0
        @inbounds for i=1:length(t)
            if t[i]
                N += 1
            end
            if N == n
                return i
            end
        end
        return length(t)
    end
    export findnth


## --- Drawing a pseudorandom array from a numerically specified distribution

    """
    ```julia
    x = draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
    ```
    Draw `n` random floating point numbers from a continuous probability distribution
    specified by a vector `dist` defining the PDF curve thereof.
    """
    function draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
        # Draw n random floating-point numbers from the distribution 'dist'
        x = Array{eltype(dist)}(undef, n)
        dist_ymax = maximum(dist)
        dist_xmax = length(dist) - 1.0

        @inbounds for i = 1:n
            while true
                # Pick random x value
                rx = rand(Float64) * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand(Float64) * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
        return x
    end
    export draw_from_distribution

    """
    ```julia
    draw_from_distribution!(dist::AbstractArray{<:AbstractFloat}, x::Array{<:AbstractFloat})
    ```
    Fill an existing variable `x` with random floating point numbers drawn from
    a continuous probability distribution specified by a vector `dist`
    defining the PDF curve thereof.
    """
    function draw_from_distribution!(dist::AbstractArray{<:AbstractFloat}, x::Array{<:AbstractFloat})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist)
        dist_xmax = length(dist) - 1.0

        @inbounds for i=1:length(x)
            while true
                # Pick random x value
                rx = rand(Float64) * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand(Float64) * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
    end
    export fill_from_distribution

## --- End of File
