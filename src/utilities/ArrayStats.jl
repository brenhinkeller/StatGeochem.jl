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

## --- Transformations of arrays with NaNs

    """
    ```julia
    nanmask(A)
    ```
    Create a Boolean mask of dimensions `size(A)` that is false wherever `A` is `NaN`
    """
    function nanmask(A)
        mask = Array{Bool}(undef,size(A))
        @avx for i=1:length(A)
            mask[i] = !isnan(A[i])
        end
        return mask
    end
    export nanmask

    """
    ```julia
    zeronan!(A)
    ```
    Replace all `NaN`s in A with zeros of the same type
    """
    function zeronan!(A)
        @inbounds @simd for i=1:length(A)
            A[i] *= !isnan(A[i])
        end
        return A
    end
    export zeronan!


## --- Min & max ignoring NaNs

    """
    ```julia
    nanmax(a,b)
    ```
    As `max(a,b)`, but if either argument is `NaN`, return the other one
    """
    function nanmax(a::AbstractFloat,b::AbstractFloat)
        ifelse(isnan(a), b, ifelse(a < b, b, a))
    end
    nanmax(a::AbstractFloat,b::Number) = nanmax(promote(a,b)...)
    nanmax(a::Number,b::AbstractFloat) = nanmax(promote(a,b)...)
    nanmax(a,b) = max(a,b) # Fallback method for non-Floats

    """
    ```julia
    nanmin(a,b)
    ```
    As `min(a,b)`, but if either argument is `NaN`, return the other one
    """
    function nanmin(a::AbstractFloat,b::AbstractFloat)
        ifelse(isnan(a), b, ifelse(a > b, b, a))
    end
    nanmin(a::AbstractFloat,b::Number) = nanmin(promote(a,b)...)
    nanmin(a::Number,b::AbstractFloat) = nanmin(promote(a,b)...)
    nanmin(a,b) = min(a,b) # Fallback method for non-Floats


## --- Percentile statistics, excluding NaNs

    """
    ```julia
    pctile(A, p; dims)
    ```
    Find the `p`th percentile of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    pctile(A, p; dims=:, dim=:) = _pctile(A, p, dims, dim)
    _pctile(A, p, region, ::Colon) = _pctile(A, p, region)
    _pctile(A, p, ::Colon, region) = _pctile(A, p, region) |> vec
    _pctile(A, p, ::Colon, ::Colon) = _pctile(A, p, :)
    function _pctile(A, p, ::Colon)
        t = .~ isnan.(A)
        return any(t) ? percentile(A[t],p) : NaN
    end
    function _pctile(A, p, region)
        s = size(A)
        if region == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif region == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            result =  _pctile(A, p, :)
        end
        return result
    end
    export pctile


    """
    ```julia
    inpctile(A, p::Number; dims)
    ```
    Return a boolean array that identifies which values of the iterable
    collection `A` fall within the central `p`th percentile, optionally along a
    dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    function inpctile(A, p)
        offset = (100 - p) / 2
        return _pctile(A, offset, :) .< A .< _pctile(A, 100-offset, :)
    end
    export inpctile


## --- Summary statistics of arrays with NaNs

    """
    ```julia
    nansum(A; dims)
    ```
    Calculate the sum of an indexable collection `A`, ignoring NaNs, optionally
    along dimensions specified by `dims`.
    """
    nansum(A; dims=:, dim=:) = _nansum(A, dims, dim)
    _nansum(A, region, ::Colon) = _nansum(A, region)
    _nansum(A, ::Colon, region) = _nansum(A, region) |> vec
    _nansum(A, ::Colon, ::Colon) = _nansum(A, :)
    function _nansum(A,::Colon)
        m = zero(eltype(A))
        @simd for x in A
            m += x * !isnan(x)
        end
        return m
    end
    function _nansum(A, region)
        sum(A.*nanmask(A), dims=region)
    end
    export nansum

    """
    ```julia
    nanminimum(A; dims)
    ```
    As `minimum` but ignoring `NaN`s: Find the smallest non-`NaN` value of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanminimum(A; dims=:, dim=:) = _nanminimum(A, dims, dim)
    _nanminimum(A, region, ::Colon) = _nanminimum(A, region)
    _nanminimum(A, ::Colon, region) = _nanminimum(A, region) |> vec
    _nanminimum(A, ::Colon, ::Colon) = _nanminimum(A, :)
    function _nanminimum(A, ::Colon)
        result = A[1]
        for x in A
            result = nanmin(x, result)
        end
        return result
    end
    _nanminimum(A, region) = reduce(nanmin, A, dims=region, init=float(eltype(A))(NaN))
    export nanminimum


    """
    ```julia
    nanmaximum(A; dims)
    ```
    Find the largest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dims`.
    """
    nanmaximum(A; dims=:, dim=:) = _nanmaximum(A, dims, dim)
    _nanmaximum(A, region, ::Colon) = _nanmaximum(A, region)
    _nanmaximum(A, ::Colon, region) = _nanmaximum(A, region) |> vec
    _nanmaximum(A, ::Colon, ::Colon) = _nanmaximum(A, :)
    function _nanmaximum(A, ::Colon)
        result = A[1]
        for x in A
            result = nanmax(x,result)
        end
        return result
    end
    _nanmaximum(A, region) = reduce(nanmax, A, dims=region, init=float(eltype(A))(NaN))
    export nanmaximum


    """
    ```julia
    nanextrema(A; dims)
    ```
    Find the extrema (maximum & minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.
    """
    nanextrema(A; dims=:) = _nanextrema(A, dims)
    _nanextrema(A, ::Colon) = (_nanminimum(A, :), _nanmaximum(A, :))
    _nanextrema(A, region) = collect(zip(_nanminimum(A, region), _nanmaximum(A, region)))
    export nanextrema


    """
    ```julia
    nanrange(A; dims)
    ```
    Calculate the range (maximum - minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.
    """
    nanrange(A; dims=:) = _nanmaximum(A, dims) - _nanminimum(A, dims)
    export nanrange


    """
    ```julia
    nanmean(A, [W]; dims)
    ```
    Ignoring NaNs, calculate the mean (optionally weighted) of an indexable
    collection `A`, optionally along dimensions specified by `dims`.
    """
    nanmean(A; dims=:, dim=:) = _nanmean(A, dims, dim)
    nanmean(A, W; dims=:, dim=:) = _nanmean(A, W, dims, dim)
    function _nanmean(A, ::Colon, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += t
            m += A[i] * t
        end
        return m / n
    end
    function _nanmean(A, region, ::Colon)
        mask = nanmask(A)
        return sum(A.*mask, dims=region) ./ sum(mask, dims=region)
    end
    _nanmean(A, ::Colon, region) = vec(_nanmean(A, region, :))
    function _nanmean(A, W, ::Colon, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += W[i] * t
            m += A[i] * W[i] * t
        end
        return m / n
    end
    function _nanmean(A, W, region, ::Colon)
        mask = nanmask(A)
        return sum(A.*W.*mask, dims=region) ./ sum(W.*mask, dims=region)
    end
    _nanmean(A, W, ::Colon, region) = vec(_nanmean(A, W, region, :))


    """
    ```julia
    nanmean(x, y, [w], xmin::Number, xmax::Number, nbins::Integer)
    ```
    Ignoring NaNs, calculate the mean (optionally weighted) of `y` values that
    fall into each of `nbins` equally spaced `x` bins between `xmin` and `xmax`,
    aligned with bin edges as `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat). If `y` is a 2-d array, then each column of `y` will be
    treated as a separate variable.
    """
    function nanmean(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        N = Array{Int}(undef, nbins, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmean!(MU, N, x, y, xmin, xmax, nbins)
    end
    function nanmean(x::AbstractVector, y::AbstractVecOrMat, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        W = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        MU = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmean!(MU, W, x, y, w, xmin, xmax, nbins)
    end
    export nanmean


    """
    ```julia
    nanmean!(MU, [N], x, y, [w], xmin::Number, xmax::Number, nbins::Integer)
    ```
    Ignoring NaNs, fill the array `MU` with the means (and optionally `N` with
    the counts) of non-NAN `y` values that fall into each of `nbins` equally
    spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If an optional array of weights [`w`] is specified, then `N` is required, and
    will be filled with the sum of weights for each bin.

    The array of `x` data should given as a one-dimensional array (any subtype
    of AbstractVector) and `y` as either a 1-d or 2-d array (any subtype of
    AbstractVecOrMat).

    The output arrays `MU` and `N` must be the same size, and must have the same
    number of columns as `y`; if `y` is a 2-d array (matrix), then each column of
    `y` will be treated as a separate variable.
    """
    function nanmean!(MU::AbstractVecOrMat, x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        N = Array{Int}(undef, size(MU))
        return nanmean!(MU, N, x, y, xmin, xmax, nbins)
    end
    # As above, but also return an array of counts, N; y, N, and MU as 1D vectors
    function nanmean!(MU::AbstractVector, N::AbstractVector, x::AbstractVector, y::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins) && !isnan(y[i])
                bin_index = ceil(Int, bin_index_float)
                N[bin_index] += 1
                MU[bin_index] += y[i]
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning; as above but with y, N, MU as 2D matrices instead of 1D vectors
    function nanmean!(MU::AbstractMatrix, N::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(N, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        N[bin_index,j] += 1
                        MU[bin_index,j] += y[i,j]
                    end
                end
            end
        end
        MU ./= N # Divide by N to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 1D vectors
    function nanmean!(MU::AbstractVector, W::AbstractVector, x::AbstractVector, y::AbstractVector, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins) && !isnan(y[i])
                bin_index = ceil(Int, bin_index_float)
                W[bin_index] += w[i]
                MU[bin_index] += y[i]*w[i]
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    # In-place binning with weights; y, W, MU as 2D matrices
    function nanmean!(MU::AbstractMatrix, W::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, w::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        # Calculate bin index from x value
        scalefactor = nbins / (xmax - xmin)

        # Calculate the means for each bin, ignoring NaNs
        fill!(W, 0)
        fill!(MU, 0) # Fill the output array with zeros to start
        @inbounds for i = 1:length(x)
            bin_index_float = (x[i] - xmin) * scalefactor
            if (0 < bin_index_float < nbins)
                bin_index = ceil(Int, bin_index_float)
                for j = 1:size(y,2)
                    if !isnan(y[i,j])
                        W[bin_index,j] += w[i]
                        MU[bin_index,j] += y[i,j]*w[i]
                    end
                end
            end
        end
        MU ./= W # Divide by sum of weights to calculate means. Empty bin = 0/0 = NaN

        return MU
    end
    export nanmean!


    """
    ```julia
    nanstd(A, [W]; dims)
    ```
    Calculate the standard deviation (optionaly weighted), ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanstd(A; dims=:, dim=:) = _nanstd(A, dims, dim)
    nanstd(A, W; dims=:, dim=:) = _nanstd(A, W, dims, dim)
    function _nanstd(A, ::Colon, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += t
            m += A[i] * t
        end
        mu = m / n
        s = zero(typeof(mu))
        @inbounds @simd for i=1:length(A)
            d = (A[i] - mu) * !isnan(A[i])
            s += d * d
        end
        return sqrt(s / (n-1))
    end
    function _nanstd(A, region, ::Colon)
        mask = nanmask(A)
        N = sum(mask, dims=region)
        s = sum(A.*mask, dims=region)./N
        d = A .- s # Subtract mean, using broadcasting
        @inbounds @simd for i = 1:length(d)
            d[i] = (d[i] * d[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        return @avx @. sqrt( s / (N - 1) )
    end
    _nanstd(A, ::Colon, region) = vec(_nanstd(A, region, :))
    function _nanstd(A, W, ::Colon, ::Colon)
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            w += W[i] * t
            m += A[i] * W[i] * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @inbounds @simd for i=1:length(A)
            d = (A[i] - mu)
            s += d * d * W[i] * !isnan(A[i])
        end
        return sqrt(s / w)
    end
    function _nanstd(A, W, region, ::Colon)
        mask = nanmask(A)
        w = sum(W.*mask, dims=region)
        s = sum(A.*W.*mask, dims=region) ./ w
        d = A .- s # Subtract mean, using broadcasting
        @inbounds @simd for i = 1:length(d)
            d[i] = (d[i] * d[i] * W[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        return @avx @. sqrt( s / w )
    end
    _nanstd(A, W, ::Colon, region) = vec(_nanstd(A, W, region, :))
    export nanstd


    """
    ```julia
    nanmedian(A; dims)
    ```
    Calculate the median, ignoring NaNs, of an indexable collection `A`,
    optionally along a dimension specified by `dims`.
    """
    nanmedian(A; dims=:, dim=:) = _nanmedian(A, dims, dim)
    _nanmedian(A, region, ::Colon) = _nanmedian(A, region)
    _nanmedian(A, ::Colon, region) = _nanmedian(A, region) |> vec
    _nanmedian(A, ::Colon, ::Colon) = _nanmedian(A, :)
    function _nanmedian(A, ::Colon)
        t = .~ isnan.(A)
        return any(t) ? median(A[t]) : float(eltype(A))(NaN)
    end
    function _nanmedian(A, region)
        s = size(A)
        if region == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? median(A[i,t]) : float(eltype(A))(NaN)
            end
        elseif region == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? median(A[t,i]) : float(eltype(A))(NaN)
            end
        else
            result = _nanmedian(A, :)
        end
        return result
    end

    """
    ```julia
    nanmedian(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
    ```
    Calculate the median, ignoring NaNs, of y values that fall into each of `nbins`
    equally spaced `x` bins between `xmin` and `xmax`, aligned with bin edges as
    `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanmedian(x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
        M = Array{float(eltype(y))}(undef, nbins, size(y)[2:end]...)
        return nanmedian!(M, x, y, xmin, xmax, nbins)
    end
    export nanmedian

    """
    ```julia
    nanmedian!(M::AbstractVecOrMat, x::AbstractVector, y::AbstractVecOrMat, xmin::Number, xmax::Number, nbins::Integer)
    ```
    Fill the array `M` with the medians of non-NaN `y` values that fall into
    each of `nbins` equally spaced `x` bins between `xmin` and `xmax`, aligned
    with bin edges as `xmin`:(`xmax`-`xmin`)/`nbins`:`xmax`

    If `y` is a 2-d array (matrix), each column will be treated as a separate variable
    """
    function nanmedian!(M::AbstractVector, x::AbstractVector, y::AbstractVector, xmin::Number, xmax::Number, nbins::Integer)
        binedges = linsp(xmin,xmax,nbins+1)
        t = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= (x.>binedges[i]) .& (x.<=binedges[i+1]) .& (.~isnan.(y))
            M[i] = any(t) ? median(y[t]) : float(eltype(A))(NaN)
        end
        return M
    end
    function nanmedian!(M::AbstractMatrix, x::AbstractVector, y::AbstractMatrix, xmin::Number, xmax::Number, nbins::Integer)
        binedges = linsp(xmin,xmax,nbins+1)
        t = Array{Bool}(undef, length(x))
        tj = Array{Bool}(undef, length(x))
        for i = 1:nbins
            t .= (x.>binedges[i]) .& (x.<=binedges[i+1])
            for j = 1:size(y,2)
                tj .= t .& (.~isnan.(y[:,j]))
                M[i,j] = any(tj) ? median(y[tj,j]) : float(eltype(A))(NaN)
            end
        end
        return M
    end
    export nanmedian!


    """
    ```julia
    nanmad(A; dims)
    ```
    Median absolute deviation from the median, ignoring NaNs, of an indexable
    collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.4826 * MAD
    """
    function nanmad(A; dims=:)
        s = size(A)
        if dims == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? median(abs.( A[i,t] .- median(A[i,t]) )) : float(eltype(A))(NaN)
            end
        elseif dims == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? median(abs.( A[t,i] .- median(A[t,i]) )) : float(eltype(A))(NaN)
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? median(abs.( A[t] .- median(A[t]) )) : float(eltype(A))(NaN)
        end
        return result
    end
    export nanmad


    """
    ```julia
    nanaad(A; dims)
    ```
    Mean (average) absolute deviation from the mean, ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    Note that for a Normal distribution, sigma = 1.253 * AAD
    """
    function nanaad(A; dims=:)
        s = size(A)
        if dims == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? mean(abs.( A[i,t] .- mean(A[i,t]) )) : float(eltype(A))(NaN)
            end
        elseif dims == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? mean(abs.( A[t,i] .- mean(A[t,i]) )) : float(eltype(A))(NaN)
            end
        else
            t = .~ isnan.(A)
            result = any(t) ? mean(abs.( A[t] .- mean(A[t]) )) : float(eltype(A))(NaN)
        end
        return result
    end
    export nanaad



## --- Normalize and standardize arrays

    """
    ```julia
    standardize!(A; dims)
    ```
    Rescale `A` to unit variance and zero mean
    """
    function standardize!(A::Array{<:AbstractFloat}, dims=:)
        A .-= _nanmean(A, dims, :)
        A ./= _nanstd(A, dims, :)
        return A
    end
    export standardize!


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
    cntr(edges::AbstractArray)
    ```
    Given an array of bin edges, return a corresponding vector of bin centers
    """
    function cntr(edges::AbstractArray)
        centers = (edges[1:end-1] + edges[2:end]) ./ 2
        return centers
    end
    export cntr

    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, outboundsval=float(eltype(y))(NaN))
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
        m = Array{float(eltype(x))}(undef,size(x))

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
