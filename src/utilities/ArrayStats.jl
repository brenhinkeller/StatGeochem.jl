## --- Weighted mean of an array

    # Calculate a weigted mean, including MSWD, but without MSWD correction to uncertainty
    function awmean(x, sigma)
        n = length(x)

        if n == 1
            wx = x[1]
            mswd = NaN
            wsigma = sigma[1]
        else
            s1 = 0.0; s2 = 0.0; s3 = 0.0;
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i])
                s2 += 1 / (sigma[i]*sigma[i])
            end
            wx = s1/s2

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i])
            end
            mswd = s3 / (n-1)
            wsigma = sqrt(1.0/s2)
        end
        return wx, wsigma, mswd
    end
    export awmean

    function gwmean(x, sigma)
        # Geochronologist's weigted mean, including MSWD, with MSWD correction to uncertainty.

        n = length(x)

        if n == 1
            wx = x[1]
            mswd = NaN
            wsigma = sigma[1]
        else
            s1 = 0.0; s2 = 0.0; s3 = 0.0;
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i])
                s2 += 1 / (sigma[i]*sigma[i])
            end
            wx = s1/s2

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i])
            end
            mswd = s3 / (n-1)
            wsigma = sqrt(mswd/s2)
        end
        return wx, wsigma, mswd
    end
    export gwmean

    # Calculate MSWD of a dataset
    function MSWD(x, sigma)

        n = length(x)

        s1 = 0.0; s2 = 0.0; s3 = 0.0;
        for i=1:n
            s1 += x[i] / (sigma[i]*sigma[i])
            s2 += 1 / (sigma[i]*sigma[i])
        end
        wx = s1/s2

        for i=1:n
            s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i])
        end

        return s3 / (n-1)
    end
    export MSWD

## --- Percentile statistics, excluding NaNs

    # Percentile of an array along a specified dimension, ignoring NaNs
    function pctile(A,p; dim=0)
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

    # Return a boolean mask for samples within the central nth percentile, optionally along a specified dimension
    function inpctile(A,p; dim=0)
        offset = (100 - p) / 2
        return (A .> pctile(A, offset, dim=dim)) .& (A .< pctile(A, 100-offset, dim=dim))
    end
    export inpctile

## --- Summary statistics of arrays with NaNs

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


    # Smallest non-NaN value of an array
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


    # Largest non-NaN value of an array
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


    # Extrema of an array, ignoring NaNs
    function nanextrema(A)
        t = .~ isnan.(A)
        return extrema(A[t])
    end
    export nanextrema


    # Range (max-min) of an array, ignoring NaNs
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


    # Mean, ignoring NaNs
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
    export nanmean


    # Standard deviation, ignoring NaNs
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
    export nanstd


    # Median, ignoring NaNs
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
    export nanmedian


    # Median absolute deviation from the median, ignoring NaNs
    # For a Normal distribution, sigma = 1.4826 * MAD
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


    # Mean (average) absolute deviation from the mean, ignoring NaNs
    # For a Normal distribution, sigma = 1.253 * AAD
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

    # Return a vector of bin centers if given a vector of bin edges
    function cntr(edges)
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


    function movmean(x, n::Number)
        # Simple moving average
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

    # Find the index of the first value in 'target' (if any) that matches
    # a given value in 'source' for each value in 'source' (else 0)
    function findmatches(source, target)
        # Allocate output array, initializing with zeros
        index = fill(0, size(source))
        # Loop through source and find first match for each (if any)
        for i = 1:length(source)
            t = source[i] .== target
            if any(t)
                index[i] = findfirst(t)
            end
        end
        return index
    end
    export findmatches

    # Return the index of the closest value in 'target' for each value in 'source'
    # If muliple values are equally close, the first one is used
    function findclosest(source, target)
        index=Array{Int64}(undef,size(source))
        for i = 1:length(source)
            index[i] = argmin((target .- source[i]).^2)
        end
        return index
    end
    export findclosest

    # Return the index of the closest value of array 'target' below (less than)
    # each value in 'source'
    function findclosestbelow(source, target)
        index=Array{Int64}(undef, size(source))
        for i = 1:length(source)
            t = findall(target .< source[i])
            ti = argmin((target[t] .- source[i]).^2)
            index[i] = t[ti]
        end
        return index
    end
    export findclosestbelow

    # Return the index of the closest value of the vector 'target' above (greater
    # than) each value in 'source'
    function findclosestabove(source, target)
        index=Array{Int64}(undef, size(source))
        for i=1:length(source)
            t = findall(target .> source[i])
            ti = argmin((target[t] .- source[i]).^2)
            index[i] = t[ti]
        end
        return index
    end
    export findclosestabove

## --- Drawing a pseudorandom array from a numerically specified distribution

    # Draw random numbers from a distribution specified by a vector of points
    # defining the PDF curve
    function draw_from_distribution(dist::Array{<:AbstractFloat}, n::Integer)
        # Draw n random floating-point numbers from the distribution 'dist'
        x = Array{eltype(dist)}(undef, n)
        dist_ymax = maximum(dist)
        dist_xmax = length(dist) - 1.0

        for i = 1:n
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

    # Fill an existing variable with random numbers from a distribution specified
    # by a vector of points defining the PDF curve
    function fill_from_distribution(dist::Array{<:AbstractFloat}, x::Array{<:AbstractFloat})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist)
        dist_xmax = length(dist) - 1.0

        for i=1:length(x)
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
