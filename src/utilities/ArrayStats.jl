## --- Weighted mean of an array

    # Calculate a weigted mean, including MSWD, but without MSWD correction to uncertainty
    function awmean(x, sigma)
        n = length(x);

        if n==1
            wx = x[1];
            mswd = 0;
            wsigma = sigma[1];
        else
            s1 = 0.0; s2 = 0.0; s3 = 0.0;
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i]);
                s2 += 1 / (sigma[i]*sigma[i]);
            end
            wx = s1/s2;

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i]);
            end
            mswd = s3 / (n-1);
            wsigma = sqrt(1.0/s2);
        end
        return wx, wsigma, mswd
    end
    export awmean


    function gwmean(x, sigma)
        # Geochronologist's weigted mean, including MSWD, with MSWD correction to uncertainty.

        n = length(x);

        if n==1
            wx = x[1];
            mswd = 0;
            wsigma = sigma[1];
        else
            s1 = 0.0; s2 = 0.0; s3 = 0.0;
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i]);
                s2 += 1 / (sigma[i]*sigma[i]);
            end
            wx = s1/s2;

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i]);
            end
            mswd = s3 / (n-1);
            wsigma = sqrt(mswd/s2);
        end
        return wx, wsigma, mswd
    end
    export gwmean

## --- Statistics of arrays with NaNs

    # Percentile of an array along a specified dimension, ignoring NaNs
    function pctile(A,p;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = percentile(A[i,t],p);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = percentile(A[t,i],p);
            end
        else
            out = percentile(A[:],p);
        end
        return out
    end
    export pctile

    function nansum(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = sum(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = sum(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = sum(A[t]);
        end
        return out
    end
    export nansum

    # Smallest non-NaN value of an array
    function nanminimum(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = minimum(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = minimum(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = minimum(A[t]);
        end
        return out
    end
    export nanminimum

    # Largest non-NaN value of an array
    function nanmaximum(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = maximum(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = maximum(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = maximum(A[t]);
        end
        return out
    end
    export nanmaximum

    # Range (max-min) of an array, ignoring NaNs
    function nanrange(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                extr = extrema(A[i,t]);
                out[i] = extr[2] - extr[1]
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                extr = extrema(A[t,i]);
                out[i] = extr[2] - extr[1]
            end
        else
            t = .~ isnan.(A)
            extr = extrema(A[t])
            out = extr[2] - extr[1];
        end
        return out
    end
    export nanrange

    # Mean, ignoring NaNs
    function nanmean(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = mean(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = mean(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = mean(A[t]);
        end
        return out
    end
    export nanmean

    # Standard deviation, ignoring NaNs
    function nanstd(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = std(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = std(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = std(A[t]);
        end
        return out
    end
    export nanstd

    # Median, ignoring NaNs
    function nanmedian(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = median(A[i,t]);
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = median(A[t,i]);
            end
        else
            t = .~ isnan.(A)
            out = median(A[t]);
        end
        return out
    end
    export nanmedian

    # Median absolute deviation from the median, ignoring NaNs
    # For a Normal distribution, sigma = 1.4826 * MAD
    function nanmad(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = median(abs.( A[i,t] - median(A[i,t]) ));
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = median(abs.( A[t,i] - median(A[t,i]) ));
            end
        else
            t = .~ isnan.(A)
            out = median(abs.( A[t] - median(A[t]) ));
        end
        return out
    end
    export nanmad

    # Mean (average) absolute deviation from the mean, ignoring NaNs
    # For a Normal distribution, sigma = 1.253 * AAD
    function nanaad(A;dim=0)
        s = size(A);
        if dim==2
            out = Array{typeof(A[1])}(s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = mean(abs.( A[i,t] - mean(A[i,t]) ));
            end
        elseif dim==1
            out = Array{typeof(A[1])}(s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = mean(abs.( A[t,i] - mean(A[t,i]) ));
            end
        else
            t = .~ isnan.(A)
            out = mean(abs.( A[t] - mean(A[t]) ));
        end
        return out
    end
    export nanaad

## --- Interpolating arrays

    # Return a vector of bin centers if given a vector of bin edges
    function cntr(edges)
        centers = (edges[1:end-1]+edges[2:end])/2;
        return centers;
    end
    export cntr

    # Interpolate y-value at xq
    function linterp1(x,y,xq)
        itp = interpolate((x,),y,Gridded(Linear()));
        yq = itp[xq]; # Interpolate value of y at queried x values
        return yq
    end
    export linterp1

    # Sort x and interpolate y-value at xq
    function linterp1s(x,y,xq)
        sI = sortperm(x); # indices to construct sorted array
        itp = interpolate((x[sI],),y[sI],Gridded(Linear()));
        yq = itp[xq]; # Interpolate value of y at queried x values
        return yq
    end
    export linterp1s


    function movmean(x, n::Number)
        # Simple moving average
        halfspan = ceil((n-1)/2);
        m = Array{Float64}(size(x))

        # 2-D case
        if length(size(x)) == 2
            iind = repmat(1:size(x,1),1,size(x,2))
            jind = repmat((1:size(x,2))',size(x,1),1)
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
        return m;
    end
    export movmean

## --- Searching arrays

    # Return the index of the closest value of Target for each value in Source
    function findclosest(source, target)
        index=Array{Int64}(size(source));
        for i=1:length(source)
            index[i]=indmin((target-source[i]).^2);
        end
        return index
    end
    export findclosest

    # Return the index of the closest value of the vector 'target' below each
    # value in 'source'
    function findclosestbelow(source, target)
        index=Array{Int64}(size(source));
        for i=1:length(source)
            t = find(target.<source[i]);
            ti = indmin((target[t]-source[i]).^2);
            index[i] = t[ti];
        end
        return index;
    end
    export findclosestbelow

    # Return the index of the closest value of the vector 'target' above each
    # value in 'source'
    function findclosestabove(source, target)
        index=Array{Int64}(size(source));
        for i=1:length(source)
            t = find(target.>source[i]);
            ti = indmin((target[t]-source[i]).^2);
            index[i] = t[ti];
        end
        return index;
    end
    export findclosestabove

## --- Drawing a pseudorandom array from a numerically specified distribution

    # Draw random numbers from a distribution specified by a vector of points
    # defining the PDF curve
    function draw_from_distribution(dist::Array{Float64}, n::Int)
        # Draw n random numbers from the distribution 'dist'
        x = Array{Float64}(n);
        dist_ymax = maximum(dist);
        dist_xmax = length(dist)-1.0;

        for i=1:n;
            while true
                # Pick random x value
                rx = rand() * dist_xmax;
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax;
                if (y > ry)
                    x[i] = rx / dist_xmax;
                    break
                end
            end
        end
        return x
    end
    export draw_from_distribution

    # Fill an existing variable with random numbers from a distribution specified
    # by a vector of points defining the PDF curve
    function fill_from_distribution(dist::Array{Float64}, x::Array{Float64})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist);
        dist_xmax = length(dist)-1.0;

        for i=1:length(x);
            while true
                # Pick random x value
                rx = rand() * dist_xmax;
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax;
                if (y > ry)
                    x[i] = rx / dist_xmax;
                    break
                end
            end
        end
    end
    export fill_from_distribution

## --- End of File
