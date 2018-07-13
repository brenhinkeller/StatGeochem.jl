## --- Bootstrap resampling

    function bsresample(data, sigma, nrows; p=0.5)
        # resampled = bsresample(data, sigma, nrows; p=0.5)
        # Bootstrap resample a variable up to size nrows. Optionally provide weights in p

        resampled = Array{Float64}(nrows,size(data,2))

        i = 1;
        while i <= nrows

            if size(data,1) > 1
                # If we have more than one sample

                # Select weighted sample of data
                t = rand(size(data,1)) .< p
                sdata = data[t,:]

                # Corresponing uncertainty (either blanket or for each datum)
                if size(sigma,1) > 1
                    serr = sigma[t,:]
                else
                    serr = ones(size(sdata)) .* sigma
                end

            else
                # If only one sample
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
    export bsresample

## --- Bin a dataset by a given independent variable

    function binmeans(x,y,min,max,nbins; oversamplingratio=1)
        binwidth=(max-min)/nbins
        binedges=linspace(min,max,nbins+1)
        bincenters=(min+binwidth/2):binwidth:(max-binwidth/2)

        averages=Array{Float64}(nbins)
        errors=Array{Float64}(nbins)
        for i=1:nbins
            t = (x.>binedges[i]) .& (x.<binedges[i+1]) .& (.~isnan.(y))
            averages[i] = mean(y[t])
            errors[i] = std(y[t]) ./ sqrt(sum(t)) .* sqrt(oversamplingratio)
        end

        return(bincenters, averages, errors)
    end
    export binmeans

## --- Bin bootstrap resampled data

    function bin_bsr_means(x,y,min,max,nbins,x_sigma,nresamples)
        data = hcat(x, y)
        sigma = hcat(x_sigma, zeros(size(x_sigma)))

        means = Array{Float64}(nbins,nresamples)
        c = Array{Float64}(nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(x))
            (c,m,s) = binmeans(dbs[:,1], dbs[:,2], min, max, nbins)
            means[:,i] = m;
        end

        m = nanmean(means,dim=2)
        el = m - pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) - m

        return (c, m, el, eu)
    end
    export bin_bsr_means

## --- Downsample an image / array

    function downsample(matrix::Array,factor::Int)
        if ndims(matrix)==2
            rows = floor(Int,size(matrix,1)/factor);
            cols = floor(Int,size(matrix,2)/factor);

            downsampled = typeof(matrix)(rows,cols)
            for i=1:rows
                for j=1:cols
                    downsampled[i,j]=matrix[i*factor,j*factor]
                end
            end
        else
            downsampled = matrix[factor:factor:end]
        end
        return downsampled
    end
    export downsample


## --- End of file
