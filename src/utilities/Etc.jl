## --- Digitize image data

    # Calculate approximate x and y positions and uncertainties for colored
    # markers in an image
    function digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)
        # (x,dx,y,dy) = digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)

        # Test for approximate equality in color to marker
        t = isapprox.(marker_color,img,atol=atol)

        # Figure out our image dimensions
        xrows = size(t,2)
        yrows = size(t,1)
        xlims = collect(xlims)
        ylims = collect(ylims)

        # Allocate index arrays
        imin = Array{Float64}(round(Int,xrows/2))
        imax = Array{Float64}(round(Int,xrows/2))
        jmin = Array{Float64}(round(Int,xrows/2))
        jmax = Array{Float64}(round(Int,xrows/2))

        # Fill index arrays
        found = false;
        markernumber = 0
        for j=1:xrows
            if sum(t[:,j])>0
                list = find(t[:,j])
                if ~found
                    markernumber += 1;
                    imin[markernumber] = minimum(list)
                    imax[markernumber] = maximum(list)
                    jmin[markernumber] = j;
                    jmax[markernumber] = j;
                else
                    imin[markernumber] = min(imin[markernumber],minimum(list))
                    imax[markernumber] = max(imax[markernumber],maximum(list))
                    jmax[markernumber] = j;
                end
                found = true;
            else
                found = false;
            end
        end

        # Return only the filled indices
        imin = imin[1:markernumber]
        imax = imax[1:markernumber]
        jmin = jmin[1:markernumber]
        jmax = jmax[1:markernumber]

        # Calculate x and y positions from indices
        y = ylims[2] - (imin+imax)/2 * nanrange(ylims) / yrows
        dy = (imax-imin)/2 * nanrange(ylims) / yrows
        x = (jmin+jmax)/2 * nanrange(xlims) / xrows + xlims[1]
        dx = (jmax-jmin)/2 * nanrange(xlims) / xrows

        return (x,dx,y,dy)
    end

## --- Geochemistry

    # Calculate the initial Hf ratio and epsilon Hf at time t Ma
    function eHf(Hf176_Hf177, Lu176_Hf177, t; eHfOnly=true)

        # Lutetium decay constant (Soderlund et al., 2004
        lambda = 1.867E-11

        # Present-day CHUR composition (Bouvier et al., 2008)
        CHUR_Hf176_Hf177 = 0.282785
        CHUR_Lu176_Hf177 = 0.0336

        # Calculate initial Hf ratio at time t
        Hf176_Hf177_t = Hf176_Hf177 - Lu176_Hf177.*(exp.(t*10^6*lambda) - 1)

        # Calculate CHUR Hf ratio at time t
        CHUR_Hf176_Hf177_t = CHUR_Hf176_Hf177 - CHUR_Lu176_Hf177.*(exp.(t .* 10^6 .* lambda) - 1);

        # Calculate corresponding epsilon Hf
        eHf=(Hf176_Hf177_t ./ CHUR_Hf176_Hf177_t - 1) * 10^4;

        if eHfOnly
            return eHf
        else
            return (eHf, Hf176_Hf177_t)
        end
    end
    export eHf

    function bin_bsr_eHf(x,Hf176_Hf177,Lu176_Hf177,age,min,max,nbins,x_sigma,Hf176_Hf177_sigma,Lu176_Hf177_sigma,age_sigma,nresamples)
        data = hcat(x,Hf176_Hf177,Lu176_Hf177,age)
        sigma = hcat(x_sigma,Hf176_Hf177_sigma,Lu176_Hf177_sigma,age_sigma)

        means = Array{Float64}(nbins,nresamples)
        c = Array{Float64}(nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(age))
            eHf_resampled = eHf(dbs[:,2], dbs[:,3], dbs[:,4])
            (c,m,s) = binmeans(dbs[:,1], eHf_resampled, min, max, nbins)
            means[:,i] = m;
        end

        m = nanmean(means,dim=2)
        el = m - pctile(means,2.5,dim=2)
        eu = pctile(means,97.5,dim=2) - m

        return (c, m, el, eu)
    end
    export bin_bsr_eHf

## --- End of File
