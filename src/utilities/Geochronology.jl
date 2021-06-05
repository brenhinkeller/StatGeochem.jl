## --- Hafnium isotopes

    # Calculate the initial Hf ratio and epsilon Hf at time t Ma
    function eHf(Hf176_Hf177, Lu176_Hf177, t; eHfOnly::Bool=true)

        # Lutetium decay constant (Soderlund et al., 2004
        lambda = 1.867E-11

        # Present-day CHUR composition (Bouvier et al., 2008)
        CHUR_Hf176_Hf177 = 0.282785
        CHUR_Lu176_Hf177 = 0.0336

        # Calculate initial Hf ratio at time t
        Hf176_Hf177_t = Hf176_Hf177 .- Lu176_Hf177.*(exp.(t .* 10^6*lambda) .- 1)

        # Calculate CHUR Hf ratio at time t
        CHUR_Hf176_Hf177_t = CHUR_Hf176_Hf177 .- CHUR_Lu176_Hf177.*(exp.(t .* 10^6*lambda) .- 1)

        # Calculate corresponding epsilon Hf
        eHf=(Hf176_Hf177_t ./ CHUR_Hf176_Hf177_t .- 1) .* 10^4

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

        means = Array{Float64}(undef,nbins,nresamples)
        c = Array{Float64}(undef,nbins)
        for i=1:nresamples
            dbs = bsresample(data,sigma,length(age))
            eHf_resampled = eHf(dbs[:,2], dbs[:,3], dbs[:,4])
            (c,m,s) = binmeans(dbs[:,1], eHf_resampled, min, max, nbins)
            means[:,i] = m
        end

        m = nanmean(means,dim=2)
        el = m - nanpctile(means,2.5,dim=2)
        eu = nanpctile(means,97.5,dim=2) - m

        return (c, m, el, eu)
    end
    export bin_bsr_eHf


## ---
