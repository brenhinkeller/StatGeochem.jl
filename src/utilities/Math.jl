## --- Fast inverse square-root

    # Float64 version
    function inv_sqrt(number::Float64)
        y = copy(number)
        y = Base.sub_int(9.603007803048109e153, Base.lshr_int(y,1)) # Floating point bit-level hacking
        y *= ( 1.5 - (number * 0.5 * y * y )) # Newton's method
        y *= ( 1.5 - (number * 0.5 * y * y ))
        return y
    end

    # Float32 version
    function inv_sqrt(number::Float32)
        y = copy(number)
        y = Base.sub_int(1.321202f19, Base.lshr_int(y,1)) # Floating point bit-level hacking
        y *= ( 1.5 - (number * 0.5 * y * y )) # Newton's method
        return y
    end

    export inv_sqrt

## --- Base-10 version of log1p

    function log10f(x::Number,from::Number=-1)
        return log10(abs(x-from)+1)*sign(x-from)
    end
    export log10f

## --- Gaussian distribution functions

    # Probability density function of the Normal (Gaussian) distribution
    function normpdf(mu::Number,sigma::Number,x::Number)
        return exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (sqrt(2*pi)*sigma)
    end
    export normpdf

    # Cumulative density function of the Normal (Gaussian) distribution
    # Not precise enough for many uses, unfortunately
    function normcdf(mu::Number,sigma::Number,x::Number)
        return 0.5 + erf((x-mu) / (sigma*sqrt(2))) / 2
    end
    export normcdf

    # How far away from the mean (in units of sigma) should we expect proportion
    # F of the samples to fall in a Normal (Gaussian) distribution
    function NormQuantile(F::Number)
        return sqrt(2)*erfinv(2*F-1)
    end
    export NormQuantile


## --- End of File
