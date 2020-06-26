## --- Fast inverse square-root
    # This is mostly just for fun: while it can be 1.5x faster than 1/sqrt(x)
    # I'm not sure there's any realistic scientific application where it's
    # worth the loss of precision. Both versions good to about 4 ppm

    # Float64 version
    function inv_sqrt(number::Float64)
        n_2 = 0.5 * number
        y = Base.sub_int(9.603007803048109e153, Base.lshr_int(number,1)) # Floating point magic
        y *= ( 1.5 - (n_2 * y * y )) # Newton's method
        y *= ( 1.5 - (n_2 * y * y )) # Newton's method (again)
        return y
    end

    # Float32 version
    function inv_sqrt(number::Float32)
        n_2 = 0.5f0 * number
        y = Base.sub_int(1.321202f19, Base.lshr_int(number,1)) # Floating point magic
        y *= ( 1.5f0 - (n_2 * y * y) ) # Newton's method
        y *= ( 1.5f0 - (n_2 * y * y) ) # Newton's method (again)
        return y
    end

    export inv_sqrt

## --- Base-10 version of log1p

    function log10f(x::Number,from::Number=-1)
        return log10(abs(x-from)+1)*sign(x-from)
    end
    export log10f

## --- Gaussian distribution functions

    """
    ```julia
    normpdf(mu,sigma,x)
    ```
    Probability density function of the Normal (Gaussian) distribution

    ``ℯ^{-(x-μ)^2 / (2σ^2)} / σ√2π``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`
    """
    function normpdf(mu,sigma,x)
        return @. exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (sqrt(2*pi)*sigma)
    end
    export normpdf

    """
    ```julia
    normpdf_ll(mu, sigma, x)
    ```
    Fast log likelihood corresponding to a Normal (Gaussian) distribution
    with mean `mu` and standard deviation `sigma`, evaluated at `x`.

    If `x`, `mu`, and `sigma` are given as arrays, the sum of the log likelihood
    over all `x` will be returned.

    See also `normpdf`
    """
    function normpdf_ll(mu::Number,sigma::Number,x::Number)
        return -(x-mu)*(x-mu) / (2*sigma*sigma)
    end
    function normpdf_ll(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = 0.0
        @avx for i=1:length(x)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) / (2*sigma[i]*sigma[i])
        end
        return ll
    end

    """
    ```julia
    normcdf(mu,sigma,x)
    ```
    Cumulative density function of the Normal (Gaussian) distribution

    ``1/2 + erf(\frac{x-μ}{σ√2})/2``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`.
    """
    function normcdf(mu,sigma,x)
        return @. 0.5 + erf((x-mu) / (sigma*sqrt(2))) / 2
    end
    function normcdf(mu::Number,sigma::Number,x::AbstractArray)
        result = Array{float(eltype(x))}(undef,length(x))
        sigma_sqrt = sigma*sqrt(2)
        @inbounds @simd for i = 1:length(x)
            result[i] = 0.5 + erf((x[i]-mu) / sigma_sqrt) / 2
        end
        return result
    end
    export normcdf

    """
    ```julia
    norm_quantile(F::Number)
    ```
    How far away from the mean (in units of sigma) should we expect proportion
    F of the samples to fall in a standard Gaussian (Normal[0,1]) distribution
    """
    function norm_quantile(F::Number)
        return sqrt(2)*erfinv(2*F-1)
    end
    export norm_quantile

    """
    ```julia
    norm_width(N::Number)
    ```
    How dispersed (in units of sigma) should we expect a sample of N numbers
    drawn from a standard Gaussian (Normal[0,1]) distribution to be?
    """
    function norm_width(N::Number)
        F = 1 - 1/(N+1)
        return 2*norm_quantile(F)
    end
    export norm_width

    """
    ```julia
    normproduct(μ1, σ1, μ2, σ2)
    ```
    The integral of the product of two normal distributions N[μ1,σ1] * N[μ2,σ2].
    This is itself just another Normal distribution! Specifically, one with
    variance σ1^2 + σ2^2, evaluated at distance |μ1-μ2| from the mean
    """
    function normproduct(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        normpdf(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end
    export normproduct

    """
    ```julia
    normproduct_ll(μ1, σ1, μ2, σ2)
    ```
    Log likelihood corresponding to the integral of N[μ1,σ1] * N[μ2,σ2]
    As `normproduct`, but using the fast log likelihood of a Normal distribution
    """
    function normproduct_ll(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        normpdf_ll(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end
    export normproduct_ll


## --- Geometry

    """
    ```julia
    inpolygon(x,y,point)
    ```
    Check if a 2D polygon defined by the arrays `x`, `y` contains a given `point`.
    Returns boolean (true or false)
    """
    function inpolygon(x,y,point)
        # Check that we have the right kind of input data
        if length(x) != length(y)
            error("polygon must have equal number of x and y points\n")
        end
        if length(x) < 3
            error("polygon must have at least 3 points\n")
        end
        if length(point) != 2
            error("point must be an ordered pair (x,y)\n")
        end

        # Extract x and y data of point
        point_x = point[1]
        point_y = point[2]
        # For first point, previous point is last
        x_here = x[end]
        y_here = y[end]
        # Ensure we are not sitting parallel to a vertex by infinitessimally moving the point
        if y_here == point_y
            point_y = nextfloat(float(point_y))
        end
        if x_here == point_x
            point_x = nextfloat(float(point_x))
        end

        # Check how many times a line projected right along x-axis from point intersects the polygon
        intersections = 0
        @inbounds for i=1:length(x)
            # Recycle our vertex
            x_last = copy(x_here)
            y_last = copy(y_here)

            # Get a new vertex
            x_here = x[i]
            y_here = y[i]

            # Ensure we are not sitting parallel to a vertex by infinitessimally moving the point
            if y_here == point_y
                point_y = nextfloat(float(point_y))
            end
            if x_here == point_x
                point_x = nextfloat(float(point_x))
            end

            if y_last > point_y && y_here > point_y
                # If both ys above point, no intersection
                continue
            elseif y_last < point_y && y_here < point_y
                # If both ys below point, no intersection
                continue
            elseif x_last < point_x && x_here < point_x
                # If both x's left of point, no intersection
                continue
            elseif x_last > point_x && x_here > point_x
                # By elimination
                # We have one y above and y below our point
                # If both y's are right of line, then definite intersection
                intersections += 1
                continue
            else
                # By elimination
                # One y above and one y below
                # One x to the right and one x to the left
                # We must project
                dy = y_here - y_last
                if abs(dy) > 0
                    dx = x_here - x_last
                    inv_slope = dx / dy
                    x_proj = x_last + (point_y - y_last) * inv_slope
                    if x_proj > point_x
                        intersections += 1
                    end
                end
            end
        end

        # If number of intersections is odd, point is in the polygon
        return Bool(mod(intersections,2))
    end
    export inpolygon


    """
    ```julia
    (rows, columns) = find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)
    ```
    Find the indexes of grid points that fall within a polygon for a grid with
    cell centers given by grid_x (j-columns of grid) and grid_y (i-rows of grid).
    Returns a list of rows and columns in the polygon
    """
    function find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)
        # Check that we have the right kind of input data
        if length(poly_x) != length(poly_y)
            error("polygon must have equal number of x and y points\n")
        end
        if length(poly_x) < 3
            error("polygon must have at least 3 points\n")
        end

        # Find maximum x and y range of polygon
        (xmin, xmax) = extrema(poly_x)
        (ymin, ymax) = extrema(poly_y)

        # Find the matrix indices within the range of the polygon (if any)
        column_inrange = findall((grid_x .>= xmin) .& (grid_x .<= xmax))
        row_inrange = findall((grid_y .>= ymin) .& (grid_y .<= ymax))

        # Keep a list of matrix indexes in the polygon
        row = Array{Int}(undef,length(column_inrange) * length(row_inrange))
        column = Array{Int}(undef,length(column_inrange) * length(row_inrange))
        n = 0
        for j = 1:length(column_inrange)
            for i = 1:length(row_inrange)
                point = [grid_x[column_inrange[j]], grid_y[row_inrange[i]]]
                if inpolygon(poly_x, poly_y, point)
                    n += 1
                    row[n] = row_inrange[i]
                    column[n] = column_inrange[j]
                end
            end
        end

        return (row[1:n], column[1:n])
    end
    export find_grid_inpolygon


    """
    ```julia
    arcdistance(latᵢ,lonᵢ,lat,lon)
    ```
    Calculate the distance on a sphere between the point (`latᵢ`,`lonᵢ`) and any
    number of points in (`lat`,`lon`).
    Latitude and Longitude should be specified in decimal degrees
    """
    function arcdistance(latᵢ,lonᵢ,lat,lon)
        # Argument for acos()
        arg = sin.(latᵢ .* pi/180) .* sin.(lat .* pi/180) .+ cos.(latᵢ*pi/180) .* cos.(lat .* pi/180).*cos.((lonᵢ .- lon) .* pi/180)

        # Avoid domain errors from imprecise sine and cosine math
        arg[arg .> 1.0] .= 1.0
        arg[arg .< -1.0] .= -1.0

        # Calculate angular distance
        theta = 180/pi .* acos.(arg)
        return theta
    end
    export arcdistance

## --- Empirical and numerical distributions



## --- End of File
