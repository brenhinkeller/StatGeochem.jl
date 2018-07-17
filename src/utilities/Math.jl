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
    function norm_quantile(F::Number)
        return sqrt(2)*erfinv(2*F-1)
    end
    export norm_quantile

## --- Geometry

    # Check if a 2D polygon defined by the arrays x, y contains a point
    function inpolygon(x::Array,y::Array,point::Array)

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

        # For first point, previous point is last
        x_here = copy(x[end])
        y_here = copy(y[end])
        # Ensure we are not sitting on a vertex
        if x_here == point[1] && y_here == point[2]
            y_here = nextfloat(float(y_here))
        end

        # Check how many times a line projected right along x-axis from point intersects the polygon
        intersections = 0
        for i=1:length(x)
            # Recycle our vertex
            x_last = copy(x_here)
            y_last = copy(y_here)

            # Get a new vertex
            x_here = x[i]
            y_here = y[i]

            # Ensure we are not sitting on a vertex by infinitessimally moving the vertex
            if x_here == point[1] && y_here == point[2]
                y_here = nextfloat(float(y_here))
            end

            if y_last > point[2] && y_here > point[2]
                # If both ys above point, no intersection
                continue
            elseif y_last < point[2] && y_here < point[2]
                # If both ys below point, no intersection
                continue
            elseif x_last < point[1] && x_here < point[1]
                # If both xs left of point, no intersection
                continue
            elseif x_last > point[1] && x_here > point[1]
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
                y_range = abs(y_here - y_last)
                if y_range > 0
                    y_fromlast = abs(point[2] - y_last)
                    y_fromhere = abs(point[2] - y_here)
                    x_proj = x_last * y_fromlast / y_range + x_here * y_fromhere / y_range
                    if x_proj > point[1]
                        intersections += 1
                    end
                end
            end
        end

        # If number of intersections is odd, point is in the polygon
        return Bool(mod(intersections,2))
    end

## --- End of File
