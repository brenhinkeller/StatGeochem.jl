## --- Linear least-squares regression

    function linreg(x,y)
        return hcat(fill!(similar(x), 1), x) \ y
    end
    export linreg

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

    # Fast Log Likelihood corresponding to a Normal (Gaussian) distribution
    function normpdf_LL(mu::Number,sigma::Number,x::Number)
        return -(x-mu)*(x-mu) / (2*sigma*sigma)
    end
    export normpdf_LL

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

    # How dispersed (in units of sigma) should we expect a sample of N numbers
    # drawn from a Normal (Gaussian) distribution to be?
    function norm_width(N::Number)
        F = 1 - 1/(N+1)
        return 2*norm_quantile(F)
    end
    export norm_width

    # Integral of the product of two normal distributions N(μ1,σ1) * N(μ2,σ2)
    function normproduct(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        # The integral of the product of two normal distributions is itself just
        # another Normal distribution! Specifically, one with variance σ1^2 + σ2^2
        normpdf(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end
    export normproduct

    # Log likelihood corresponding to the integral of N(μ1,σ1) * N(μ2,σ2)
    function normproduct_LL(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        # As above, but using the fast log likelihood of a Normal distribution
        normpdf_LL(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end
    export normproduct_LL


## --- Geometry

    # Check if a 2D polygon defined by the arrays x, y contains a point.
    # Returns boolean (true or false)
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
    export inpolygon

    # Find the indexes of grid points that fall within a polygon for a grid /
    # matrix with cell centers given by grid_x (j-columns of matrix) and
    # grid_y (i-rows of matrix).
    # Returns a list of rows and columns in the polygon
    function find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)
        # (row, column) = find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)

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


    # Calculate the distance between two (lat,lon) points on a sphere.
    # Lat, Lon in decimal degrees
    function arcdistance(lat_i,lon_i,lat,lon)
        # Argument for acos()
        arg = sin.(lat_i .* pi/180) .* sin.(lat .* pi/180) .+ cos.(lat_i*pi/180) .* cos.(lat .* pi/180).*cos.((lon_i .- lon) .* pi/180)

        # Avoid domain errors from imprecise sine and cosine math
        arg[arg .> 1.0] .= 1.0
        arg[arg .< -1.0] .= -1.0

        # Calculate angular distance
        theta = 180/pi .* acos.(arg)
        return theta
    end
    export arcdistance

## --- Empirical and numerical distributions

    # Two-sided linear exponential distribution joined by an atan sigmoid.
    function bilinear_exponential(x,p)
        # If to a normal-esque PDF, parameters p roughly correspond to:
        # p[1] = pre-exponential (normaliation constant)
        # p[2] = mean (central moment)
        # p[3] = standard deviation
        # p[4] = sharpness
        # p[5] = skew
        xs = (x-p[2])./p[3]; # X scaled by mean and variance
        v = 1/2-atan.(xs)/pi; # Sigmoid (positive on LHS)
        f = p[1] .* exp.((p[4].^2).*(p[5].^2).*xs.*v - (p[4].^2)./(p[5].^2).*xs.*(1-v))
        return f
    end
    export bilinear_exponential

    # Log of two-sided linear exponential distribution joined by an atan sigmoid.
    function bilinear_exponential_LL(x,p)
        # If to a normal-esque PDF, parameters p roughly correspond to:
        # p[1] = pre-exponential (normaliation constant)
        # p[2] = mean (central moment)
        # p[3] = standard deviation
        # p[4] = sharpness
        # p[5] = skew
        xs = (x-p[2,:])./p[3,:]; # X scaled by mean and variance
        v = 1/2-atan.(xs)/pi; # Sigmoid (positive on LHS)
        f = log.(p[1,:]) + (p[4,:].^2).*(p[5,:].^2).*xs.*v - (p[4,:].^2)./(p[5,:].^2).*xs.*(1-v)
        return f
    end
    export bilinear_exponential_LL

    # Interpolate log likelihood from an array
    function interpolate_LL(x,p)
        ll = Array{Float64}(undef, size(x))
        for i = 1:length(x)
            ll[i] = linterp_at_index(p[:,i], x[i], -Inf)
        end
        return ll
    end
    export interpolate_LL

    UniformDistribution = [1.0, 1.0]
    export UniformDistribution

    TriangularDistribution = [2.0,1.0,0.0]
    export TriangularDistribution

    TruncatedNormalDistribution =
    [1.15708, 1.2038, 1.25037, 1.29662, 1.34239, 1.38751, 1.43181, 1.47511, 1.51724, 1.55802, 1.5973, 1.63489, 1.67064, 1.70439, 1.73598, 1.76527, 1.79213, 1.81643, 1.83805, 1.8569, 1.87289, 1.88593, 1.89596, 1.90294, 1.90682, 1.9076, 1.90527, 1.89983, 1.89133, 1.87978, 1.86526, 1.84783, 1.82758, 1.80461, 1.77901, 1.75092, 1.72046, 1.68777, 1.65301, 1.61632, 1.57786, 1.53781, 1.49633, 1.45359, 1.40977, 1.36504, 1.31958, 1.27355, 1.22712, 1.18045, 1.13371, 1.08704, 1.04059, 0.9945, 0.948905, 0.903921, 0.859665, 0.816241, 0.773749, 0.732271, 0.691887, 0.652662, 0.614655, 0.577921, 0.542496, 0.508411, 0.475689, 0.444348, 0.414395, 0.38583, 0.358649, 0.332839, 0.308382, 0.285256, 0.263435, 0.242885, 0.223573, 0.205461, 0.188509, 0.172673, 0.157909, 0.144172, 0.131416, 0.119592, 0.108655, 0.0985572, 0.0892521, 0.0806935, 0.0728368, 0.0656377, 0.0590535, 0.0530436, 0.0475671, 0.0425867, 0.0380655, 0.0339688, 0.0302636, 0.0269185, 0.0239041, 0.0211926]
    export TruncatedNormalDistribution

    HalfNormalDistribution =
    [2.65307, 2.6516, 2.64718, 2.63984, 2.6296, 2.61648, 2.60054, 2.58182, 2.5604, 2.53633, 2.5097, 2.48059, 2.4491, 2.41532, 2.37936, 2.34133, 2.30135, 2.25954, 2.21603, 2.17095, 2.12442, 2.07657, 2.02755, 1.97749, 1.92653, 1.87479, 1.82242, 1.76954, 1.71629, 1.66279, 1.60917, 1.55555, 1.50205, 1.44878, 1.39584, 1.34335, 1.29139, 1.24006, 1.18946, 1.13965, 1.09071, 1.04272, 0.995728, 0.9498, 0.904985, 0.861327, 0.818864, 0.777631, 0.737653, 0.698954, 0.66155, 0.625452, 0.590667, 0.557197, 0.52504, 0.494189, 0.464635, 0.436363, 0.409356, 0.383594, 0.359054, 0.335711, 0.313537, 0.292503,0.272576, 0.253725, 0.235916, 0.219112, 0.20328, 0.188382, 0.174383, 0.161244, 0.14893, 0.137403, 0.126628, 0.116568,0.107188, 0.0984536, 0.0903304, 0.0827854, 0.0757864, 0.069302, 0.0633021, 0.0577574, 0.0526399, 0.0479225, 0.0435794, 0.0395859, 0.0359185, 0.0325546, 0.029473, 0.0266534, 0.0240769, 0.0217252, 0.0195815, 0.0176297, 0.0158548, 0.0142428, 0.0127805, 0.0114555, 0.0102566]
    export HalfNormalDistribution

    EllisDistribution =
    [6.80942, 5.36792, 4.45867, 3.83457, 3.28267, 2.77244, 2.33403, 1.98717, 1.72219, 1.52427, 1.37856, 1.27023, 1.18445,1.10697, 1.03176, 0.958823, 0.888329, 0.820435, 0.755302, 0.693089, 0.633956, 0.578063, 0.525569, 0.476635, 0.431421,0.390065, 0.35249, 0.31849, 0.287855, 0.260378, 0.235851, 0.214064, 0.194811, 0.177882, 0.163069, 0.150163, 0.138957,0.129243, 0.120811, 0.113453, 0.106962, 0.101129, 0.095775, 0.0908504, 0.0863393, 0.0822255, 0.0784932, 0.0751262, 0.0721086, 0.0694244, 0.0670576]
    export EllisDistribution

    MeltsZirconDistribution =
    [0.282361, 0.28919, 0.296019, 0.302849, 0.30968, 0.316567, 0.323614, 0.33064, 0.337727, 0.344848, 0.352146, 0.359642,0.367482, 0.375622, 0.384052, 0.392828, 0.401968, 0.411518, 0.421442, 0.43171, 0.44213, 0.45295, 0.464036, 0.47539, 0.486938, 0.498644, 0.51075, 0.523026, 0.535688, 0.548764, 0.562124, 0.575927, 0.590363, 0.604879, 0.620415, 0.636022, 0.652333, 0.669112, 0.686441, 0.704341, 0.72283, 0.742036, 0.761964, 0.782541, 0.803718, 0.825707, 0.848386, 0.871895,0.896139, 0.920462, 0.946071, 0.972964, 0.999905, 1.02776, 1.05664, 1.08637, 1.11731, 1.14919, 1.18202, 1.21582, 1.24956, 1.28342, 1.31828, 1.35427, 1.39153, 1.43006, 1.46879, 1.50812, 1.5477, 1.58888, 1.63149, 1.6748, 1.71724, 1.76126, 1.80668, 1.85101, 1.89546, 1.94144, 1.98379, 2.02785, 2.06738, 2.10669, 2.1377, 2.16306, 2.17843, 2.17924, 2.16073, 2.11744, 2.04444, 1.93323, 1.7923, 1.62527, 1.44425, 1.25401, 1.0528, 0.843628, 0.632687, 0.421876, 0.211064, 0.000252985]
    export MeltsZirconDistribution

    MeltsVolcanicZirconDistribution =
    [0.54933, 0.556409, 0.563488, 0.570567, 0.577653, 0.584759, 0.591912, 0.599251, 0.606793, 0.614519, 0.622425, 0.630421, 0.63852, 0.646681, 0.654972, 0.663533, 0.672274, 0.681233, 0.690399, 0.699787, 0.709334, 0.719174, 0.729157, 0.739423, 0.749935, 0.760644, 0.771726, 0.782974, 0.794507, 0.806296, 0.818297, 0.830517, 0.842957, 0.855411, 0.866744, 0.878127, 0.889792, 0.901792, 0.914121, 0.926689, 0.939557, 0.952834, 0.966425, 0.980333, 0.994521, 1.00914, 1.02403, 1.03928, 1.05487, 1.0705, 1.08587, 1.10097, 1.11608, 1.13153, 1.1474, 1.16353, 1.18025, 1.19743, 1.21504, 1.23312, 1.25034, 1.26711, 1.28441, 1.30212, 1.32024, 1.33892, 1.35769, 1.37491, 1.3923, 1.41046, 1.42924, 1.44775, 1.46432, 1.48171, 1.49969, 1.51516, 1.53001, 1.54571, 1.5566, 1.56814, 1.57522, 1.58168, 1.58206, 1.57869, 1.56907, 1.55064, 1.51982, 1.4737, 1.40944, 1.32047, 1.21218, 1.09157, 0.965488, 0.834108, 0.697552, 0.558304, 0.418827, 0.279262, 0.139695, 0.000127237]
    export MeltsVolcanicZirconDistribution


## --- End of File
