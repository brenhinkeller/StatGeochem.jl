## --- Read ESRI Arc/Info ASCII grid files

    function importAAIGrid(fname, T=Float64; undefval=NaN)
        # Open the file
        fid = open(fname)

        metadata = Dict{String,Number}()
        metadata["ncols"] = parse(Int64, first(match(r"  *(.*?)$", readline(fid))))
        metadata["nrows"] = parse(Int64, first(match(r"  *(.*?)$", readline(fid))))
        metadata["xll_corner"] = parse(Float64, first(match(r"  *(.*?)$", readline(fid))))
        metadata["yll_corner"] = parse(Float64, first(match(r"  *(.*?)$", readline(fid))))
        metadata["cellsize"] = parse(Float64, first(match(r"  *(.*?)$", readline(fid))))
        metadata["nodata"] = parse(Float64, first(match(r"  *(.*?)$", readline(fid))))

        nrows = metadata["nrows"]
        ncols = metadata["ncols"]

        data = Array{T}(undef,ncols,nrows)
        for i = 1:nrows
            l = readline(fid)
            delim_string_parse!(data, l, ' ', T, offset=(i-1)*ncols, undefval=undefval)
        end

        # Close the file
        close(fid)

        return (data', metadata)
    end
    export importAAIGrid

## --- Calculate slope from a DEM

    function maxslope(matrix, x_lon_cntr, y_lat_cntr, cellsize, T=UInt16; minmatval=-12000, km_per_lat=111.1)
        # Returns slope in units/kilometer given a latitude-longitude grid of z-values

        # Allocate output array
        slope = Array{T}(undef,size(matrix))

        # Fill in the center first
        distNS = 2 * cellsize * km_per_lat
        for i = 2:(size(matrix,1)-1)
            # Distance between grid cell centers
            km_per_lon = cos(y_lat_cntr[i] * pi/180) * km_per_lat
            distEW = 2 * cellsize * km_per_lon
            distDiag = sqrt(distNS^2 + distEW^2)

            for j = 2:(size(matrix,2)-1)
                # Gradients, in matrix units per km
                if matrix[i,j] < minmatval
                    slope[i,j] = 0
                else
                    if (matrix[i+1,j] < minmatval) || (matrix[i-1,j] < minmatval)
                        NS = 0
                    else
                        NS = abs(matrix[i+1,j] - matrix[i-1,j]) / distNS
                    end
                    if (matrix[i,j+1] < minmatval) || (matrix[i,j-1]<minmatval)
                        EW = 0
                    else
                        EW = abs(matrix[i,j+1] - matrix[i,j-1]) / distEW
                    end
                    if (matrix[i+1,j-1] < minmatval) || (matrix[i-1,j+1] < minmatval)
                        NESW = 0
                    else
                        NESW = abs(matrix[i+1,j-1] - matrix[i-1,j+1]) / distDiag
                    end
                    if (matrix[i+1,j+1] < minmatval) || (matrix[i-1,j-1] < minmatval)
                        NWSE = 0
                    else
                        NWSE = abs(matrix[i+1,j+1] - matrix[i-1,j-1]) / distDiag
                    end

                    # Record the steepest slope
                    slope[i,j] = nearest(T, max(NS,EW,NESW,NWSE))
                end
            end

            # Fill in edges too
            distEW = cellsize * km_per_lon
            distDiag = sqrt((distNS/2)^2 + distEW^2)

            # Left edge
            if (matrix[i+1,1] < minmatval) || (matrix[i-1,1] < minmatval)
                NS = 0
            else
                NS = abs(matrix[i+1,1] - matrix[i-1,1]) / distNS
            end
            if (matrix[i,2] < minmatval) || (matrix[i,1] < minmatval)
                EW = 0
            else
                EW = abs(matrix[i,2] - matrix[i,1]) / distEW
            end
            if (matrix[i+1,1] < minmatval) || (matrix[i-1,2] < minmatval)
                NESW = 0
            else
                NESW = abs(matrix[i+1,1] - matrix[i-1,2]) / distDiag
            end
            if (matrix[i+1,2] < minmatval) || (matrix[i-1,1] < minmatval)
                NWSE = 0
            else
                NWSE = abs(matrix[i+1,2] - matrix[i-1,1]) / distDiag
            end
            slope[i,1] = nearest(T, max(NS,EW,NESW,NWSE))

            # Right edge
            if (matrix[i+1,end] < minmatval) || (matrix[i-1,end] < minmatval)
                NS = 0
            else
                NS = abs(matrix[i+1,end] - matrix[i-1,end]) / distNS
            end
            if matrix[i,end]<minmatval || matrix[i,end-1]<minmatval
                EW = 0
            else
                EW = abs(matrix[i,end] - matrix[i,end-1]) / distEW
            end
            if (matrix[i+1,end-1] < minmatval) || (matrix[i-1,end] < minmatval)
                NEWS = 0
            else
                NESW = abs(matrix[i+1,end-1] - matrix[i-1,end]) / distDiag
            end
            if (matrix[i+1,end] < minmatval) || (matrix[i-1,end-1] < minmatval)
                NWSE = 0
            else
                NWSE = abs(matrix[i+1,end] - matrix[i-1,end-1]) / distDiag
            end
            slope[i,end] = nearest(T, max(NS,EW,NESW,NWSE))
        end

        # Fill in the top and bottom row
        distNS = cellsize * km_per_lat

        # Top row
        km_per_lon = cos(y_lat_cntr[1]*pi/180) * km_per_lat
        distEW = 2*cellsize*km_per_lon
        distDiag = sqrt(distNS^2+(distEW/2)^2)

        for j = 2:(size(matrix,2)-1)
            # Gradients, in meters per km
            if (matrix[2,j] < minmatval) || (matrix[1,j] < minmatval)
                NS = 0
            else
                NS = abs(matrix[2,j] - matrix[1,j]) / distNS
            end
            if (matrix[1,j+1] < minmatval) || (matrix[1,j-1] < minmatval)
                EW = 0
            else
                EW = abs(matrix[1,j+1] - matrix[1,j-1]) / distEW
            end
            if (matrix[2,j-1] < minmatval) || (matrix[1,j] < minmatval)
                NESW = 0
            else
                NESW = abs(matrix[2,j-1] - matrix[1,j]) / distDiag
            end
            if (matrix[2,j+1] < minmatval) || (matrix[1,j] < minmatval)
                NWSE = 0
            else
                NWSE = abs(matrix[2,j+1] - matrix[1,j]) / distDiag
            end
            slope[1,j] = nearest(T, max(NS,EW,NESW,NWSE))
        end
        slope[1,1] = 0
        slope[1,end] = 0

        # Bottom row
        km_per_lon = cos(y_lat_cntr[end] * pi/180) * km_per_lat
        distEW = 2 * cellsize * km_per_lon
        distDiag = sqrt(distNS^2 + (distEW/2)^2)
        for j = 2:(size(matrix,2)-1)
            # Gradients, in meters per Km
            if (matrix[end-1,j] < minmatval) || (matrix[end,j] < minmatval)
                NS = 0
            else
                NS = abs(matrix[end-1,j] - matrix[end,j]) / distNS
            end
            if (matrix[end,j+1] < minmatval) || (matrix[end,j-1] < minmatval)
                EW = 0
            else
                EW = abs(matrix[end,j+1] - matrix[end,j-1]) / distEW
            end
            if (matrix[end-1,j-1] < minmatval) || (matrix[end,j] < minmatval)
                NESW = 0
            else
                NESW = abs(matrix[end-1,j-1] - matrix[end,j]) / distDiag
            end
            if (matrix[end-1,j+1] < minmatval) || (matrix[end,j] < minmatval)
                NWSE = 0
            else
                NWSE = abs(matrix[end-1,j+1] - matrix[end,j]) / distDiag
            end
            slope[end,j] = nearest(T, max(NS,EW,NESW,NWSE))
        end
        slope[end,1] = 0
        slope[end,end] = 0

        return slope
    end
    export maxslope

    function aveslope(matrix, x_lon_cntr, y_lat_cntr, cellsize, T=UInt16; minmatval=-12000, maxmatval=9000, km_per_lat=111.1)
        # Returns slope in units/kilometer given a latitude-longitude grid of z-values

        # Allocate intermediate and output arrays
        distance = Array{Float64}(undef,8)
        local_slopes = Array{Float64}(undef,8)
        slope = Array{T}(undef,size(matrix))

        # Index offsets to cycle through:
        #         [N,NE,E,SE,S,SW,W,NW]
        ioffset = [-1,-1,0,1,1,1,0,-1]
        joffset = [0,1,1,1,0,-1,-1,-1]
        #
        # i.e. Layout:
        # 8 1 2
        # 7 x 3
        # 6 5 4

        # Distance between grid cell centers
        # N, S
        distance[[1,5]] .= cellsize * km_per_lat

        # Fill in the center first
        for i = 2:(size(matrix,1)-1)
            # Distance between grid cell centers
            km_per_lon = cos(y_lat_cntr[i]*pi/180) * km_per_lat
            distance[[3,7]] .= cellsize*km_per_lon; #E, W
            distance[[2,4,6,8]] .= sqrt(distance[1]^2+distance[3]^2)  # Diagonals

            # Center
            for j = 2:(size(matrix,2)-1)
                # Gradients, in matrix z-units per km
                here = matrix[i,j]
                if (here < minmatval) || (here > maxmatval)
                    slope[i,j] = 0
                else
                    for k = 1:8
                        there = matrix[i+ioffset[k], j+joffset[k]]
                        if (there < minmatval) || (there > maxmatval)
                            local_slopes[k] = 0
                        else
                            local_slopes[k] = abs(there-here) / distance[k]
                        end
                    end
                    # Record the average slope
                    slope[i,j] = nearest(T, nanmean(local_slopes))
                end
            end

            # Left edge
            here = matrix[i,1]
            if (here < minmatval) || (here > maxmatval)
                slope[i,1] = 0
            else
                for k = 1:5
                    there = matrix[i+ioffset[k], 1+joffset[k]]
                    if (there < minmatval) || (there > maxmatval)
                        local_slopes[k] = 0
                    else
                        local_slopes[k] = abs(there-here) / distance[k]
                    end
                end
                slope[i,1] = nearest(T, nanmean(local_slopes[1:5]))
            end

            # Right edge
            here = matrix[i,end]
            if (here < minmatval) || (here > maxmatval)
                slope[i,end] = 0
            else
                for k = [5,6,7,8,1]
                    there = matrix[i+ioffset[k], end+joffset[k]]
                    if (there < minmatval) || (there > maxmatval)
                        local_slopes[k] = 0
                    else
                        local_slopes[k] = abs(there-here) / distance[k]
                    end
                end
                slope[i,end] = nearest(T, nanmean(view(local_slopes, [5,6,7,8,1])))
            end
        end

        # Top row
        km_per_lon = cos(y_lat_cntr[1] * pi/180) * km_per_lat
        distance[[3,7]] .= cellsize * km_per_lon #E, W
        distance[[2,4,6,8]] .= sqrt(distance[1]^2 + distance[3]^2)  # Diagonals
        for j = 2:(size(matrix,2)-1)
            # Gradients, in matrix units per km
            here = matrix[1,j]
            if (here < minmatval) || (here > maxmatval)
                slope[1,j] = 0
            else
                for k=3:7
                    there = matrix[1+ioffset[k], j+joffset[k]]
                    if (there < minmatval) || (there > maxmatval)
                        local_slopes[k] = 0
                    else
                        local_slopes[k] = abs(there-here) / distance[k]
                    end
                end
                slope[1,j] = nearest(T, nanmean(view(local_slopes, 3:7)))
            end
        end
        slope[1,1] = 0
        slope[1,end] = 0

        # Bottom row
        km_per_lon = cos(y_lat_cntr[end] *pi/180) * km_per_lat
        distance[[3,7]] .= cellsize * km_per_lon #E, W
        distance[[2,4,6,8]] .= sqrt(distance[1]^2+distance[3]^2)  # Diagonals
        for j = 2:(size(matrix,2)-1)
            # Gradients, in matrix units per km
            here = matrix[end,j]
            if (here < minmatval) || (here > maxmatval)
                slope[end,j] = 0
            else
                for k = [7,8,1,2,3]
                    there = matrix[end+ioffset[k], j+joffset[k]]
                    if (there < minmatval) || (there > maxmatval)
                        local_slopes[k] = 0
                    else
                        local_slopes[k] = abs(there-here) / distance[k]
                    end
                end
                slope[end,j] = nearest(T, nanmean(view(local_slopes, [7,8,1,2,3])))
            end
        end
        slope[end,1] = 0
        slope[end,end] = 0

        return slope
    end
    export aveslope

## --- Generate random latitude and longitude pairs uniformly distributed across the globe

    CONST_180_PI = 180/pi

    function randlatlon(n::Integer; land=false)
        if land
            c = 0
            while c < n
                lats, lons = 90 .- CONST_180_PI*acos.(2*rand(5n) .- 1), rand(5n)*360 .- 180
                notland = .!find_land(lats, lons)
                @inbounds for i in eachindex(notland)
                    c += !notland[i]
                    if c > n
                        notland[i] = true
                    end
                end
                c < n && (c=0)
                deleteat!(lats, notland)
                deleteat!(lons, notland)
            end
        else
            lats, lons = 90 .- CONST_180_PI*acos.(2*rand(n) .- 1), rand(n)*360 .- 180
        end
        return lats, lons
    end
    function randlatlon()
        90 - CONST_180_PI*acos(2*rand() - 1), rand()*360 - 180
    end
    export randlatlon

## --- Calculate distance uncertainty in arc degrees
    """
    ```julia
    haversine(lat₁, lon₁, lat₂, lon₂)
    ```
    Calculate the arc degree distance between two decimal degree points (lat₁, lon₁) and
    (lat₂, lon₂).

    """
    function haversine(lat₁, lon₁, lat₂, lon₂)
        lat₁ᵣ, lon₁ᵣ, lat₂ᵣ, lon₂ᵣ = (lat₁, lon₁, lat₂, lon₂) .* (pi/180)
        dist = acos(sin(lat₁ᵣ) * sin(lat₂ᵣ) + cos(lat₁ᵣ) * cos(lat₂ᵣ) * cos(lon₁ᵣ - lon₂ᵣ))
        return dist * 180/pi
    end
    export haversine

## --- Calculate maximum arc-degree distance between a series of points
    """
    ```julia
    dist_uncert(lats, lons)
    ```

    Find the decimal degree center and associated uncertainty (in arc degrees) from
    lists `lats` and `lons` of decimal degree coordinates.

    ### Examples
    ```julia
    (lat_ctr, lon_ctr, uncertainty) = dist_uncert(lats, lons)
    ```
    """
    function dist_uncert(lats, lons)
        @assert eachindex(lats) == eachindex(lons)
        latc, lonc = centroid(lats, lons)
        maxdist = zero(float(eltype(lats)))
        for i in eachindex(lats)
            for j in 1+firstindex(lats):lastindex(lats)
                # If a point is compared to itself, distance is 0; comparison is susceptible to roundoff error
                if i != j
                    dist = haversine(lats[i], lons[i], lats[j], lons[j])
                    dist > maxdist && (maxdist = dist)
                end
            end
        end
        return latc, lonc, maxdist/2
    end
    export dist_uncert

## --- Other lat and lon conversions

    """
    ```julia
    centroid(lats, lons)
    ```
    Return the centroid of a set of latitudes and longitudes on a sphere.
    """
    function centroid(lats::AbstractArray{T1}, lons::AbstractArray{T2}) where {T1,T2}
        T = float(promote_type(T1, T2))
        x, y, z = similar(lats, T), similar(lats, T), similar(lats, T)
        @inbounds for i in eachindex(lats, lons)
            φ = deg2rad(90 - lats[i])
            θ = deg2rad(lons[i])
            x[i], y[i], z[i] = cartesian(one(T), φ, θ)
        end
        x₀ = nanmean(x)
        y₀ = nanmean(y)
        z₀ = nanmean(z)
        ρ, φ, θ = spherical(x₀, y₀, z₀)
        latc = 90 - rad2deg(φ)
        lonc = rad2deg(θ)
        return latc, lonc
    end
    export centroid

    """
    ```julia
    x, y, z = cartesian(ρ, φ, θ)
    ```
    Convert from spherical coordinates (`ρ`,`φ`,`θ`) to cartesian coordinates (`x`,`y`,`z`).
    """
    function cartesian(ρ::Number, φ::Number, θ::Number)
        x = ρ * sin(φ) * cos(θ)
        y = ρ * sin(φ) * sin(θ)
        z = ρ * cos(φ)
        return x, y, z
    end
    export cartesian

    """
    ```julia
    ρ, θ, φ = cartesian(x, y, z)
    ```
    Convert from cartesian coordinates (`x`,`y`,`z`) to spherical coordinates (`ρ`,`φ`,`θ`).
    """
    function spherical(x::Number, y::Number, z::Number)
        ρ = sqrt(x^2 + y^2 + z^2)
        φ = acos(z/ρ)
        θ = atan(y,x)
        return ρ, φ, θ
    end
    export spherical


## --- End of File
