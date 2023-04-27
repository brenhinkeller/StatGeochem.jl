## --- CRUST 1.0

    function crust1layer(layer::Number)
        if isinteger(layer) && (1 <= layer <= 9)
            Int(layer)
        else
            @error """crust1 layer $layer not found.
            Available layers include:
                1 | :water
                2 | :ice
                3 | :upper_sediments
                4 | :middle_sediments
                5 | :lower_sediments
                6 | :upper_crust
                7 | :middle_crust
                8 | :lower_crust
                9 | :mantle
            """
        end
    end
    function crust1layer(layer::Symbol)
        if layer===:water
            1
        elseif layer===:ice
            2
        elseif layer===:upper_sediments
            3
        elseif layer===:middle_sediments
            4
        elseif layer===:lower_sediments
            5
        elseif layer===:upper_crust
            6
        elseif layer===:middle_crust
            7
        elseif layer===:lower_crust
            8
        elseif layer===:mantle
            9
        else
            @error """crust1 layer $layer not found.
            Available layers include:
                1 | :water
                2 | :ice
                3 | :upper_sediments
                4 | :middle_sediments
                5 | :lower_sediments
                6 | :upper_crust
                7 | :middle_crust
                8 | :lower_crust
                9 | :mantle
            """
        end
    end

    """
    ```julia
    find_crust1_layer(lat,lon,layer)
    ```
    Return all point data (Vp, Vs, Rho, layer thickness) for a given `lat`itude,
    `lon`gitude, and crustal `layer`.

    Accepts `lat` and `lon` both as `Numbers` and as `AbstractArray`s, but given
    the overhead of opening and reading the crust1 files, you should generally
    aim to provide large arrays with as many values in a single query as possible.

    Available `layer`s:
    ```
        1 | :water
        2 | :ice
        3 | :upper_sediments
        4 | :middle_sediments
        5 | :lower_sediments
        6 | :upper_crust
        7 | :middle_crust
        8 | :lower_crust
        9 | :mantle
    ```
    Results are returned in form `(Vp, Vs, Rho, thickness)`

    ## Examples
    ```julia
    julia> vp, vs, rho, thickness = find_crust1_layer([43.702245], [-72.0929], 8)
    ([7.0], [3.99], [2950.0], [7.699999999999999])
    ```
    """
    function find_crust1_layer(lat,lon,layer)
        # Get Vp, Vs, Rho, and thickness for a given lat, lon, and crustal layer.
        @assert eachindex(lat) == eachindex(lon)
        layerindex = crust1layer(layer)

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        vp = Array{Float64,3}(undef,nlayers,nlat,nlon)
        vs = Array{Float64,3}(undef,nlayers,nlat,nlon)
        rho = Array{Float64,3}(undef,nlayers,nlat,nlon)
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        vpfile = open(artifact"crust1/crust1.vp", "r")
        vsfile = open(artifact"crust1/crust1.vs", "r")
        rhofile = open(artifact"crust1/crust1.rho", "r")
        bndfile = open(artifact"crust1/crust1.bnds", "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              vp[:,j,i] = delim_string_parse(readline(vpfile), ' ', Float64, merge=true)
              vs[:,j,i] = delim_string_parse(readline(vsfile), ' ', Float64, merge=true)
              rho[:,j,i] = delim_string_parse(readline(rhofile), ' ', Float64, merge=true) * 1000 # convert to kg/m3
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(vpfile)
        close(vsfile)
        close(rhofile)
        close(bndfile)

        # Allocate output arrays
        vpout = Array{Float64}(undef,size(lat))
        vsout = Array{Float64}(undef,size(lat))
        rhoout = Array{Float64}(undef,size(lat))
        thkout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        @inbounds for j ∈ eachindex(lat)
            # Avoid edge cases at lat = -90.0, lon = 180.0
            lonⱼ = mod(lon[j] + 180, 360) - 180
            latⱼ = lat[j]

            if -90 < latⱼ < 90 && -180 < lonⱼ < 180
                # Convert lat and lon to index
                ilat = 91 - ceil(Int,latⱼ)
                ilon = 181 + floor(Int,lonⱼ)

                vpout[j] = vp[layerindex,ilat,ilon]
                vsout[j] = vs[layerindex,ilat,ilon]
                rhoout[j] = rho[layerindex,ilat,ilon]
                thkout[j] = bnd[layerindex,ilat,ilon] - bnd[layerindex+1,ilat,ilon]
            else
                vpout[j] = NaN
                vsout[j] = NaN
                rhoout[j] = NaN
                thkout[j] = NaN
            end
        end

        # The end
        return (vpout, vsout, rhoout, thkout)
    end
    export find_crust1_layer

    """
    ```julia
    find_crust1_seismic(lat,lon,layer)
    ```
    Return all seismic data (Vp, Vs, Rho) for a given `lat`itude, `lon`gitude,
    and crustal `layer`.

    Accepts `lat` and `lon` both as `Numbers` and as `AbstractArray`s, but given
    the overhead of opening and reading the crust1 files, you should generally
    aim to provide large arrays with as many values in a single query as possible.

    Available `layer`s:
    ```
        1 | :water
        2 | :ice
        3 | :upper_sediments
        4 | :middle_sediments
        5 | :lower_sediments
        6 | :upper_crust
        7 | :middle_crust
        8 | :lower_crust
        9 | :mantle
    ```
    Results are returned in form `(Vp, Vs, Rho, thickness)`

    ## Examples
    ```julia
    julia> vp, vs, rho = find_crust1_seismic([43.702245], [-72.0929], 8)
    ([7.0], [3.99], [2950.0])

    julia> vp, vs, rho = find_crust1_seismic([43.702245], [-72.0929], :lower_crust)
    ([7.0], [3.99], [2950.0])
    ```
    """
    function find_crust1_seismic(lat,lon,layer)
        # Vp, Vs, and Rho for a given lat, lon, and crustal layer.
        @assert eachindex(lat) == eachindex(lon)
        layerindex = crust1layer(layer)

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        vp = Array{Float64,3}(undef,nlayers,nlat,nlon)
        vs = Array{Float64,3}(undef,nlayers,nlat,nlon)
        rho = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        vpfile = open(artifact"crust1/crust1.vp", "r")
        vsfile = open(artifact"crust1/crust1.vs", "r")
        rhofile = open(artifact"crust1/crust1.rho", "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              vp[:,j,i] = delim_string_parse(readline(vpfile), ' ', Float64, merge=true)
              vs[:,j,i] = delim_string_parse(readline(vsfile), ' ', Float64, merge=true)
              rho[:,j,i] = delim_string_parse(readline(rhofile), ' ', Float64, merge=true) * 1000 # convert to kg/m3
          end
        end

        # Close data files
        close(vpfile)
        close(vsfile)
        close(rhofile)

        # Allocate output arrays
        vpout = Array{Float64}(undef,size(lat))
        vsout = Array{Float64}(undef,size(lat))
        rhoout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        @inbounds for j ∈ eachindex(lat)
            # Avoid edge cases at lat = -90.0, lon = 180.0
            lonⱼ = mod(lon[j] + 180, 360) - 180
            latⱼ = lat[j]

            if -90 < latⱼ < 90 && -180 < lonⱼ < 180
                # Convert lat and lon to index
                ilat = 91 - ceil(Int,latⱼ)
                ilon = 181 + floor(Int,lonⱼ)

                vpout[j] = vp[layerindex,ilat,ilon]
                vsout[j] = vs[layerindex,ilat,ilon]
                rhoout[j] = rho[layerindex,ilat,ilon]
            else
                vpout[j] = NaN
                vsout[j] = NaN
                rhoout[j] = NaN
            end
        end

        # The end
        return (vpout, vsout, rhoout)
    end
    export find_crust1_seismic

    """
    ```julia
    find_crust1_thickness(lat,lon,layer)
    ```
    Return layer thickness for a crust 1.0 `layer` at a given `lat`itude and
    `lon`gitude.

    Accepts `lat` and `lon` both as `Numbers` and as `AbstractArray`s, but given
    the overhead of opening and reading the crust1 files, you should generally
    aim to provide large arrays with as many values in a single query as possible.

    Available `layer`s:
    ```
        1 | :water
        2 | :ice
        3 | :upper_sediments
        4 | :middle_sediments
        5 | :lower_sediments
        6 | :upper_crust
        7 | :middle_crust
        8 | :lower_crust
        9 | :mantle
    ```
    Results are returned in form `(Vp, Vs, Rho, thickness)`

    ## Examples
    ```julia
    julia> find_crust1_thickness([43.702245], [-72.0929], 8)
    1-element Vector{Float64}:
     7.699999999999999

    julia> find_crust1_thickness([43.702245], [-72.0929], :lower_crust)
    1-element Vector{Float64}:
    7.699999999999999
    ```
    """
    function find_crust1_thickness(lat,lon,layer)
        # Layer thickness for a given lat, lon, and crustal layer.
        @assert eachindex(lat) == eachindex(lon)
        layerindex = crust1layer(layer)

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        bndfile = open(artifact"crust1/crust1.bnds", "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(bndfile)

        # Allocate output arrays
        thkout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        @inbounds for j ∈ eachindex(lat)
            # Avoid edge cases at lat = -90.0, lon = 180.0
            lonⱼ = mod(lon[j] + 180, 360) - 180
            latⱼ = lat[j]

            if -90 < latⱼ < 90 && -180 < lonⱼ < 180
                # Convert lat and lon to index
                ilat = 91 - ceil(Int,latⱼ)
                ilon = 181 + floor(Int,lonⱼ)

                thkout[j] = bnd[layerindex,ilat,ilon]-bnd[layerindex+1,ilat,ilon]
            else
                thkout[j] = NaN
            end
        end

        # The end
        return thkout
    end
    export find_crust1_thickness

    """
    ```julia
    find_crust1_base(lat,lon,layer)
    ```
    Return elevation (relative to sea level) of the layer base for a crust 1.0
    `layer` at a given `lat`itude and `lon`gitude.

    Accepts `lat` and `lon` both as `Numbers` and as `AbstractArray`s, but given
    the overhead of opening and reading the crust1 files, you should generally
    aim to provide large arrays with as many values in a single query as possible.

    Available `layer`s:
    ```
        1 | :water
        2 | :ice
        3 | :upper_sediments
        4 | :middle_sediments
        5 | :lower_sediments
        6 | :upper_crust
        7 | :middle_crust
        8 | :lower_crust
        9 | :mantle
    ```
    Results are returned in form `(Vp, Vs, Rho, thickness)`

    ## Examples
    ```julia
    julia> find_crust1_base([43.702245], [-72.0929], 8)
    1-element Vector{Float64}:
     -36.26

    julia> find_crust1_base([43.702245], [-72.0929], :lower_crust)
    1-element Vector{Float64}:
    -36.26
    ```
    """
    function find_crust1_base(lat,lon,layer)
        # Depth to layer base for a given lat, lon, and crustal layer.
        @assert eachindex(lat) == eachindex(lon)
        layerindex = crust1layer(layer)

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        bndfile = open(artifact"crust1/crust1.bnds", "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(bndfile)

        # Allocate output arrays
        baseout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        @inbounds for j ∈ eachindex(lat)
            # Avoid edge cases at lat = -90.0, lon = 180.0
            lonⱼ = mod(lon[j] + 180, 360) - 180
            latⱼ = lat[j]

            if -90 < latⱼ < 90 && -180 < lonⱼ < 180
                # Convert lat and lon to index
                ilat = 91 - ceil(Int,latⱼ)
                ilon = 181 + floor(Int,lonⱼ)

                baseout[j] = bnd[layerindex+1,ilat,ilon]
            else
                baseout[j] = NaN
            end
        end

        # The end
        return baseout
    end
    export find_crust1_base

## --- End of File
