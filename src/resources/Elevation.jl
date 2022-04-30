## --- Geolcont

    continentcolors = parse.(Color, ["#333399","#0066CC","#06A9C1","#66CC66","#FFCC33","#FFFF00","#FFFFFF"])
    export continentcolors

    continents = ["Africa","Eurasia","North America","South America","Australia","Antarctica","NA"]
    export continents

    """
    ```julia
    find_geolcont(lat,lon)
    ```
    Find which continent a sample originates from.

    Continents:
      1: "Africa"
      2: "Eurasia"
      3: "North America"
      4: "South America"
      5: "Australia"
      6: "Antarctica"
      7: "NA"

    ## Examples
    ```julia
    julia> find_geolcont(43.702245, -72.0929)
    0-dimensional Array{Int64, 0}:
    3

    julia> continents[find_geolcont(43.702245, -72.0929)]
    0-dimensional Array{String, 0}:
    "North America"
    ```
    """
    function find_geolcont(lat,lon)
        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Construct file path
        filedir = joinpath(resourcepath,"maps")
        filepath = joinpath(filedir,"geolcontwshelf.png")

        # Download png map from Google Cloud if necessary
        if ~isfile(filepath)
            @info "Downloading map to $filedir"
            run(`mkdir -p $filedir`)
            Downloads.download("https://storage.googleapis.com/statgeochem/geolcontwshelf.png", filepath)
        end

        img = load(filepath)

        ind = fill(7,size(img))
        for i=1:6
         ind[img .== continentcolors[i]] .= i
        end

        # Create and fill output vector
        contindex = Array{Int}(undef,size(lat))
        for i=1:length(lat)
         if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
             # Result is unknown if either input is NaN or out of bounds
             contindex[i] = 7
         else
             # Convert latitude and longitude into indicies of the elevation map array
             # Note that STRTM15 plus has N+1 columns where N = 360*sf
             row = 1 + trunc(Int,(90-lat[i])*512/180)
             col = 1 + trunc(Int,(180+lon[i])*512/180)
             # Find result by indexing
             contindex[i] = ind[row,col]
         end
        end

        return contindex
    end
    export find_geolcont

## --- ETOPO1 (1 arc minute topography)

    """
    ```julia
    get_etopo([varname])
    ```
    Read ETOPO1 (1 arc minute topography) file from HDF5 storage, downloading
    from cloud if necessary.

    Available `varname`s (ETOPO variable names) include:
      "elevation"
      "y_lat_cntr"
      "x_lon_cntr"
      "cellsize"
      "scalefactor"
      "reference"

    Units are meters of elevation and decimal degrees of latitude and longitude.

    ## Examples
    ```julia
    julia> get_etopo()
    Dict{String, Any} with 6 entries:
      "cellsize"    => 0.0166667
      "scalefactor" => 60
      "x_lon_cntr"  => [-179.992, -179.975, -179.958, -179.942, -179.925, -1…
      "reference"   => "Amante, C. and B.W. Eakins, 2009. ETOPO1 1 Arc-Minut…
      "y_lat_cntr"  => [-89.9917, -89.975, -89.9583, -89.9417, -89.925, -89.…
      "elevation"   => [-58.0 -58.0 … -58.0 -58.0; -61.0 -61.0 … -61.0 -61.0…

    julia> get_etopo("elevation")
    10800×21600 Matrix{Float64}:
      -58.0    -58.0    -58.0  …    -58.0    -58.0    -58.0
      -61.0    -61.0    -61.0       -61.0    -61.0    -61.0
      -62.0    -63.0    -63.0       -63.0    -63.0    -62.0
      -61.0    -62.0    -62.0       -62.0    -62.0    -61.0
        ⋮                      ⋱
     -4226.0  -4226.0  -4227.0     -4227.0  -4227.0  -4227.0
     -4228.0  -4228.0  -4229.0     -4229.0  -4229.0  -4229.0
     -4229.0  -4229.0  -4229.0     -4229.0  -4229.0  -4229.0
    ```
    """
    function get_etopo(varname="")
        # Available variable names: "elevation", "y_lat_cntr", "x_lon_cntr",
        # "cellsize", "scalefactor", and "reference". Units are meters of
        # elevation and decimal degrees of latitude and longitude

        # Construct file path
        filedir = joinpath(resourcepath,"etopo")
        filepath = joinpath(filedir,"etopo1.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            @info "Downloading etopo1.h5 from google cloud storage to $filedir"
            run(`mkdir -p $filedir`)
            Downloads.download("https://storage.googleapis.com/statgeochem/etopo1.references.txt", joinpath(filedir,"etopo1.references.txt"))
            Downloads.download("https://storage.googleapis.com/statgeochem/etopo1.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath, "vars/"*varname)
    end
    export get_etopo


    """
    ```julia
    find_etopoelev([etopo], lat, lon, [T=Float64])
    ```
    Find the elevation of points at position (`lat`, `lon`) on the surface of the
    Earth, using the ETOPO1 one-arc-degree elevation model.

    Units are meters of elevation and decimal degrees of latitude and longitude.

    ## Examples
    ```julia
    julia> etopo=get_etopo("elevation")
    10800×21600 Matrix{Float64}:
       -58.0    -58.0    -58.0  …    -58.0    -58.0    -58.0
       -61.0    -61.0    -61.0       -61.0    -61.0    -61.0
       -62.0    -63.0    -63.0       -63.0    -63.0    -62.0
       -61.0    -62.0    -62.0       -62.0    -62.0    -61.0
         ⋮                      ⋱
     -4226.0  -4226.0  -4227.0     -4227.0  -4227.0  -4227.0
     -4228.0  -4228.0  -4229.0     -4229.0  -4229.0  -4229.0
     -4229.0  -4229.0  -4229.0     -4229.0  -4229.0  -4229.0

    julia> find_etopoelev(etopo, 43.702245, -72.0929)
    0-dimensional Array{Float64, 0}:
    294.0
    ```
    """
    find_etopoelev(lat,lon) = find_etopoelev(get_etopo(),lat,lon)
    find_etopoelev(etopo::Dict, lat, lon) = find_etopoelev(etopo["elevation"], lat, lon)
    function find_etopoelev(etopo::AbstractArray, lat, lon, T=Float64)
        # Interpret user input
        length(lat) != length(lon) && error("lat and lon must be of equal length\n")

        # Scale factor (cells per degree) = 60 = arc minutes in an arc degree
        sf = 60
        maxrow = 180 * sf
        maxcol = 360 * sf

        # Create and fill output vector
        result = Array{T}(undef,size(lat))
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                result[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                row = 1 + trunc(Int,(90+lat[i])*sf)
                if row == (maxrow+1) # Edge case
                    row = maxrow
                end

                col = 1 + trunc(Int,(180+lon[i])*sf)
                if col == (maxcol+1) # Edge case
                    col = maxcol
                end

                # Find result by indexing
                result[i] = etopo[row,col]
            end
        end

        return result
    end
    export find_etopoelev


## --- SRTM15_PLUS (15 arc second topography)

    # Read srtm15plus file from HDF5 storage, downloading from cloud if necessary
    function get_srtm15plus(varname="")
        # Available variable names: "elevation", "y_lat_cntr", "x_lon_cntr",
        # "nanval", "cellsize", "scalefactor", and "reference". Units are
        # meters of elevation and decimal degrees of latitude and longitude

        # Construct file path
        filedir = joinpath(resourcepath,"srtm15plus")
        filepath = joinpath(filedir,"srtm15plus.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            @info "Downloading srtm15plus.h5 from google cloud storage to $filedir"
            run(`mkdir -p $filedir`)
            Downloads.download("https://storage.googleapis.com/statgeochem/srtm15plus.references.txt", joinpath(filedir,"srtm15plus.references.txt"))
            Downloads.download("https://storage.googleapis.com/statgeochem/srtm15plus.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/"*varname)
    end
    export get_srtm15plus

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the SRTM15plus 15-arc-second elevation model.
    find_srtm15plus(lat,lon) = find_srtm15plus(get_srtm15plus(),lat,lon)
    find_srtm15plus(srtm::Dict, lat, lon) = find_srtm15plus(srtm["elevation"], lat, lon)
    function find_srtm15plus(srtm::AbstractArray, lat, lon, T=Float64)
        # Interpret user input
        length(lat) != length(lon) && error("lat and lon must be of equal length")

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out = Array{T}(undef,size(lat))
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90+lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)
                # Find result by indexing
                res = srtm[row,col]
                if res == -32768
                    out[i] = NaN
                else
                    out[i] = res
                end
            end
        end
        return out
    end
    export find_srtm15plus

## --- End of File
