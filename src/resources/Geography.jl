## --- Geolcont

    continentcolors = parse.(Color, ["#333399","#0066CC","#06A9C1","#66CC66","#FFCC33","#FFFF00","#FFFFFF"])
    export continentcolors

    continents = ["Africa","Eurasia","North America","South America","Australia","Antarctica","NA"]
    export continents

    """
    ```julia
    find_geolcont(lat,lon)
    ```
    Find which geographic continent a sample originates from.

    Continents:
      1: "Africa"
      2: "Eurasia"
      3: "North America"
      4: "South America"
      5: "Australia"
      6: "Antarctica"
      7: "NA"

    See also: `continents`, `continentcolors`.

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
        @assert eachindex(lat) == eachindex(lon)

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
        for i ∈ eachindex(lat)
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


## --- geolprov

    """
    ```julia
    find_geolprov(lat,lon)
    ```
    Find which tectonic setting a sample originates from.

    Settings:
      10: Accreted Arc
      11: Island Arc
      12: Continental Arc
      13: Collisional orogen
      20: Extensional
      21: Rift
      22: Plume
      31: Shield
      32: Platform
      33: Basin
      00: No data

    ## Examples
    ```julia
    julia> find_geolprov(43.702245, -72.0929)
    0-dimensional Array{Int64, 0}:
    10

    julia> lat = rand(4)*180 .- 90
    4-element Vector{Float64}:
     -28.352224011759773
      14.521710123066882
      43.301961981794335
      79.26368353708557

    julia> lon = rand(4)*360 .- 180
    4-element Vector{Float64}:
       5.024149409750521
     161.04362679392233
     123.21726489255786
     -54.34797401313695

    julia> find_geolprov(lat, lon)
    4-element Vector{Int64}:
      0
      0
     32
      0
    ```
    """
    function find_geolprov(lat, lon)
        @assert eachindex(lat) == eachindex(lon)
        filepath = joinpath(moduleresourcepath,"geolprov.h5")
        geolprov = h5read(filepath, "geolprov")

        result = zeros(Int, size(lat))
        for i ∈ eachindex(lat)
            if -180 < lon[i] <= 180 && -90 <= lat[i] < 90
                x = ceil(Int, (lon[i]+180) * 2161/360)
                y = ceil(Int, (90-lat[i]) * 1801/180)
                result[i] = geolprov[y,x]
            end
        end
        return result
    end
    export find_geolprov
