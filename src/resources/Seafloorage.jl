# --- Müller et al. seafloor age and spreading rate

    """
    ```julia
    get_seafloorage(varname="")
    ```
    Read seafloor age file from HDF5 storage, downloading from cloud if necessary.

    Available `varname`s (variable names) include:
      "seafloorage"
      "seafloorage_sigma",
      "seafloorrate"
      "information"
      "reference"

    Units are millions of years for age and mm/yr for rate.

    Reference:
    M⁠⁠üller, R. D., M. Sdrolias, C. Gaina, and W. R. Roest (2008). Age, spreading
    rates, and spreading asymmetry of the world's ocean crust, Geochem. Geophys. Geosyst., 9,
    Q04006, doi:10.1029/2007GC001743. ftp://ftp.es.usyd.edu.au/pub/agegrid/2008/Grids/

    See also: `find_seafloorage`.

    ## Examples
    ```julia
    julia> get_seafloorage()
    Dict{String, Any} with 5 entries:
      "seafloorrate"      => [6.83 6.79 … 6.91 6.87; 6.82 6.78 … 6.91 6.…
      "reference"         => "Muller, R. D., M. Sdrolias, C. Gaina, and …
      "seafloorage"       => [8.02 8.13 … 7.8 7.91; 8.01 8.12 … 7.81 7.9…
      "information"       => "Mercator projection, from 80.738 to -80.73…
      "seafloorage_sigma" => [2.15 2.19 … 2.08 2.11; 2.15 2.19 … 2.08 2.…

    julia> get_seafloorage("seafloorage")
    8640×10800 Matrix{Float64}:
       8.02    8.13    8.24    8.36  …    7.6     7.7     7.8     7.91
       8.01    8.12    8.23    8.34       7.61    7.71    7.81    7.91
       8.01    8.11    8.22    8.32       7.62    7.72    7.81    7.91
       8.01    8.11    8.2     8.3        7.64    7.73    7.82    7.91
       ⋮                             ⋱
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN
    ```
    """
    function get_seafloorage(varname="")
        # Available variable names: "seafloorage", "seafloorage_sigma",
        # "seafloorrate", "information", and "reference".

        # Construct file path
        filedir = joinpath(resourcepath,"seafloorage")
        filepath = joinpath(filedir,"seafloorage.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            @info "Downloading seafloorage.h5 from google cloud storage to $filedir"
            run(`mkdir -p $filedir`)
            Downloads.download("https://storage.googleapis.com/statgeochem/seafloorage.references.txt", joinpath(filedir,"seafloorage.references.txt"))
            Downloads.download("https://storage.googleapis.com/statgeochem/seafloorage.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/"*varname)
    end
    export get_seafloorage


    """
    ```julia
    find_seafloorage([sfdata], lat, lon)
    ```
    Find the age of the seafloor at positions (`lat`, `lon`) on the ocean crust,
    using the Müller et al. (2008) dataset.

    Units are millions of years for age and mm/yr for rate.

    Reference:
    M⁠⁠üller, R. D., M. Sdrolias, C. Gaina, and W. R. Roest (2008). Age, spreading
    rates, and spreading asymmetry of the world's ocean crust, Geochem. Geophys. Geosyst., 9,
    Q04006, doi:10.1029/2007GC001743. ftp://ftp.es.usyd.edu.au/pub/agegrid/2008/Grids/

    See also: `get_seafloorage`.

    ## Examples
    ```julia
    julia> sfdata = get_seafloorage("seafloorage")
    8640×10800 Matrix{Float64}:
       8.02    8.13    8.24    8.36  …    7.6     7.7     7.8     7.91
       8.01    8.12    8.23    8.34       7.61    7.71    7.81    7.91
       8.01    8.11    8.22    8.32       7.62    7.72    7.81    7.91
       8.01    8.11    8.2     8.3        7.64    7.73    7.82    7.91
       ⋮                             ⋱
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN
     NaN     NaN     NaN     NaN        NaN     NaN     NaN     NaN

    julia> find_seafloorage(sfdata, 43.702245, -40)
    0-dimensional Array{Float64, 0}:
    75.95
    ```
    """
    find_seafloorage(lat, lon) = find_seafloorage(get_seafloorage("seafloorage"), lat, lon)
    find_seafloorage(sfdata::Dict, lat, lon) = find_seafloorage(sfdata["seafloorage"], lat, lon)
    function find_seafloorage(sfdata::AbstractArray, lat, lon)

        # Interpret user input
        length(lat) != length(lon) && error("`lat` and `lon` must be of equal length\n")

        # Find the column numbers (using mod to convert lon from -180:180 to 0:360
        x = floor.(Int, mod.(lon, 360) * 10800/360) .+ 1

        # find the y rows, converting from lat to Mercator (lat -80.738:80.738)
        y = 4320 .- floor.(Int, 8640 * asinh.(tan.(lat*pi/180)) / asinh.(tan.(80.738*pi/180)) / 2 ) .+ 1

        # Make and fill output array
        out=Array{Float64}(undef,size(x))
        for i=1:length(x)
            # If there is out data for row(i), col(i)
            if isnan(x[i]) || isnan(y[i]) || x[i]<1 || x[i]>10800 || y[i]<1 || y[i]>8640
                out[i] = NaN
            else
                # Then fill in the output data (Age, Age_Min, Age_Max)
                out[i] = sfdata[y[i], x[i]]
            end
        end

        return out
    end
    export find_seafloorage

## --- End of File
