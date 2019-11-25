# --- Müller et al. seafloor age and spreading rate

    # Read seafloorage file from HDF5 storage, downloading from cloud if necessary
    function get_seafloorage(varname = "")
        # Available variable names: "seafloorage", "seafloorage_sigma",
        # "seafloorrate", "information", and "reference".

        # Construct file path
        filedir = joinpath(resourcepath,"seafloorage")
        filepath = joinpath(filedir,"seafloorage.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading seafloorage.h5 from google cloud storage to $filedir\n")
            run(`mkdir -p $filedir`)
            download("https://storage.googleapis.com/statgeochem/seafloorage.references.txt", joinpath(filedir,"seafloorage.references.txt"))
            download("https://storage.googleapis.com/statgeochem/seafloorage.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/"*varname)
    end
    export get_seafloorage

    # Parse seafloorage, seafloorage_sigma, or seafloorrate from file
    # data = find_seafloorage(sfdata, lat, lon)
    function find_seafloorage(sfdata,lat,lon)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(sfdata,Dict)
            data = sfdata["seafloorage"]
        elseif isa(srtm15plus, Array)
            data = sfdata
        else
            error("wrong sfdata variable")
        end

        # Find the column numbers (using mod to convert lon from -180:180 to 0:360
        x = floor.(Int, mod.(lon, 360) * 10800/360) + 1

        # find the y rows, converting from lat to Mercator (lat -80.738:80.738)
        y = 4320 - floor.(Int, 8640 * asinh.(tan.(lat*pi/180)) / asinh.(tan.(80.738*pi/180)) / 2 ) + 1

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
