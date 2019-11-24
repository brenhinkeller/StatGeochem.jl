## --- Geolcont

    continentcolors = [RGB4{N0f8}(0.2,0.2,0.6)
              RGB4{N0f8}(0.0,0.4,0.8)
              RGB4{N0f8}(0.024,0.663,0.757)
              RGB4{N0f8}(0.4,0.8,0.4)
              RGB4{N0f8}(1.0,0.8,0.2)
              RGB4{N0f8}(1.0,1.0,0.0)
              RGB4{N0f8}(1.0,1.0,1.0)
    ]
    export continentcolors

    continents = ["Africa","Eurasia","North America","South America","Australia","Antarctica","NA"]
    export continents

    # Find which continent a sample originates from
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
            print("Downloading map to $filedir\n")
            run(`mkdir -p $filedir`)
            download("https://storage.googleapis.com/statgeochem/geolcontwshelf.png", filepath)
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

## --- TC1

    tc1_550 = readdlm(joinpath(moduleresourcepath,"tc1_550.csv"), ',')
    tc1_1300 = readdlm(joinpath(moduleresourcepath,"tc1_1300.csv"), ',')
    tc1_age = readdlm(joinpath(moduleresourcepath,"tc1_age.csv"), ',', Int)

    function find_tc1_crust(lat::Number,lon::Number)
        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            crust=tc1_550[i,j]
        else
            crust = NaN
        end
        return crust
    end
    function find_tc1_crust(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Query the tc1_550 array for our lat and lon
        crust = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                crust[n]=tc1_550[i,j]
            end
        end
        return crust
    end
    export find_tc1_crust

    function find_tc1_lith(lat::Number,lon::Number)
        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            lith=tc1_1300[i,j]
        else
            lith = NaN
        end
        return lith
    end
    function find_tc1_lith(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Query the tc1_1300 array for our lat and lon
        lith = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                lith[n]=tc1_1300[i,j]
            end
        end
        return lith
    end
    export find_tc1_lith

    function find_tc1_age(lat::Number,lon::Number)
        ages=[ NaN  NaN  NaN
                25    0   50
               150   50  250
               395  250  540
               695  540  850
               975  850 1100
              1400 1100 1700
              2100 1700 2500
              2750 2500 3000
              3250 3000 3500]

        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            result = (ages[tc1_age[i,j],:]...,)
        else
            result = (NaN, NaN, NaN,)
        end
        return result
    end
    function find_tc1_age(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        ages=[ NaN  NaN  NaN
                25    0   50
               150   50  250
               395  250  540
               695  540  850
               975  850 1100
              1400 1100 1700
              2100 1700 2500
              2750 2500 3000
              3250 3000 3500]

        # Query the tc1_1300 array for our lat and lon
        age = fill(NaN, size(lat))
        minage = fill(NaN, size(lat))
        maxage = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                ageindex = tc1_age[i,j]
                age[n] = ages[ageindex,1]
                minage[n] = ages[ageindex,2]
                maxage[n] = ages[ageindex,3]
            end
        end
        return (age, minage, maxage)
    end
    export find_tc1_age

## --- CRUST 1.0

    # Download CRUST 1.0 data and references from cloud
    function get_crust1()
        # Available variable names: "seafloorage", "seafloorage_sigma",
        # "seafloorrate", "information", and "reference".

        # Construct file paths
        filedir = joinpath(resourcepath,"crust1")
        referencepath = joinpath(filedir,"crust1.references.txt")
        vppath = joinpath(filedir,"crust1.vp")
        vspath = joinpath(filedir,"crust1.vs")
        rhopath = joinpath(filedir,"crust1.rho")
        bndpath = joinpath(filedir,"crust1.bnds")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(referencepath)
            print("Downloading crust1 files from google cloud storage to $filedir\n")
            run(`mkdir -p $filedir`)
            download("https://storage.googleapis.com/statgeochem/crust1.references.txt", referencepath)
            download("https://storage.googleapis.com/statgeochem/crust1.vp", vppath)
            download("https://storage.googleapis.com/statgeochem/crust1.vs", vspath)
            download("https://storage.googleapis.com/statgeochem/crust1.rho", rhopath)
            download("https://storage.googleapis.com/statgeochem/crust1.bnds", bndpath)
        end

        return 0 # Success
    end
    export get_crust1

    # Get all point data (Vp, Vs, Rho, layer thickness) from Crust 1.0 layer
    function find_crust1_layer(lat,lon,layer)
        # Get Vp, Vs, Rho, and thickness for a given lat, lon, and crustal layer.

        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        if ~isa(layer,Integer) || layer < 1 || layer > 8
            error("""Error: layer must be an integer between 1 and 8.
            Available layers:
            1) water
            2) ice
            3) upper sediments   (VP, VS, rho not defined in all cells)
            4) middle sediments  "
            5) lower sediments   "
            6) upper crystalline crust
            7) middle crystalline crust
            8) lower crystalline crust
            Results are returned in form (Vp, Vs, Rho, thickness)
            """)
        end

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        vp = Array{Float64,3}(undef,nlayers,nlat,nlon)
        vs = Array{Float64,3}(undef,nlayers,nlat,nlon)
        rho = Array{Float64,3}(undef,nlayers,nlat,nlon)
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
        vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
        rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              vp[:,j,i] = delim_string_parse(readline(vpfile), ' ', Float64, merge=true)
              vs[:,j,i] = delim_string_parse(readline(vsfile), ' ', Float64, merge=true)
              rho[:,j,i] = delim_string_parse(readline(rhofile), ' ', Float64, merge=true)
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(vpfile)
        close(vsfile)
        close(rhofile)
        close(bndfile)

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon .+ 180, 360) .- 180
        ilat = max.(lat, -90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        vpout = Array{Float64}(undef,size(lat))
        vsout = Array{Float64}(undef,size(lat))
        rhoout = Array{Float64}(undef,size(lat))
        thkout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        for j=1:length(lat)
            if isnan(lat[j]) || isnan(lon[j]) || lat[j] > 90 || lat[j] < -90 || lon[j] > 180 || lat[j] < -180
                vpout[j] = NaN
                vsout[j] = NaN
                rhoout[j] = NaN
                thkout[j] = NaN
            else
                vpout[j] = vp[layer,ilat[j],ilon[j]]
                vsout[j] = vs[layer,ilat[j],ilon[j]]
                rhoout[j] = rho[layer,ilat[j],ilon[j]]
                thkout[j] = bnd[layer,ilat[j],ilon[j]] - bnd[layer+1,ilat[j],ilon[j]]
            end
        end

        # The end
        return (vpout, vsout, rhoout, thkout)
    end
    export find_crust1_layer

    # Get seismic data (Vp, Vs, Rho) for crust 1.0 layer
    function find_crust1_seismic(lat,lon,layer)
        # Get Vp, Vs, and Rho for a given lat, lon, and crustal layer.

        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        if ~isa(layer,Integer) || layer < 1 || layer > 9
            error("""Error: layer must be an integer between 1 and 9.
            Available layers:
            1) water
            2) ice
            3) upper sediments   (VP, VS, rho not defined in all cells)
            4) middle sediments  "
            5) lower sediments   "
            6) upper crystalline crust
            7) middle crystalline crust
            8) lower crystalline crust
            9) Top of mantle below crust
            Results are returned in form (Vp, Vs, Rho)
            """)
        end

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        vp = Array{Float64,3}(undef,nlayers,nlat,nlon)
        vs = Array{Float64,3}(undef,nlayers,nlat,nlon)
        rho = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
        vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
        rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")

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

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon .+ 180, 360) .- 180
        ilat = max.(lat, -90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        vpout = Array{Float64}(undef,size(lat))
        vsout = Array{Float64}(undef,size(lat))
        rhoout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        for j=1:length(lat)
            if isnan(lat[j]) || isnan(lon[j]) || lat[j] > 90 || lat[j] < -90 || lon[j] > 180 || lat[j] < -180
                vpout[j] = NaN
                vsout[j] = NaN
                rhoout[j] = NaN
            else
                vpout[j] = vp[layer,ilat[j],ilon[j]]
                vsout[j] = vs[layer,ilat[j],ilon[j]]
                rhoout[j] = rho[layer,ilat[j],ilon[j]]
            end
        end

        # The end
        return (vpout, vsout, rhoout)
    end
    export find_crust1_seismic

    # Get layer thickness for crust 1.0 layer
    function find_crust1_thickness(lat,lon,layer)
        # Layer thickness for a given lat, lon, and crustal layer.

        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        if ~isa(layer,Integer) || layer < 1 || layer > 8
            error("""Error: layer must be an integer between 1 and 8.
            Available layers:
            1) water
            2) ice
            3) upper sediments   (VP, VS, rho not defined in all cells)
            4) middle sediments  "
            5) lower sediments   "
            6) upper crystalline crust
            7) middle crystalline crust
            8) lower crystalline crust
            Result is thickness of the requested layer
            """)
        end

        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(bndfile)

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon+180, 360) - 180
        ilat = max.(lat,-90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        thkout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        for j=1:length(lat)
            if isnan(lat[j]) || isnan(lon[j]) || lat[j] > 90 || lat[j] < -90 || lon[j] > 180 || lat[j] < -180
                thkout[j] = NaN
            else
                thkout[j] = bnd[layer,ilat[j],ilon[j]]-bnd[layer+1,ilat[j],ilon[j]]
            end
        end

        # The end
        return thkout
    end
    export find_crust1_thickness

    # Get detph to layer base for crust 1.0 layer
    function find_crust1_base(lat,lon,layer)
        # Layer thickness for a given lat, lon, and crustal layer.

        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        if ~isa(layer,Integer) || layer < 1 || layer > 8
            error("""layer must be an integer between 1 and 8.
            Available layers:
            1) water
            2) ice
            3) upper sediments   (VP, VS, rho not defined in all cells)
            4) middle sediments  "
            5) lower sediments   "
            6) upper crystalline crust
            7) middle crystalline crust
            8) lower crystalline crust
            Result is depth from sea level to base of the requested layer
            """)
        end
        nlayers=9
        nlon=360
        nlat=180

        # Allocate data arrays
        bnd = Array{Float64,3}(undef,nlayers,nlat,nlon)

        # Open data files
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = delim_string_parse(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(bndfile)

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon+180, 360) - 180
        ilat = max.(lat,-90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        baseout = Array{Float64}(undef,size(lat))

        # Fill output arrays
        for j=1:length(lat)
            if isnan(lat[j]) || isnan(lon[j]) || lat[j] > 90 || lat[j] < -90 || lon[j] > 180 || lat[j] < -180
                baseout[j] = NaN
            else
                baseout[j] = bnd[layer+1,ilat[j],ilon[j]]
            end
        end

        # The end
        return baseout
    end
    export find_crust1_base

## --- ETOPO1 (1 arc minute topography)

    # Read etopo file from HDF5 storage, downloading from cloud if necessary
    function get_etopo(varname="")
        # Available variable names: "elevation", "y_lat_cntr", "x_lon_cntr",
        # "cellsize", "scalefactor", and "reference". Units are meters of
        # elevation and decimal degrees of latitude and longitude

        # Construct file path
        filedir = joinpath(resourcepath,"etopo")
        filepath = joinpath(filedir,"etopo1.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading etopo1.h5 from google cloud storage to $filedir\n")
            run(`mkdir -p $filedir`)
            download("https://storage.googleapis.com/statgeochem/etopo1.references.txt", joinpath(filedir,"etopo1.references.txt"))
            download("https://storage.googleapis.com/statgeochem/etopo1.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath, "vars/"*varname)
    end
    export get_etopo

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the ETOPO elevation model.
    function find_etopoelev(etopo,lat,lon)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(etopo,Dict)
            data = etopo["elevation"]
        elseif isa(etopo, Array)
            data = etopo
        else
            error("wrong etopo variable")
        end


        # Scale factor (cells per degree) = 60 = arc minutes in an arc degree
        sf = 60
        maxrow = 180 * sf
        maxcol = 360 * sf

        # Create and fill output vector
        out = Array{Float64}(undef,size(lat))
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
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
                out[i] = data[row,col]
            end
        end

        return out
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
            print("Downloading srtm15plus.h5 from google cloud storage to $filedir\n")
            run(`mkdir -p $filedir`)
            download("https://storage.googleapis.com/statgeochem/srtm15plus.references.txt", joinpath(filedir,"srtm15plus.references.txt"))
            download("https://storage.googleapis.com/statgeochem/srtm15plus.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/"*varname)
    end
    export get_srtm15plus

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the SRTM15plus 15-arc-second elevation model.
    function find_srtm15plus(srtm15plus,lat,lon)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(srtm15plus,Dict)
            data = srtm15plus["elevation"]
        elseif isa(srtm15plus, Array)
            data = srtm15plus
        else
            error("wrong srtm15plus variable")
        end

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out = Array{Float64}(undef,size(lat))
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
                res = data[row,col]
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

## --- Müller et al. seafloor age and spreading rate

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
