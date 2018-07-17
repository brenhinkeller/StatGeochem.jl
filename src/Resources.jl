## -- Establish path
    resourcepath = joinpath(Pkg.dir("StatGeochem"),"resources");
    export resourcepath

## --- CRUST 1.0

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

        nlayers=9;
        nlon=360;
        nlat=180;

        # Allocate data arrays
        vp = Array{Float64,3}(nlayers,nlat,nlon)
        vs = Array{Float64,3}(nlayers,nlat,nlon)
        rho = Array{Float64,3}(nlayers,nlat,nlon)
        bnd = Array{Float64,3}(nlayers,nlat,nlon)

        # Open data files
        vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
        vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
        rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              vp[:,j,i] = parse_delim_string(readline(vpfile), ' ', Float64, merge=true)
              vs[:,j,i] = parse_delim_string(readline(vsfile), ' ', Float64, merge=true)
              rho[:,j,i] = parse_delim_string(readline(rhofile), ' ', Float64, merge=true)
              bnd[:,j,i] = parse_delim_string(readline(bndfile), ' ', Float64, merge=true)
          end
        end

        # Close data files
        close(vpfile)
        close(vsfile)
        close(rhofile)
        close(bndfile)

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon+180, 360) - 180
        ilat = max.(lat,-90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        vpout = Array{Float64}(size(lat));
        vsout = Array{Float64}(size(lat));
        rhoout = Array{Float64}(size(lat));
        thkout = Array{Float64}(size(lat));

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
                thkout[j] = bnd[layer,ilat[j],ilon[j]]-bnd[layer+1,ilat[j],ilon[j]]
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

        nlayers=9;
        nlon=360;
        nlat=180;

        # Allocate data arrays
        vp = Array{Float64,3}(nlayers,nlat,nlon)
        vs = Array{Float64,3}(nlayers,nlat,nlon)
        rho = Array{Float64,3}(nlayers,nlat,nlon)

        # Open data files
        vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
        vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
        rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              vp[:,j,i] = parse_delim_string(readline(vpfile), ' ', Float64, merge=true)
              vs[:,j,i] = parse_delim_string(readline(vsfile), ' ', Float64, merge=true)
              rho[:,j,i] = parse_delim_string(readline(rhofile), ' ', Float64, merge=true) * 1000 # convert to kg/m3
          end
        end

        # Close data files
        close(vpfile)
        close(vsfile)
        close(rhofile)

        # Avoid edge cases at lat = -90.0, lon = 180.0
        ilon = mod.(lon+180, 360) - 180
        ilat = max.(lat,-90+1e-9)

        # Convert lat and lon to index
        ilat = 90 - ceil.(Int,ilat) + 1
        ilon = 180 + floor.(Int,ilon) + 1

        # Allocate output arrays
        vpout = Array{Float64}(size(lat));
        vsout = Array{Float64}(size(lat));
        rhoout = Array{Float64}(size(lat));

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

        nlayers=9;
        nlon=360;
        nlat=180;

        # Allocate data arrays
        bnd = Array{Float64,3}(nlayers,nlat,nlon)

        # Open data files
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = parse_delim_string(readline(bndfile), ' ', Float64, merge=true)
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
        thkout = Array{Float64}(size(lat));

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
        nlayers=9;
        nlon=360;
        nlat=180;

        # Allocate data arrays
        bnd = Array{Float64,3}(nlayers,nlat,nlon)

        # Open data files
        bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

        # Read data files into array
        for j=1:nlat
           for i=1:nlon
              bnd[:,j,i] = parse_delim_string(readline(bndfile), ' ', Float64, merge=true)
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
        baseout = Array{Float64}(size(lat));

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

    # Read etopoelev file from HDF5 storage, downloading from cloud if necessary
    function get_etopoelev()
        # Construct file path
        filepath = joinpath(resourcepath,"etopo","etopoelev.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading etopoelev.h5 from google cloud storage\n")
            download("https://storage.googleapis.com/statgeochem/etopoelev.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/etopoelev")
    end
    export get_etopoelev

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the ETOPO elevation model.
    function find_etopoelev(etopoelev,lat,lon)
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Scale factor = 60 = arc minutes in an arc degree
        sf = 60
        maxrow = 180 * sf
        maxcol = 360 * sf

        # Create and fill output vector
        elev=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                elev[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                row = 1 + trunc(Int,(90+lat[i])*sf)
                if row == (maxrow+1) # Edge case
                    row = maxrow;
                end

                col = 1 + trunc(Int,(180+lon[i])*sf)
                if col == (maxcol+1) # Edge case
                    col = maxcol;
                end

                # Find result by indexing
                elev[i] = etopoelev[row,col]
            end
        end

        return elev
    end
    export find_etopoelev

## --- SRTM15_PLUS (15 arc second topography)

    # Read srtm15plus file from HDF5 storage, downloading from cloud if necessary
    function get_srtm15plus()
        # Construct file path
        filepath = joinpath(resourcepath,"SRTM15plus","SRTM15plus.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading SRTM15plus.h5 from google cloud storage\n")
            download("https://storage.googleapis.com/statgeochem/SRTM15plus.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/elev")
    end
    export get_srtm15plus

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the SRTM15plus 15-arc-second elevation model.
    function find_srtm15plus(srtm15plus,lat,lon)
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Scale factor = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        elev=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                elev[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90-lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)
                # Find result by indexing
                res = srtm15plus[row,col]
                if res == -32768
                    elev[i] = NaN
                else
                    elev[i] = res
                end
            end
        end

        return elev
    end
    export find_srtm15plus

## --- Müller et al. seafloor age and spreading rate

    # Read seafloorage file from HDF5 storage, downloading from cloud if necessary
    function get_seafloorage()
        # Construct file path
        filepath = joinpath(resourcepath,"seafloorage","seafloorage.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading seafloorage.h5 from google cloud storage\n")
            download("https://storage.googleapis.com/statgeochem/seafloorage.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/seafloorage")
    end
    export get_seafloorage

    function get_seafloorage_sigma()
        # Construct file path
        filepath = joinpath(resourcepath,"seafloorage","seafloorage.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading seafloorage.h5 from google cloud storage\n")
            download("https://storage.googleapis.com/statgeochem/seafloorage.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/seafloorage_sigma")
    end
    export get_seafloorage_sigma

    # Read seafloorrate file from HDF5 storage, downloading from cloud if necessary
    function get_seafloorrate()
        # Construct file path
        filepath = joinpath(resourcepath,"seafloorage","seafloorrate.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading seafloorrate.h5 from google cloud storrate\n")
            download("https://storrate.googleapis.com/statgeochem/seafloorrate.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath,"vars/seafloorrate")
    end
    export get_seafloorrate

    # Parse seafloorage, seafloorage_sigma, or seafloorrate from file
    # data = find_seafloorage(sfdata, lat, lon)
    function find_seafloorage(sfdata,lat,lon)
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Find the column numbers (using mod to convert lon from -180:180 to 0:360
        x = floor.(Int, mod.(lon, 360) * 10800/360) + 1

        # find the y rows, converting from lat to Mercator (lat -80.738:80.738)
        y = 4320 - floor.(Int, 8640 * asinh.(tan.(lat*pi/180)) / asinh.(tan.(80.738*pi/180)) / 2 ) + 1

        # Make and fill output array
        out=Array{Float64}(size(x));
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
