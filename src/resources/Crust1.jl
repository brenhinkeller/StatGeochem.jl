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
        ilat = 91 .- ceil.(Int,ilat)
        ilon = 181 .+ floor.(Int,ilon)

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
        ilat = 91 .- ceil.(Int,ilat)
        ilon = 181 .+ floor.(Int,ilon)

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
        ilon = mod.(lon .+ 180, 360) .- 180
        ilat = max.(lat, -90+1e-9)

        # Convert lat and lon to index
        ilat = 91 .- ceil.(Int,ilat)
        ilon = 181 .+ floor.(Int,ilon)

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
        ilon = mod.(lon .+ 180, 360) .- 180
        ilat = max.(lat, -90+1e-9)

        # Convert lat and lon to index
        ilat = 91 .- ceil.(Int,ilat)
        ilon = 181 .+ floor.(Int,ilon)

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

## --- End of File
