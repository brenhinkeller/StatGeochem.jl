## -- Establish path
resourcepath = joinpath(Pkg.dir("StatGeochem"),"resources");
export resourcepath


# Get points from Crust 1.0.
function get_crust1_point(lat,lon,layer)
    # Get Vp, Vs, Rho, and thickness for a given lat, lon, and crustal layer.
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

    np=9;
    nlo=360;
    nla=180;

    # Allocate data arrays
    vp = Array{Float64,3}(np,nla,nlo)
    vs = Array{Float64,3}(np,nla,nlo)
    rho = Array{Float64,3}(np,nla,nlo)
    bnd = Array{Float64,3}(np,nla,nlo)

    # Open data files
    vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
    vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
    rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")
    bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

    # Read data files into array
    for j=1:nla
       for i=1:nlo
          vp[:,j,i] = parseDelimString(readline(vpfile), ' ', Float64, merge=true)
          vs[:,j,i] = parseDelimString(readline(vsfile), ' ', Float64, merge=true)
          rho[:,j,i] = parseDelimString(readline(rhofile), ' ', Float64, merge=true)
          bnd[:,j,i] = parseDelimString(readline(bndfile), ' ', Float64, merge=true)
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
        if lat[j] <= 90 && lat[j]>=-90 && lon[j] <= 180 && lat[j]>=-180
            vpout[j] = vp[layer,ilat[j],ilon[j]]
            vsout[j] = vs[layer,ilat[j],ilon[j]]
            rhoout[j] = rho[layer,ilat[j],ilon[j]]
            thkout[j] = bnd[layer,ilat[j],ilon[j]]-bnd[layer+1,ilat[j],ilon[j]]
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
export get_crust1_point

function get_crust1_seismic(lat,lon,layer)
    # Get Vp, Vs, and Rho for a given lat, lon, and crustal layer.
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

    np=9;
    nlo=360;
    nla=180;

    # Allocate data arrays
    vp = Array{Float64,3}(np,nla,nlo)
    vs = Array{Float64,3}(np,nla,nlo)
    rho = Array{Float64,3}(np,nla,nlo)

    # Open data files
    vpfile = open(joinpath(resourcepath,"crust1","crust1.vp"), "r")
    vsfile = open(joinpath(resourcepath,"crust1","crust1.vs"), "r")
    rhofile = open(joinpath(resourcepath,"crust1","crust1.rho"), "r")

    # Read data files into array
    for j=1:nla
       for i=1:nlo
          vp[:,j,i] = parseDelimString(readline(vpfile), ' ', Float64, merge=true)
          vs[:,j,i] = parseDelimString(readline(vsfile), ' ', Float64, merge=true)
          rho[:,j,i] = parseDelimString(readline(rhofile), ' ', Float64, merge=true) * 1000 # convert to kg/m3
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
        if lat[j] <= 90 && lat[j]>=-90 && lon[j] <= 180 && lat[j]>=-180
            vpout[j] = vp[layer,ilat[j],ilon[j]]
            vsout[j] = vs[layer,ilat[j],ilon[j]]
            rhoout[j] = rho[layer,ilat[j],ilon[j]]
        else
            vpout[j] = NaN
            vsout[j] = NaN
            rhoout[j] = NaN
        end
    end

    # The end
    return (vpout, vsout, rhoout)
end
export get_crust1_seismic

function get_crust1_thickness(lat,lon,layer)
    # Layer thickness for a given lat, lon, and crustal layer.
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

    np=9;
    nlo=360;
    nla=180;

    # Allocate data arrays
    bnd = Array{Float64,3}(np,nla,nlo)

    # Open data files
    bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

    # Read data files into array
    for j=1:nla
       for i=1:nlo
          bnd[:,j,i] = parseDelimString(readline(bndfile), ' ', Float64, merge=true)
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
        if lat[j] <= 90 && lat[j]>=-90 && lon[j] <= 180 && lat[j]>=-180
            thkout[j] = bnd[layer,ilat[j],ilon[j]]-bnd[layer+1,ilat[j],ilon[j]]
        else
            thkout[j] = NaN
        end
    end

    # The end
    return thkout
end
export get_crust1_thickness

function get_crust1_base(lat,lon,layer)
    # Layer thickness for a given lat, lon, and crustal layer.
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
        Result is depth to base of the requested layer
        """)
    end
    np=9;
    nlo=360;
    nla=180;

    # Allocate data arrays
    bnd = Array{Float64,3}(np,nla,nlo)

    # Open data files
    bndfile = open(joinpath(resourcepath,"crust1","crust1.bnds"), "r")

    # Read data files into array
    for j=1:nla
       for i=1:nlo
          bnd[:,j,i] = parseDelimString(readline(bndfile), ' ', Float64, merge=true)
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
        if lat[j] <= 90 && lat[j]>=-90 && lon[j] <= 180 && lat[j]>=-180
            baseout[j] = bnd[layer+1,ilat[j],ilon[j]]
        else
            baseout[j] = NaN
        end
    end

    # The end
    return baseout
end
export get_crust1_base
