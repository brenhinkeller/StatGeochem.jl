## -- Establish path
resourcepath = joinpath(Pkg.dir("StatGeochem"),"resources");
export resourcepath

function getCN1point(lat,lon,layer)
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

    # make sure longitudes go from -180 to 180
    ilon = mod.(lon+180, 360) - 180
    ilat = mod.(lat+90, 180) - 90

    # Convert lat and lon to index
    ilat = 90 - floor(Int,ilat) + 1
    ilon = 180 + floor(Int,ilon) + 1

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
export getCN1point
