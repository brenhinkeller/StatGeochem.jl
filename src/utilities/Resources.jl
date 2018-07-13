## --- Find ETOPO elevation

# load data as:
#
# using MAT
# include("Utilities.jl")
# etopoelev = matread("/Users/cbkeller/Documents/MATLAB/statistical geochemistry/resources/etopo/etopoelev.mat")["etopoelev"]

function findelevation(etopoelev,lat,lon)
    # elev=findelevation(lat,lon,[etopoelev])
    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the ETOPO elevation model.

    # Scale factor used in map (=nrows/180=ncols/360)
    sf=60;
    maxrow = 180*sf;
    maxcol = 360*sf;

    # Create and fill output vector
    elev=Array{Float64}(length(lat));
    for i=1:length(lat)
        if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
            elev[i]=NaN; # Result is NaN if either input is NaN
        else
            # Convert latitude and longitude into indicies of the elevation map array
            row=round(Int,(90+lat[i])*sf+0.5);
            if row == (maxrow+1)
                row = maxrow;
            end

            col=round(Int,(180+lon[i])*sf+0.5);
            if col == (maxcol+1)
                col = maxcol;
            end

            elev[i]=etopoelev[row,col]; # Otherwise, find result
        end
    end

    return elev
end
export findelevation
