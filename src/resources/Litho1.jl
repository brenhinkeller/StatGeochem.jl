litho1layer(layer::Number) = Int(layer)
function litho1layer(layer::Symbol)
    layernumber = 0
    if layer===:ice
        1
    elseif layer===:water
        2
    elseif layer===:upper_sediments
        3
    elseif layer===:middle_sediments
        4
    elseif layer===:lower_sediments
        5
    elseif layer===:upper_crust
        6
    elseif layer===:middle_crust
        7
    elseif layer===:lower_crust
        8
    elseif layer===:lithosphere || layer===:sclm || layer===:lid
        9
    elseif layer===:asthenosphere
        10
    else
        @warn "litho1 layer $layer not found"
        0
    end
end

function find_litho1_property(lat, lon, layer, property::Symbol)
    @assert eachindex(lat)==eachindex(lon)
    layerindex = litho1layer(layer)::Int
    @assert 0 < layerindex < 11
    litho1path = artifact"litho1-gridded/litho1-gridded.h5"

    result = fill(NaN, size(lat))
    if property===:vp
        data = h5read(litho1path, "vp")::Array{Float64,3}
    elseif property===:vs
        data = h5read(litho1path, "vs")::Array{Float64,3}
    elseif property===:rho || property===:density
        data = h5read(litho1path, "rho")::Array{Float64,3}
    elseif property===:base || property===:bottom
        data = h5read(litho1path, "bottom")::Array{Float64,3}
    elseif property===:thickness
        data = h5read(litho1path, "thickness")::Array{Float64,3}
    else
        @warn """litho1 property `$property` not found. Available options include:
                    :vp, :vs, :rho, :bottom, :thickness
        """
        return result
    end
    # For reference:
    # lats = -90:0.5:90
    # lons = -179.5:0.5:180
    scale = 2
    @inbounds for i in eachindex(lat)
        r = scale*(lat[i]+90) + 1
        c = scale*(lon[i]+180)
        if 1 <= r <= size(data,1) && 0 <= c <= size(data,2)
            row, col = trunc(Int, r), trunc(Int, c)
            col < 1 && (col = size(data,2))
            result[i] = data[row, col, layerindex]
        end
    end
    return result
end
export find_litho1_property
