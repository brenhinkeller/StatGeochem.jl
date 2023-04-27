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

"""
```julia
find_litho1_property(lat, lon, layer::Symbol, property::Symbol)
```
Return values for a LITHO1.0 `property` of a given `layer` at one or more given
`lat`itudes and `lon`gitudes, to the nearest 0.5-arc-degree grid point.

Accepts `lat` and `lon` both as `Numbers` and as `AbstractArray`s, but given
the overhead of opening and reading the LITHO1.0 files, you should generally
aim to provide large arrays with as many values in a single query as possible.

Available properties include:
```
    :vp             | p-wave velocity [m/s]
    :vs             | s-wave velocity [m/s]
    :rho            | density [kg/m^3]
    :bottom         | depth to bottom of the layer [km] (above sea level = negative)
    :thickness      | layer thickness [km]
```
while avialble `layer`s are:
```
    1 | :ice
    2 | :water
    3 | :upper_sediments
    4 | :middle_sediments
    5 | :lower_sediments
    6 | :upper_crust
    7 | :middle_crust
    8 | :lower_crust
    9 | :sclm (or :lithosphere)
    10 | :asthenosphere
```

## Examples
```julia
julia> find_litho1_property([43.702245, 44], [-72.0929, -73], :upper_crust, :vp)
2-element Vector{Float64}:
 6219.99
 6253.16
```
"""
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
