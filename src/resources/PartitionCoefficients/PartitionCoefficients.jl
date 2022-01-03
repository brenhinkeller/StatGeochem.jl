function claiborne_zircon_kd(elem, T)
    # Convert temperature to Kelvin
    T += 273.15
    # Calculate partition coefficient
    if elem=="Hf"
        0.0965*exp(10206/T)
    elseif elem=="Th"
        0.0126*exp(6696/T)
    elseif elem=="U"
        0.0465*exp(7167/T)
    elseif elem=="Y"
        0.0036*exp(9806/T)
    elseif elem=="Nb"
        0.0003*exp(7241/T)
    elseif elem=="Nd"
        0.0001*exp(5867/T)
    elseif elem=="Sm"
        0.0009*exp(6636/T)
    elseif elem=="Tb"
        0.0021*exp(9160/T)
    elseif elem=="Eu"
        0.0032*exp(6026/T)
    elseif elem=="Dy"
        0.0041*exp(9090/T)
    elseif elem=="Gd"
        0.0005*exp(9436/T)
    elseif elem=="Ho"
        0.0038*exp(9948/T)
    elseif elem=="Er"
        0.0052*exp(10088/T)
    elseif elem=="Yb"
        0.0044*exp(10784/T)
    elseif elem=="Tm"
        0.0086*exp(9990/T)
    elseif elem=="Lu"
        0.0032*exp(11358/T)
    else
        NaN
    end
end
export claiborne_zircon_kd

export germ_kd