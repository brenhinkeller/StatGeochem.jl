## --- Resize and interpolate colormaps

    # Linearly interpolate cmap at positions xq
    function linterp1(x::AbstractArray, cmap::AbstractArray{<:Color}, xq)
        # Interpolate red, green, and blue vectors separately
        r_interp = linterp1(x, cmap .|> c -> c.r, xq)
        g_interp = linterp1(x, cmap .|> c -> c.g, xq)
        b_interp = linterp1(x, cmap .|> c -> c.b, xq)
        # Convert back to a color
        return RGB.(r_interp,g_interp,b_interp)
    end

    function resize_colormap(cmap::AbstractArray{<:Color}, n::Integer)
        cNum = length(cmap)
        if n<2
            cmap[1:1]
        else
            linterp1(1:cNum,cmap,collect(range(1,cNum,length=n)))
        end
    end
    export resize_colormap


## --- Map colormaps to images

    # Convert matrix to image using colormap
    function imsc(matrix::AbstractArray,colormap::AbstractArray=viridis,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap)
        if cmin>=cmax
            cmin = nanminimum(matrix)
            cmax = nanmaximum(matrix)
        end
        crange = cmax - cmin
        return  matrix .|> x -> colormap[isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc))]
    end
    export imsc

    # Convert matrix to indirect array image using colormap
    function imsci(matrix::AbstractArray,colormap::AbstractArray=viridis,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap)
        if cmin>=cmax
            cmin = nanminimum(matrix)
            cmax = nanmaximum(matrix)
        end
        crange = cmax - cmin
        return IndirectArray(matrix .|> x -> isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc)), colormap)
    end
    export imsci


## -- End of File
