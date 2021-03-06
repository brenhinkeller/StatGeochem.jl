## --- Digitize image data

    """
    ```julia
    (x,dx,y,dy) = digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)
    ```
    Calculate approximate x and y positions and uncertainties for colored
    markers in an image
    """
    function digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)

        # Test for approximate equality in color to marker
        t = isapprox.(marker_color,img,atol=atol)

        # Figure out our image dimensions
        xrows = size(t,2)
        yrows = size(t,1)
        xlims = collect(xlims)
        ylims = collect(ylims)

        # Allocate index arrays
        imin = Array{Float64}(undef,round(Int,xrows/2))
        imax = Array{Float64}(undef,round(Int,xrows/2))
        jmin = Array{Float64}(undef,round(Int,xrows/2))
        jmax = Array{Float64}(undef,round(Int,xrows/2))

        # Fill index arrays
        found = false
        markernumber = 0
        for j=1:xrows
            if any(t[:,j])
                list = findall(t[:,j])
                if ~found
                    markernumber += 1
                    imin[markernumber] = minimum(list)
                    imax[markernumber] = maximum(list)
                    jmin[markernumber] = j
                    jmax[markernumber] = j
                else
                    imin[markernumber] = min(imin[markernumber],minimum(list))
                    imax[markernumber] = max(imax[markernumber],maximum(list))
                    jmax[markernumber] = j
                end
                found = true
            else
                found = false
            end
        end

        # Return only the filled indices
        imin = imin[1:markernumber]
        imax = imax[1:markernumber]
        jmin = jmin[1:markernumber]
        jmax = jmax[1:markernumber]

        # Calculate x and y positions from indices
        y = ylims[2] - (imin+imax)/2 * nanrange(ylims) / yrows
        dy = (imax-imin)/2 * nanrange(ylims) / yrows
        x = (jmin+jmax)/2 * nanrange(xlims) / xrows + xlims[1]
        dx = (jmax-jmin)/2 * nanrange(xlims) / xrows

        return x, dx, y, dy
    end
    export digitize_plotmarkers

    """
    ```julia
    (x,y) = digitize_plotline(img, line_color, xlims, ylims; atol=0.16)
    ```
    Calculate approximate x and y positions for a colored line in an image
    """
    function digitize_plotline(img, line_color, xlims, ylims; atol=0.16)
        # Test for approximate equality in color to marker
        t = isapprox.(line_color,img,atol=atol)

        # Figure out our image dimensions
        xrows = size(t,2)
        yrows = size(t,1)
        xlims = collect(xlims)
        ylims = collect(ylims)

        # Calculate x for each column
        x = cntr(range(xlims[1], xlims[2], length=xrows+1))

        # y as a function of i-position in image
        # (note: images are typically flipped)
        yᵢ(i) = ylims[2] - i * nanrange(ylims) / yrows

        # Calculate y, defaulting to NaN if no matches
        y = fill(NaN, xrows)
        for j = 1:xrows
            if any(t[:,j])
                list = findall(t[:,j])
                y[j] = yᵢ(nanmean(list))
            else
                y[j] = NaN
            end
        end

        return x, y
    end
    export digitize_plotline

## --- End of File
