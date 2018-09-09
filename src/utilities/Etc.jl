## --- Digitize image data

    # Calculate approximate x and y positions and uncertainties for colored
    # markers in an image
    function digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)
        # (x,dx,y,dy) = digitize_plotmarkers(img, marker_color, xlims, ylims; atol=0.16)

        # Test for approximate equality in color to marker
        t = isapprox.(marker_color,img,atol=atol)

        # Figure out our image dimensions
        xrows = size(t,2)
        yrows = size(t,1)
        xlims = collect(xlims)
        ylims = collect(ylims)

        # Allocate index arrays
        imin = Array{Float64}(round(Int,xrows/2))
        imax = Array{Float64}(round(Int,xrows/2))
        jmin = Array{Float64}(round(Int,xrows/2))
        jmax = Array{Float64}(round(Int,xrows/2))

        # Fill index arrays
        found = false;
        markernumber = 0
        for j=1:xrows
            if sum(t[:,j])>0
                list = find(t[:,j])
                if ~found
                    markernumber += 1;
                    imin[markernumber] = minimum(list)
                    imax[markernumber] = maximum(list)
                    jmin[markernumber] = j;
                    jmax[markernumber] = j;
                else
                    imin[markernumber] = min(imin[markernumber],minimum(list))
                    imax[markernumber] = max(imax[markernumber],maximum(list))
                    jmax[markernumber] = j;
                end
                found = true;
            else
                found = false;
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

        return (x,dx,y,dy)
    end

## --- End of File
