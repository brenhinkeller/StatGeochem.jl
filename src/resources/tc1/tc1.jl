## --- Query TC1 properties

    tc1_550 = readdlm(joinpath(moduleresourcepath,"tc1_550.csv"), ',')
    tc1_1300 = readdlm(joinpath(moduleresourcepath,"tc1_1300.csv"), ',')
    tc1_age = readdlm(joinpath(moduleresourcepath,"tc1_age.csv"), ',', Int)

    function find_tc1_crust(lat::Number,lon::Number)
        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            crust=tc1_550[i,j]
        else
            crust = NaN
        end
        return crust
    end
    function find_tc1_crust(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Query the tc1_550 array for our lat and lon
        crust = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                crust[n]=tc1_550[i,j]
            end
        end
        return crust
    end
    export find_tc1_crust

    function find_tc1_lith(lat::Number,lon::Number)
        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            lith=tc1_1300[i,j]
        else
            lith = NaN
        end
        return lith
    end
    function find_tc1_lith(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        # Query the tc1_1300 array for our lat and lon
        lith = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                lith[n]=tc1_1300[i,j]
            end
        end
        return lith
    end
    export find_tc1_lith

    function find_tc1_age(lat::Number,lon::Number)
        ages=[ NaN  NaN  NaN
                25    0   50
               150   50  250
               395  250  540
               695  540  850
               975  850 1100
              1400 1100 1700
              2100 1700 2500
              2750 2500 3000
              3250 3000 3500]

        if !isnan(lat) && !isnan(lon)
            i = round(Int, 91-lat)
            j = round(Int, lon+181)
            result = (ages[tc1_age[i,j],:]...,)
        else
            result = (NaN, NaN, NaN,)
        end
        return result
    end
    function find_tc1_age(lat::AbstractArray,lon::AbstractArray)
        # Check input dimensions match
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        end

        ages=[ NaN  NaN  NaN
                25    0   50
               150   50  250
               395  250  540
               695  540  850
               975  850 1100
              1400 1100 1700
              2100 1700 2500
              2750 2500 3000
              3250 3000 3500]

        # Query the tc1_1300 array for our lat and lon
        age = fill(NaN, size(lat))
        minage = fill(NaN, size(lat))
        maxage = fill(NaN, size(lat))
        for n=1:length(lat)
            if !isnan(lat[n]) && !isnan(lon[n])
                i = round(Int, 91-lat[n])
                j = round(Int, lon[n]+181)
                ageindex = tc1_age[i,j]
                age[n] = ages[ageindex,1]
                minage[n] = ages[ageindex,2]
                maxage[n] = ages[ageindex,3]
            end
        end
        return (age, minage, maxage)
    end
    export find_tc1_age

## --- End of File
