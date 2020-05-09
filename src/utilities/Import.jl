## --- Parse a delimited string

    # Parse a delimited string, return the results in a pre-allocated array provided as input
    function delim_string_parse!(result::Array, str::AbstractString, delim::Char, parseType::Type=eltype(result); offset::Integer=0, merge::Bool=false, undefval=NaN)

        # Make sure the output data type allows our chosen value for undefined data
        undefval = convert(parseType, undefval)

        # Ignore initial delimiter
        last_delim_pos = 0
        if ~isempty(str) && str[1] == delim
            last_delim_pos = 1
        end

        # Cycle through string parsing text betweeen delims
        delim_pos = 0
        n = offset
        if merge
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos+1
                        n += 1
                        parsed = nothing
                        if delim_pos > last_delim_pos+1
                            parsed = tryparse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? undefval : parsed
                    end
                    last_delim_pos = delim_pos
                end
            end
        else
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos
                        n += 1
                        parsed = nothing
                        if delim_pos > last_delim_pos+1
                            parsed = tryparse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? undefval : parsed
                        last_delim_pos = delim_pos
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str) > last_delim_pos
            n += 1
            parsed = tryparse(parseType, str[(last_delim_pos+1):length(str)])
            result[n] = isnothing(parsed) ? undefval : parsed
        end

        # Return the number of result values
        return n-offset
    end
    export delim_string_parse!

    # Parse a delimited string, return an array as output
    function delim_string_parse(str::AbstractString, delim::Char, parseType::Type=Float64; merge::Bool=false, undefval=NaN)

        # Allocate an array to hold our parsed results
        result = Array{parseType}(undef,ceil(Int,length(str)/2))

        # Make sure the output data type allows our chosen value for undefined data
        undefval = convert(parseType, undefval)

        # Ignore initial delimiter
        last_delim_pos = 0
        if ~isempty(str) && str[1] == delim
            last_delim_pos = 1
        end

        # Cycle through string parsing text betweeen delims
        delim_pos = 0
        n = 0
        if merge
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos+1
                        n += 1
                        parsed = nothing
                        if delim_pos > last_delim_pos+1
                            parsed = tryparse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? undefval : parsed
                    end
                    last_delim_pos = delim_pos
                end
            end
        else
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos
                        n += 1
                        parsed = nothing
                        if delim_pos>last_delim_pos+1
                            parsed = tryparse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? undefval : parsed
                        last_delim_pos = delim_pos
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1
            parsed = tryparse(parseType, str[(last_delim_pos+1):length(str)])
            result[n] = isnothing(parsed) ? undefval : parsed
        end

        # Return the result values
        return result[1:n]
    end
    export delim_string_parse

    function delim_string_function(f::Function, str::AbstractString, delim::Char, outType::Type; merge::Bool=false)
        # result = delim_string_function(f::Function, str::AbstractString, delim::Char, outType::Type; merge=false)

        # Max number of delimted values
        ndelims = 2
        for i = 1:length(str)
            if str[i] == delim
                ndelims += 1
            end
        end

        # Allocate output array
        result = Array{outType}(undef,ceil(Int,ndelims))

        # Ignore initial delimiter
        last_delim_pos = 0
        if str[1] == delim
            last_delim_pos = 1
        end

        # Cycle through string parsing text betweeen delims
        delim_pos = 0
        n = 0
        if merge
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos+1
                        n += 1
                        if delim_pos > last_delim_pos+1
                            result[n] = f(str[(last_delim_pos+1):(delim_pos-1)])
                        end
                    end
                    last_delim_pos = delim_pos
                end
            end
        else
            for i = 1:length(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos
                        n += 1
                        if delim_pos > last_delim_pos+1
                            result[n] = f(str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        last_delim_pos = delim_pos
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1
            result[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the result values
        return result[1:n]
    end
    export delim_string_function

## --- Classifying imported datasets

    if VERSION >= v"0.7"
        # Return true for numbers and strings that can be parsed as numbers
        function plausiblynumeric(x)
            if isa(x,Number)
                return true
            elseif isa(x,AbstractString) && tryparse(Float64,x) != nothing
                return true
            else
                return false
            end
        end
    else
        # Return true for numbers and strings that can be parsed as numbers
        function plausiblynumeric(x)
            if isa(x,Number)
                return true
            elseif isa(x,AbstractString) && ~isnull(tryparse(Float64,x))
                return true
            else
                return false
            end
        end
    end
    export plausiblynumeric

    if VERSION >= v"0.7"
        # Return true for values that are not missing and cannot be parsed as numbers
        function nonnumeric(x)
            if isa(x,Number)
                return false
            elseif isa(x,AbstractString) && (tryparse(Float64,x) != nothing || x == "")
                return false
            else
                return true
            end
        end
    else
        # Return true for values that are not missing and cannot be parsed as numbers
        function nonnumeric(x)
            if isa(x,Number)
                return false
            elseif isa(x,AbstractString) && (~isnull(tryparse(Float64,x)) || x == "")
                return false
            else
                return true
            end
        end
    end
    export nonnumeric

## --- Transforming imported datasets

    if VERSION >= v"0.7"
        # Convert to a Float64 if possible, or a Float64 NaN if not.
        function floatify(x)
            if isa(x,Number)
                return Float64(x)
            elseif isa(x,AbstractString) && tryparse(Float64,x) != nothing
                return parse(Float64,x)
            else
                return NaN
            end
        end
    else
        # Convert to a Float64 if possible, or a Float64 NaN if not.
        function floatify(x)
            if isa(x,Number)
                return Float64(x)
            elseif isa(x,AbstractString) && ~isnull(tryparse(Float64,x))
                return parse(Float64,x)
            else
                return NaN
            end
        end
    end
    export floatify

    # Convert a flat array into a dict with each column as a variable
    function elementify(dataset::Array, elements::Array=dataset[1,:]; floatout::Bool=true, skipstart::Integer=1, skipnameless::Bool=true)
        # Output as dictionary
        result = Dict()
        if skipnameless
            result["elements"] = elements[elements .!= ""]
        else
            result["elements"] = elements
        end

        # Parse the input array, minus empty-named columns
        for i = 1:length(elements)
            thiscol = dataset[(1+skipstart):end,i]
            floatcol = floatout && ( sum(plausiblynumeric.(thiscol)) >= sum(nonnumeric.(thiscol)) )

            if haskey(result,elements[i])
                # If key already exists
                if floatcol || (floatout && (sum(plausiblynumeric.(result[elements[i]])) >= sum(nonnumeric.(result[elements[i]]))))
                    # If either this column or the existing one is plausibly numeric, average the two
                    result[elements[i]] = nanmean( hcat(floatify.(result[elements[i]]), floatify.(thiscol)), dim=2 )
                else
                    # If neither is plausibly numeric, just contatenate the columns and move on
                    result[elements[i]] = hcat(result[elements[i]], thiscol)
                end
            elseif floatcol
                # If column is numeric
                result[elements[i]] = floatify.(thiscol)
            else
                # If column is non-numeric
                result[elements[i]] = thiscol
            end
        end

        # Return only unique elements, since dictionary keys must be unique
        result["elements"] = unique(elements)

        return result
    end
    export elementify

    # Convert a dict into a flat array with variables as columns
    function unelementify(dataset::Dict, elements::Array=sort(collect(keys(dataset))); floatout::Bool=false, findnumeric::Bool=false, skipnan::Bool=false)

        # Find the elements in the input dict if they exist and aren't otherwise specified
        if any(elements .== "elements")
            elements = dataset["elements"]
        end

        # Figure out how many are numeric (if necessary), so we can export only
        # those if `findnumeric` is set
        if findnumeric
            is_numeric_element = Array{Bool}(undef,length(elements))
            for i = 1:length(elements)
                is_numeric_element = sum(plausiblynumeric.(dataset[elements[i]])) > sum(nonnumeric.(dataset[elements[i]]))
            end
            elements = elements[is_numeric_element]
        end

        # Generate output array
        if floatout
            # Allocate output Array{Float64}
            result = Array{Float64}(undef,length(dataset[elements[1]]),length(elements))

            # Parse the input dict. No column names if `floatout` is set
            for i = 1:length(elements)
                result[:,i] = floatify.(dataset[elements[i]])
            end
        else
            # Allocate output Array{Any}
            result = Array{Any}(undef,length(dataset[elements[1]])+1,length(elements))

            # Parse the input dict
            for i = 1:length(elements)
                # Column name goes in the first row, everything else after that
                result[1,i] = elements[i]
                result[2:end,i] .= dataset[elements[i]]

                # if `skipnan` is set, replace each NaN in the output array with
                # an empty string ("") such that it is empty when printed to file
                # with dlmwrite or similar
                if skipnan
                    for n = 2:length(result[:,i])
                        if isa(result[n,i], AbstractFloat) && isnan(result[n,i])
                            result[n,i] = ""
                        end
                    end
                end
            end
        end
        return result
    end
    export unelementify

## --- Renormalization of imported datasets

    function renormalize!(dataset::Dict, elements::Array=sort(collect(keys(dataset))); total::Number=1)
        current_sum = zeros(size(dataset[elements[1]]))
        for e in elements
            current_sum .+= dataset[e] .|> x -> isnan(x) ? 0 : x
        end
        current_sum[current_sum .== 0] .= NaN

        for e in elements
            dataset[e] .*= total ./ current_sum
        end
    end
    function renormalize!(A::Array{<:AbstractFloat}; dim::Number=0, total::Number=1)
        current_sum = nansum(A, dim=dim)
        A .*= total ./ current_sum
    end
    export renormalize!

## --- High-level import/export functions

    function importdataset(filepath::AbstractString, delim::AbstractChar; floatout::Bool=true, skipstart::Integer=0, skipnameless::Bool=true, mindefinedcolumns::Integer=0)
        data = readdlm(filepath, delim, skipstart=skipstart)
        if mindefinedcolumns > 0
            definedcolumns = vec(sum(.~ isempty.(data), dims=2))
            t = definedcolumns .>= mindefinedcolumns
            data = data[t,:]
        end
        return elementify(data, floatout=floatout, skipnameless=skipnameless)
    end
    export importdataset

    function exportdataset(dataset::Dict, filepath::AbstractString, delim::AbstractChar, skipnan::Bool=true)
        if skipnan
            result = writedlm(filepath, unelementify(dataset, skipnan=true), delim)
        else
            result = writedlm(filepath, unelementify(dataset), delim)
        end
        return result
    end
    export exportdataset

## --- End of File
