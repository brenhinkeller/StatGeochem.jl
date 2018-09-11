## --- Parse a delimited string

    # Parse a delimited string, return the results in a pre-allocated array provided as input
    function delim_string_parse!(parsed::Array, str::AbstractString, delim::Char, parseType::Type; offset=0, merge=false)

        # Ignore initial delimiter
        last_delim_pos=0;
        if str[1]==delim
            last_delim_pos=1;
        end

        # Cycle through string parsing text betweeen delims
        delim_pos=0;
        n = offset;
        if merge
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos+1
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = parse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                    end
                    last_delim_pos = delim_pos;
                end
            end
        else
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = parse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        last_delim_pos = delim_pos;
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1;
            parsed[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the number of parsed values
        return n-offset
    end
    export delim_string_parse!

    # Parse a delimited string, return an array as output
    function delim_string_parse(str::AbstractString, delim::Char, parseType::Type; merge=false)
        parsed = Array{parseType}(ceil(Int,length(str)/2));

        # Ignore initial delimiter
        last_delim_pos=0;
        if str[1]==delim
            last_delim_pos=1;
        end

        # Cycle through string parsing text betweeen delims
        delim_pos=0;
        n = 0;
        if merge
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos+1
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = parse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                    end
                    last_delim_pos = delim_pos;
                end
            end
        else
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = parse(parseType, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        last_delim_pos = delim_pos;
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1;
            parsed[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the parsed values
        return parsed[1:n]
    end
    export delim_string_parse

    function delim_string_function(f::Function, str::AbstractString, delim::Char, outType::Type; merge=false)
        # parsed = delim_string(f::Function, str::AbstractString, delim::Char, outType::Type; merge=false)

        # Max number of delimted values
        ndelims = 2;
        for i=1:length(str)
            if str[i] == delim
                ndelims +=1
            end
        end

        # Allocate output array
        parsed = Array{outType}(ceil(Int,ndelims));

        # Ignore initial delimiter
        last_delim_pos=0;
        if str[1]==delim
            last_delim_pos=1;
        end

        # Cycle through string parsing text betweeen delims
        delim_pos=0;
        n = 0;
        if merge
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos+1
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = f(str[(last_delim_pos+1):(delim_pos-1)])
                        end
                    end
                    last_delim_pos = delim_pos;
                end
            end
        else
            for i=1:length(str)
                if str[i] == delim
                    delim_pos = i;
                    if delim_pos>last_delim_pos
                        n += 1;
                        if delim_pos>last_delim_pos+1
                            parsed[n] = f(str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        last_delim_pos = delim_pos;
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1;
            parsed[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the parsed values
        return parsed[1:n]
    end
    export delim_string_function

## --- Classifying imported datasets

    # Return true for numbers and strings that can be parsed as numbers
    function isnumeric(x)
        if isa(x,Number)
            return true
        elseif isa(x,AbstractString) && ~isnull(tryparse(Float64,x))
            return true
        else
            return false
        end
    end
    export isnumeric

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
    export nonnumeric

## --- Transforming imported datasets

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
    export floatify

    # Convert a flat array into a dict with each column as a variable
    function elementify(in::Array, elements::Array=in[1,:]; floatout::Bool=true, skipstart::Int=1, skipnameless::Bool=true)
        # Output as dictionary
        out = Dict()
        out["elements"] = elements

        # Parse the input array, minus empty-named columns
        for i=1:length(elements)
            thiscol = in[(1+skipstart):end,i]
            floatcol = floatout && ( sum(isnumeric.(thiscol)) >= sum(nonnumeric.(thiscol)) )
            includecol = ~skipnameless || (elements[i] .!= "")

            if haskey(out,elements[i]) && ( sum(isnumeric.(out[elements[i]])) >= sum(nonnumeric.(thiscol)) ) && floatcol && includecol
                # If key already exists and is numeric
                out[elements[i]] = nanmean.( hcat( floatify.(out[elements[i]]), floatify.(thiscol) ) )
            elseif floatcol && includecol
                # If key is numeric
                out[elements[i]] = floatify.(thiscol)
            elseif includecol
                # If key is non-numeric
                out[elements[i]] = thiscolumn
            end
        end
        return out
    end
    export elementify

    # Convert a dict into a flat array with variables as columns
    function unelementify(in::Dict, elements::Array=sort(collect(keys(in))); floatout::Bool=false, findnumeric::bool=false)

        # Find the elements in the input dict
        if any(elements .== "elements")
            elements = in["elements"]
        end

        # Figure out how many are numeric (if necessary)
        if findnumeric
            numericelements = Array{Bool}(length(elements))
            for i=1:length(elements)
                numericelements = sum(isnumeric.(in[elements[i]])) > sum(nonnumeric.(in[elements[i]]))
            end
            elements = elements[numericelements]
        end

        if floatout
            # Allocate output Array{Float64}
            out=Array{Float64}(length(in[elements[1]]),length(elements))

            # Parse the input dict
            for i=1:length(elements)
                out[:,i] = floatify.(in[elements[i]])
            end
        else
            # Allocate output Array{Any}
            out=Array{Any}(length(in[elements[1]])+1,length(elements))

            # Parse the input dict
            for i=1:length(elements)
                out[1,i] = elements[i]
                out[2:end,i] = in[elements[i]]
            end
        end
        return out
    end
    export unelementify

## --- End of File
