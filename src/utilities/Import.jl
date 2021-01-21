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
            result[n] = f(str[(last_delim_pos+1):length(str)])
        end

        # Return the result values
        return result[1:n]
    end
    export delim_string_function

    """
    ```julia
    parsedlm(str::AbstractString, delimiter::Char, parsetype::Type=Float64; rowdelimiter::Char='\n')
    ```
    Parse a string delimited by both row and column into a single (2-D) matrix. Default column delimiter is newline.
    Similar to `readdlm`, but operating on a string instead of a file.
    """
    function parsedlm(str::AbstractString, delimiter::Char, parsetype::Type=Float64; rowdelimiter::Char='\n')
    	if parsetype <: AbstractFloat
    		emptyval = parsetype(NaN)
    	else
    		emptyval = zero(parsetype)
    	end

    	# Count rows, and find maximum number of delimiters per row
    	numcolumns = maxcolumns = maxrows = 0
    	cₗ = delimiter
    	for c in str
    		(c == delimiter) && (numcolumns += 1)
    		if c == rowdelimiter
    			maxrows += 1
    			numcolumns += 1
    			# See if we have a new maximum, and reset the counters
    			(numcolumns > maxcolumns) && (maxcolumns = numcolumns)
    			numcolumns=0
    		end
    		cₗ = c
    	end
    	# If the last line isn't blank, add one more to the row counter
    	(cₗ != rowdelimiter) && (maxrows += 1)
    	maxchars = length(str)


    	# Allocate space for the imported array and fill with emptyval
    	parsedmatrix = fill(emptyval, maxrows, maxcolumns)

    	k = kₗ = 0 # Last delimiter position
    	for i = 1:maxrows
    		for j = 1:maxcolumns
    			while (str[k+1] != delimiter) && (str[k+1] != rowdelimiter) && (k+1 < maxchars)
    				k += 1
    			end
    			# If we're at the end of the string, move on
    			(k+1 == maxchars) && break

    			# Otherwise, parse the string
    			parsed = tryparse(parsetype, str[kₗ+1:k])
    			isnothing(parsed) || (parsedmatrix[i,j] = parsed)

    			# Step over the delimiter
    			kₗ = k += 1
    			# If we've hit a row delimiter, move to next row
    			(str[k] == rowdelimiter) && break
    		end
    	end
    	return parsedmatrix
    end


## --- Classifying imported datasets

    """
    ```julia
    plausiblynumeric(x)
    ```
    Return `true` if `x` can be parsed as a number, else `false`
    """
    function plausiblynumeric(x)
        if isa(x,Number)
            return true
        elseif isa(x,AbstractString) && tryparse(Float64,x) != nothing
            return true
        else
            return false
        end
    end
    export plausiblynumeric

    """
    ```julia
    nonnumeric(x)
    ```
    Return true for if `x` is not missing but cannot be parsed as a number
    """
    function nonnumeric(x)
        if isa(x,Number)
            false
        elseif isa(x,AbstractString) && (tryparse(Float64,x) != nothing || x == "")
            false
        else
            true
        end
    end
    export nonnumeric

## --- Transforming imported datasets

    """
    ```julia
    floatify(x, T::Type=Float64)
    ```
    Convert `x` to a floating-point number (default `Float64`) by any means necessary
    """
    function floatify(x, T::Type=Float64)
        if isa(x, Number)
            T(x)
        elseif isa(x,AbstractString) && tryparse(T,x) != nothing
            parse(T, x)
        else
            T(NaN)
        end
    end

    function _columnformat(x, standardize=true, floattype=Float64)
        if standardize
            if sum(plausiblynumeric.(x)) >= sum(nonnumeric.(x))
                return floatify.(x, floattype)
            else
                return string.(x)
            end
        else
            if all(xi -> isa(xi, AbstractString), x)
                return string.(x)
            else
                return x
            end
        end
    end

    function _sanitizevarname(s::AbstractString)
        s = replace(s, r"[\[\](){}]" => "") # Remove parentheses entirely
        s = replace(s, r"^([0-9])" => s"_\1") # Can't begin with a number
        s = replace(s, r"([\0-\x1F -/:-@\[-`{-~])" => s"_") # Everything else becomes an underscore
        return s
    end


    """
    ```julia
    elementify(data::Array, elements::Array=data[1,:];
        \timportas=:Dict,
        \tstandardize::Bool=true,
        \tfloattype=Float64,
        \tskipstart::Integer=1,
        \tskipnameless::Bool=true
    )
    ```
    Convert a flat array `data` into a dictionary (`importas=:Dict`) or named
    tuple (`importas=:Tuple`) with each column as a variable.
    Tuples are substantially more efficient, so should be favored where possible.
    """
    function elementify(data::Array, elements::Array=data[1,:];
            importas=:Dict,
            skipstart::Integer=1,
            standardize::Bool=true,
            floattype=Float64,
            skipnameless::Bool=true
        )
        if importas == :Dict || importas==:dict
            # Output as dictionary
            if standardize
                # Constrain types somewhat for a modicum of type-stability
                if 1+skipstart == size(data,1)
                    result = Dict{String,Union{Vector{String}, String, Float64}}()
                else
                    result = Dict{String,Union{Vector{String}, Vector{Float64}}}()
                end
            else
                result = Dict{String, Any}()
            end

            # Process elements array
            elements = string.(elements)
            if skipnameless
                elements = elements[elements .!= ""]
            end
            result["elements"] = elements

            # Parse the input array, minus empty-named columns
            for i = 1:length(elements)
                if 1+skipstart == size(data,1)
                    column = data[end,i]
                else
                    column = data[(1+skipstart):end,i]
                end
                column_is_numeric = sum(plausiblynumeric.(column)) >= sum(nonnumeric.(column))

                if haskey(result,elements[i])
                    # If key already exists
                    if column_is_numeric && (standardize || (sum(plausiblynumeric.(result[elements[i]])) >= sum(nonnumeric.(result[elements[i]]))) )
                        # If either this column or the existing one is plausibly numeric, average the two
                        result[elements[i]] = nanmean( hcat(floatify.(result[elements[i]], floattype), floatify.(column, floattype)), dim=2 )
                    elseif standardize
                        # If neither is numeric, but standardize is set, must return a string
                        result[elements[i]] = string.(result[elements[i]]) .* "|" .* string(lastcol)
                    else
                        # If neither is plausibly numeric, contatenate the columns and move on
                        result[elements[i]] = hcat(result[elements[i]], column)
                    end
                else
                    result[elements[i]] = _columnformat(column, standardize, floattype)
                end
            end

            # Return only unique elements, since dictionary keys must be unique
            result["elements"] = unique(elements)
            return result
        elseif importas==:Tuple || importas==:tuple || importas==:NamedTuple
            # Import as NamedTuple (more efficient future default)
            t = skipnameless ? elements .!= "" : isa.(elements, AbstractString)
            elements = _sanitizevarname.(elements[t])
            symbols = ((Symbol(e) for e ∈ elements)...,)
            values = [_columnformat(data[1+skipstart:end,i], standardize, floattype) for i in findall(t)]
            return NamedTuple{symbols}(values)
        end
    end
    export elementify


    """
    ```julia
    unelementify(dataset, elements;
        \tfloatout::Bool=false,
        \tfloattype=Float64,
        \tfindnumeric::Bool=false,
        \tskipnan::Bool=false,
        \trows=:
    )
    ```
    Convert a dict or named tuple of vectors into a 2-D array with variables as columns
    """
    function unelementify(dataset::Dict, elements::Array=sort(collect(keys(dataset)));
            floatout::Bool=false,
            floattype=Float64,
            findnumeric::Bool=false,
            skipnan::Bool=false,
            rows=:
        )

        # Find the elements in the input dict if they exist and aren't otherwise specified
        if any(elements .== "elements")
            elements = dataset["elements"]
        end

        # Figure out how many are numeric (if necessary), so we can export only
        # those if `findnumeric` is set
        if findnumeric
            is_numeric_element = Array{Bool}(undef,length(elements))
            for i = 1:length(elements)
                is_numeric_element[i] = sum(plausiblynumeric.(dataset[elements[i]])) > sum(nonnumeric.(dataset[elements[i]]))
            end
            elements = elements[is_numeric_element]
        end

        # Generate output array
        if floatout
            # Allocate output Array{Float64}
            result = Array{Float64}(undef, length(dataset[elements[1]][rows]), length(elements))

            # Parse the input dict. No column names if `floatout` is set
            for i = 1:length(elements)
                result[:,i] = floatify.(dataset[elements[i]][rows], floattype)
            end
        else
            # Allocate output Array{Any}
            result = Array{Any}(undef, length(dataset[elements[1]][rows])+1, length(elements))

            # Parse the input dict
            for i = 1:length(elements)
                # Column name goes in the first row, everything else after that
                result[1,i] = elements[i]
                result[2:end,i] .= dataset[elements[i]][rows]

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
    function unelementify(dataset::NamedTuple, elements=keys(dataset);
            floatout::Bool=false,
            floattype=Float64,
            findnumeric::Bool=false,
            skipnan::Bool=false,
            rows=:
        )
        # Figure out how many are numeric (if necessary), so we can export only
        # those if `findnumeric` is set
        if findnumeric
            elements = filter(x -> sum(plausiblynumeric.(dataset[x])) > sum(nonnumeric.(dataset[x])), elements)
        end

        # Generate output array
        if floatout
            # Allocate output Array{Float64}
            result = Array{floattype}(undef,length(dataset[elements[1]][rows]),length(elements))

            # Parse the input dict. No column names if `floatout` is set
            for i = 1:length(elements)
                result[:,i] = floatify.(dataset[elements[i]][rows], floattype)
            end
        else
            # Allocate output Array{Any}
            result = Array{Any}(undef,length(dataset[elements[1]][rows])+1,length(elements))

            # Parse the input dict
            for i = 1:length(elements)
                # Column name goes in the first row, everything else after that
                result[1,i] = string(elements[i])
                result[2:end,i] .= dataset[elements[i]][rows]

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

## --- Concatenating / stacking datasets

    """
    ```julia
    concatenatedatasets(d1::AbstractDict, d2::AbstractDict)
    ```
    Vertically concatenate two Dict-based datasets, variable-by-variable
    """
    function concatenatedatasets(d1::AbstractDict, d2::AbstractDict)
        if isempty(keys(d1))
            result = d2
        elseif isempty(keys(d2))
            result = d1
        else
            # If 'elements' field doesn't exist, populate it
            if ~haskey(d1,"elements")
                d1["elements"] = sort(collect(keys(d1)))
            end
            if ~haskey(d2,"elements")
                d2["elements"] = sort(collect(keys(d2)))
            end

            # Find variable size
            s1 = size(d1[d1["elements"][1]])
            s2 = size(d2[d2["elements"][1]])

            # Combine datasets
            result = Dict()
            result["elements"] = d1["elements"] ∪ d2["elements"]
            for e in result["elements"]
                # Make any missing fields
                if haskey(d1, e) && ~haskey(d2, e)
                    if eltype(d1[e]) <: Number
                        d2[e] = fill(float(eltype(d1[e]))(NaN), s2)
                    else
                        d2[e] = fill("", s2)
                    end
                elseif haskey(d2, e) && ~haskey(d1, e)
                    if eltype(d2[e]) <: Number
                        d1[e] = fill(float(eltype(d2[e]))(NaN), s1)
                    else
                        d1[e] = fill("", s1)
                    end
                end

                # Combine fields
                result[e] = vcat(d1[e], d2[e])
            end
        end
        return result
    end
    export concatenatedatasets


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
    function renormalize!(A::AbstractArray{<:Number}; dim::Number=0, total::Number=1)
        current_sum = nansum(A, dim=dim)
        A .*= total ./ current_sum
    end
    export renormalize!

## --- High-level import/export functions

    """
    ```julia
    function importdataset(filepath, delim;
        \timportas=:Dict,
        \tstandardize::Bool=true,
        \tfloattype=Float64,
        \tskipstart::Integer=0,
        \tskipnameless::Bool=true,
        \tmindefinedcolumns::Integer=0
    )
    ```
    Import a delimited file specified by `filepath` with delimiter `delim` as a
    dataset in the form of either a `Dict` or a `NamedTuple`.

    Possible keyword arguments include:

        \importas
    Specify the format of the imported dataset. Options include `:Dict` and `:Tuple`

        \tstandardize
    Convert columns to uniform type wherever possible. Boolean; `true` by default.

        \tfloattype
    Preferred floating-point type for numerical data. `Float64` by default

        \tskipstart
    Ignore this many rows at the start of the input file (useful if input file has
    a header or other text before the column names). `0` by default

        \tskipnameless
    Skip columns with no column name. Boolean; `true` by default

        \tmindefinedcolumns
    Skip rows with fewer than this number of delimiters. `0` by default.
    """
    function importdataset(filepath::AbstractString, delim::AbstractChar;
            importas=:Dict,
            standardize::Bool=true,
            floattype=Float64,
            skipstart::Integer=0,
            skipnameless::Bool=true,
            mindefinedcolumns::Integer=0
        )

        # Read file
        data = readdlm(filepath, delim, skipstart=skipstart)
        # Exclude rows with fewer than `mindefinedcolumns` columns
        if mindefinedcolumns > 0
            definedcolumns = vec(sum(.~ isempty.(data), dims=2))
            t = definedcolumns .>= mindefinedcolumns
            data = data[t,:]
        end

        return elementify(data,
            importas=importas,
            standardize=standardize,
            floattype=floattype,
            skipstart=1,
            skipnameless=skipnameless
        )
    end
    export importdataset

    """
    ```julia
    exportdataset(dataset, [elements], filepath, delim;
        \tfloatout::Bool=false,
        \tfindnumeric::Bool=false,
        \tskipnan::Bool=true,
        \tdigits::Integer,
        \tsigdigits::Integer
        \trows=:
    )
    ```
    Convert a dict or named tuple of vectors into a 2-D array with variables as columns
    Export a `dataset` (in the form of either a `Dict` or a `NamedTuple`),
    optionally specifying which `elements` to export, as a delimited ASCII text file
    with the name specified by `filepath` and delimiter `delim`.

    Possible keyword arguments include:

        \tdigits
        \tsigdigits
    Specify a number of absolute or significant digits to which to round the printed output

        \tskipnan
    Leave `NaN`s as empty cells in the delimited output file

        \tfloatout
    Force all output to be represented as a floating-point number, or else `NaN`

        \tfindnumeric
    Export only numeric columns

        \trows
    specify which rows of the dataset to export (default `:` exports all rows)
    """
    function exportdataset(dataset::Union{Dict,NamedTuple}, filepath::AbstractString, delim::AbstractChar;
            floatout::Bool=false,
            findnumeric::Bool=false,
            skipnan::Bool=true,
            digits::Integer=0,
            sigdigits::Integer=0,
            rows=:
        )

        # Convert dataset to flat 2d array
        data = unelementify(dataset,
            floatout=floatout,
            findnumeric=findnumeric,
            skipnan=skipnan,
            rows=rows
        )

        # Round output if applicable
        if sigdigits > 0
            data .= data .|> x -> isa(x, Number) ? round(x, sigdigits=sigdigits) : x
        elseif digits > 0
            data .= data .|> x -> isa(x, Number) ? round(x, digits=digits) : x
        end

        # Write to file
        return writedlm(filepath, data, delim)
    end
    # As above, but specifying which elements to export
    function exportdataset(dataset::Union{Dict,NamedTuple}, elements::Array, filepath::AbstractString, delim::AbstractChar;
            floatout::Bool=false,
            findnumeric::Bool=false,
            skipnan::Bool=true,
            digits::Integer=0,
            sigdigits::Integer=0,
            rows=:
        )

        # Convert dataset to flat 2d array
        data = unelementify(dataset, elements,
            floatout=floatout,
            findnumeric=findnumeric,
            skipnan=skipnan,
            rows=rows
        )

        # Round output if applicable
        if sigdigits > 0
            data .= data .|> x -> isa(x, Number) ? round(x, sigdigits=sigdigits) : x
        elseif digits > 0
            data .= data .|> x -> isa(x, Number) ? round(x, digits=digits) : x
        end

        # Write to file
        return writedlm(filepath, data, delim)
    end
    export exportdataset

## --- End of File
