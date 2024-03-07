## --- Parse a delimited string


    """
    ```julia
    delim_string_parse!(result, str, delim, [T];
        \toffset::Integer=0,
        \tmerge::Bool=false,
        \tundefval=NaN)
    ```
    Parse a delimited string `str` with delimiter `delim` into values of type `T`
    and return the answers in a pre-allocated `result` array provided as input.

    If `T` is not specified explicitly, the `eltype` of the `result` array will
    be used by default.

    Optional keyword arguments and defaults:

        offset::Integer=0

    Start writing the parsed results into `result` at index `1+offset`

        merge::Bool=false

    Merge repeated delimiters?

        undefval=NaN

    A value to subsitute for any value that cannot be `parse`d to type `T`.

    See also `delim_string_parse` for a non-in-place version that will automatically
    allocate a result array.

    ### Examples
    ```julia
    julia> A = zeros(100);

    julia> n = delim_string_parse!(A, "1,2,3,4,5", ',', Float64)
    5

    julia> A[1:n]
    5-element Vector{Float64}:
     1.0
     2.0
     3.0
     4.0
     5.0
    ```
    """
    function delim_string_parse!(result::Array, str::AbstractString, delim::Char, T::Type=eltype(result); offset::Integer=0, merge::Bool=false, undefval=NaN)

        # Ignore initial delimiter
        last_delim_pos = 0
        if ~isempty(str) && first(str) == delim
            last_delim_pos = 1
        end

        # Cycle through string parsing text betweeen delims
        delim_pos = 0
        n = offset
        if merge
            for i ∈ eachindex(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos+1
                        n += 1
                        parsed = nothing
                        if delim_pos > last_delim_pos+1
                            parsed = tryparse(T, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? T(undefval) : parsed
                    end
                    last_delim_pos = delim_pos
                end
            end
        else
            for i ∈ eachindex(str)
                if str[i] == delim
                    delim_pos = i
                    if delim_pos > last_delim_pos
                        n += 1
                        parsed = nothing
                        if delim_pos > last_delim_pos+1
                            parsed = tryparse(T, str[(last_delim_pos+1):(delim_pos-1)])
                        end
                        result[n] = isnothing(parsed) ? T(undefval) : parsed
                        last_delim_pos = delim_pos
                    end
                end
            end
        end

        # Check for final value after last delim
        if length(str) > last_delim_pos
            n += 1
            parsed = tryparse(T, str[(last_delim_pos+1):length(str)])
            result[n] = isnothing(parsed) ? T(undefval) : parsed
        end

        # Return the number of result values
        return n-offset
    end
    export delim_string_parse!

    """
    ```julia
    delim_string_parse(str, delim, T;
        \tmerge::Bool=false,
        \tundefval=NaN)
    ```
    Parse a delimited string `str` with delimiter `delim` into values of type `T`
    and return the answers as an array with eltype `T`

    Optional keyword arguments and defaults:

        merge::Bool=false

    Merge repeated delimiters?

        undefval=NaN

    A value to subsitute for any value that cannot be `parse`d to type `T`.

    See also `delim_string_parse!` for an in-place version.

    ### Examples
    ```julia
    julia> delim_string_parse("1,2,3,4,5", ',', Float64)
    5-element Vector{Float64}:
     1.0
     2.0
     3.0
     4.0
     5.0
    ```
    """
    function delim_string_parse(str::AbstractString, delim::Char, T::Type=Float64; merge::Bool=false, undefval=NaN)

        # Allocate an array to hold our parsed results
        result = Array{T}(undef,ceil(Int,length(str)/2))

        # Parse the string
        n = delim_string_parse!(result, str, delim, T; merge=merge, undefval=undefval)

        # Return the result values
        return result[1:n]
    end
    export delim_string_parse

    """
    ```julia
    delim_string_function(f, str, delim, T;
        \tmerge::Bool=false,
    ```
    Parse a delimited string `str` with delimiter `delim` into substrings that will
    then be operated upon by function `f`. The results of `f` will be returned
    in an array with eltype `T`.

    ### Examples
    ```julia
    julia> delim_string_function(x -> delim_string_parse(x, ',', Int32, undefval=0), "1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", '\n', Array{Int32,1})
    4-element Vector{Vector{Int32}}:
     [1, 2, 3, 4]
     [5, 6, 7, 8]
     [9, 10, 11, 12]
     [13, 14, 15, 16]
    ```
    """
    function delim_string_function(f::Function, str::AbstractString, delim::Char, T::Type; merge::Bool=false)

        # Max number of delimted values
        ndelims = 2
        for i ∈ eachindex(str)
            if str[i] == delim
                ndelims += 1
            end
        end

        # Allocate output array
        result = Array{T}(undef,ceil(Int,ndelims))

        # Ignore initial delimiter
        last_delim_pos = 0
        if first(str) == delim
            last_delim_pos = 1
        end

        # Cycle through string parsing text betweeen delims
        delim_pos = 0
        n = 0
        if merge
            for i ∈ eachindex(str)
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
            for i ∈ eachindex(str)
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
    parsedlm(str::AbstractString, delimiter::Char, T::Type=Float64; rowdelimiter::Char='\\n')
    ```
    Parse a string delimited by both row and column into a single (2-D) matrix. Default column delimiter is newline.
    Similar to `readdlm`, but operating on a string instead of a file.

    ### Examples
    ```julia
    julia> parsedlm("1,2,3\n4,5,6\n7,8,9\n", ',', Float64)
    3×3 Matrix{Float64}:
     1.0  2.0  3.0
     4.0  5.0  6.0
     7.0  8.0  9.0

    julia> parsedlm("1,2,3,4\n5,6,7,8\n9,10,11,12\n13,14,15,16", ',', Int64)
    4×4 Matrix{Int64}:
      1   2   3   4
      5   6   7   8
      9  10  11  12
     13  14  15  16
    ```
    """
    function parsedlm(str::AbstractString, delimiter::Char, ::Type{T}=Float64; rowdelimiter::Char='\n') where {T}

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

    	# Allocate space for the imported array and fill with emptyval
    	parsedmatrix = emptys(T, maxrows, maxcolumns)

        maxchars = length(str)
    	kₗ = kₙ = firstindex(str) # Last and next delimiter position
    	@inbounds for i = 1:maxrows
    		for j = 1:maxcolumns
                c = str[kₙ]
    			while (kₙ < maxchars) && (c !== delimiter) && (c !== rowdelimiter)
    				kₙ = nextind(str, kₙ)
                    c = str[kₙ]
    			end

                if kₙ>kₗ
                    # Parse the string
                    k = (c===delimiter || c===rowdelimiter) ? prevind(str,kₙ) : kₙ
        			parsed = tryparse(T, str[kₗ:k])
        			isnothing(parsed) || (parsedmatrix[i,j] = parsed)
                end

                # If we're at the end of the string, move on
                (kₙ == maxchars) && break

    			# Step over the delimiter
    			kₗ = kₙ = nextind(str, kₙ)

                # If we've hit a row delimiter, move to next row
                (str[kₙ] == rowdelimiter) && break
    		end
    	end
    	return parsedmatrix
    end
    export parsedlm


## --- Classifying imported datasets

    """
    ```julia
    isnumeric(x)
    ```
    Return `true` if `x` can be parsed as a number, else `false`

    ### Examples
    ```julia
    julia> StatGeochem.isnumeric(1)
    true

    julia> StatGeochem.isnumeric("1")
    true

    julia> StatGeochem.isnumeric("0.5e9")
    true

    julia> StatGeochem.isnumeric("foo")
    false
    ```
    """
    isnumeric(x) = false
    isnumeric(x::Number) = true
    isnumeric(x::AbstractString) = tryparse(Float64,x) !== nothing

    """
    ```julia
    nonnumeric(x)
    ```
    Return true for if `x` is not missing but cannot be parsed as a number

    ### Examples
    ```julia
    julia> StatGeochem.nonnumeric(1)
    false

    julia> StatGeochem.nonnumeric("1")
    false

    julia> StatGeochem.nonnumeric("0.5e9")
    false

    julia> StatGeochem.nonnumeric("foo")
    true
    ```
    """
    nonnumeric(x) = true
    nonnumeric(x::Number) = false
    nonnumeric(x::Missing) = false
    nonnumeric(x::AbstractString) = (tryparse(Float64,x) === nothing) && (x != "")


## --- Transforming imported datasets

    """
    ```julia
    floatify(x, T::Type=Float64)
    ```
    Convert `x` to a floating-point number (default `Float64`) by any means necessary

    ### Examples
    ```julia
    julia> StatGeochem.floatify(5)
    5.0

    julia> StatGeochem.floatify("5")
    5.0

    julia> StatGeochem.floatify("0x05")
    5.0

    julia> StatGeochem.floatify("0.5e1")
    5.0
    ```
    """
    floatify(x, T::Type{<:AbstractFloat}=Float64) = T(NaN)
    floatify(x::Number, T::Type{<:AbstractFloat}=Float64) = T(x)
    floatify(x::AbstractString, T::Type{<:AbstractFloat}=Float64) = (n = tryparse(T,x)) !== nothing ? n : T(NaN)


    columnformat(x, standardize::Bool=true, floattype=Float64) = _columnformat(x, Val(standardize), floattype)
    function _columnformat(x, ::Val{true}, floattype)
        if sum(isnumeric.(x)) >= sum(nonnumeric.(x))
            floatify.(x, floattype)
        else
            string.(x)
        end
    end
    function _columnformat(x, ::Val{false}, floattype)
        if all(xi -> isa(xi, AbstractString), x)
            string.(x)
        elseif all(xi -> isa(xi, AbstractFloat), x)
            float.(x)
        elseif all(xi -> isa(xi, Integer), x)
            Integer.(x)
        else
            x
        end
    end


    """
    ```julia
    sanitizevarname(s::AbstractString)
    ```
    Modify an input string `s` to transform it into an acceptable variable name.

    ### Examples
    ```julia
    julia> StatGeochem.sanitizevarname("foo")
    "foo"

    julia> StatGeochem.sanitizevarname("523foo")
    "_523foo"

    julia> StatGeochem.sanitizevarname("Length (μm)")
    "Length_μm"
    ```
    """
    function sanitizevarname(s::AbstractString)
        s = replace(s, r"[\[\](){}]" => "") # Remove parentheses entirely
        s = replace(s, r"^([0-9])" => s"_\1") # Can't begin with a number
        s = replace(s, r"([\0-\x1F -/:-@\[-`{-~])" => s"_") # Everything else becomes an underscore
        return s
    end
    sanitizevarname(s::Symbol) = s

    symboltuple(x::NTuple{N, Symbol}) where {N} = x
    symboltuple(x::NTuple{N}) where {N} = ntuple(i->Symbol(x[i]), N)
    symboltuple(x) = ((Symbol(s) for s in x)...,)

    stringarray(x::Vector{String}) = x
    stringarray(x::NTuple{N, String}) where {N} = [s for s in x]
    stringarray(x) = [String(s) for s in x]


    """
    ```julia
    TupleDataset(d::Dict, elements=keys(d))
    ```
    Convert a dict-based dataset to a tuple-based dataset. 

    See also `DictDataset`

    ### Examples
    ```julia
    julia> d
    Dict{String, Vector{Float64}} with 2 entries:
      "Yb" => [0.823733, 0.0531003, 0.47996, 0.560998, 0.001816, 0.455064, 0.694017, 0.737816, 0.0755015, 0.46098 …
      "La" => [0.440947, 0.937551, 0.464318, 0.694184, 0.253974, 0.521292, 0.857979, 0.0545946, 0.716639, 0.597616…
    
    julia> t = TupleDataset(d)
    NamedTuple with 2 elements:
      Yb  = Vector{Float64}(100,)   [0.8237334494155881 ... 0.012863893327602627]
      La  = Vector{Float64}(100,)   [0.44094669199955616 ... 0.5371416189174069]
    ```
    """
    function TupleDataset(d::Dict, elements=haskey(d,"elements") ? d["elements"] : keys(d))
        symbols = symboltuple(sanitizevarname.(elements))
        return NamedTuple{symbols}(d[e] for e in elements)
    end
    export TupleDataset


    """
    ```julia
    DictDataset(t::NamedTuple, elements=keys(t))
    ```
    Convert a tuple-based dataset to a dict-based dataset. 

    See also `TupleDataset`

    ### Examples
    ```julia
    julia> t 
    NamedTuple with 2 elements:
      La  = Vector{Float64}(100,)   [0.6809734028326375 ... 0.30665937715972313]
      Yb  = Vector{Float64}(100,)   [0.8851029525168138 ... 0.866246147690925]
    
    julia> d = DictDataset(t)
    Dict{String, Vector{Float64}} with 2 entries:
      "Yb" => [0.885103, 0.284384, 0.351527, 0.643542, 0.631274, 0.653966, 0.968414, 0.00204819, 0.0655173, 0.5343…
      "La" => [0.680973, 0.35098, 0.0198742, 0.139642, 0.0703337, 0.0328973, 0.639431, 0.245205, 0.424142, 0.48889…    
    ```
    """
    function DictDataset(t::NamedTuple, elements=keys(t))
        d = Dict(String(e) => t[Symbol(e)] for e in elements)
    end
    export DictDataset

    """
    ```julia
    elementify(data::AbstractArray, [elements=data[1,:]];
        \timportas=:Dict,
        \tstandardize::Bool=true,
        \tfloattype=Float64,
        \tskipstart::Integer=1,
        \tskipnameless::Bool=true
    )
    ```
    Convert a flat array `data` into a Named Tuple (`importas=:Tuple`) or
    Dictionary (`importas=:Dict`) with each column as a variable.
    Tuples are substantially more efficient, so should be favored where possible.

    ### Examples
    ```julia
    julia> A = ["La" "Ce" "Pr"; 1.5 1.1 1.0; 3.7 2.9 2.5]
    3×3 Matrix{Any}:
      "La"   "Ce"   "Pr"
     1.5    1.1    1.0
     3.7    2.9    2.5

    julia> elementify(A, importas=:Tuple)
    NamedTuple with 3 elements:
    La  = Vector{Float64}(2,) [1.5 ... 3.7]
    Ce  = Vector{Float64}(2,) [1.1 ... 2.9]
    Pr  = Vector{Float64}(2,) [1.0 ... 2.5]

    julia> elementify(A, importas=:Dict)
    Dict{String, Union{Vector{Float64}, Vector{String}}} with 4 entries:
      "Ce"       => [1.1, 2.9]
      "Pr"       => [1.0, 2.5]
      "elements" => ["La", "Ce", "Pr"]
      "La"       => [1.5, 3.7]
    ```
    """
    function elementify(data::AbstractArray;
            importas=:Tuple,
            skipstart::Integer=1,
            standardize::Bool=true,
            floattype=Float64,
            skipnameless::Bool=true,
            sumduplicates::Bool=false
        )
        elementify(data, data[firstindex(data),:];
            importas=importas,
            skipstart=skipstart,
            standardize=standardize,
            floattype=floattype,
            skipnameless=skipnameless,
            sumduplicates=sumduplicates)
    end
    function elementify(data::AbstractArray, elements;
            importas=:Tuple,
            skipstart::Integer=0,
            standardize::Bool=true,
            floattype=Float64,
            skipnameless::Bool=true,
            sumduplicates::Bool=false
        )
        if importas === :Dict || importas === :dict
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
            elements = stringarray(elements)
            if skipnameless
                elements = filter(!isempty, elements)
            end
            result["elements"] = isa(elements, Vector) ? elements : collect(elements)

            # Parse the input array, minus empty-named columns
            i₀ = firstindex(data) + skipstart
            for j ∈ eachindex(elements)
                if skipstart == size(data,1)-1
                    column = data[end,j]
                else
                    column = data[i₀:end,j]
                end

                if !haskey(result, elements[j])
                    result[elements[j]] = columnformat(column, standardize, floattype)
                else
                    lastcol = result[elements[j]]
                    treat_as_numbers = ((sum(isnumeric.(column)) >= sum(nonnumeric.(column))) || (sum(isnumeric.(lastcol)) >= sum(nonnumeric.(lastcol))))
                    if treat_as_numbers 
                        if sumduplicates
                            @info "Duplicate key $(elements[j]) found, summing"
                            result[elements[j]] = nanadd(floatify.(lastcol, floattype), floatify.(column, floattype))
                        else
                            @info "Duplicate key $(elements[j]) found, averaging"
                            result[elements[j]] = nanadd(floatify.(lastcol, floattype), floatify.(column, floattype)) ./ 2.0
                        end
                    else
                        n = 1
                        while haskey(result, elements[j]*string(n))
                            n+=1
                        end
                        @info "Duplicate key $(elements[j]) found, replaced with $(elements[j]*string(n))"
                        elements[j] = elements[j]*string(n)
                        result[elements[j]] = columnformat(column, standardize, floattype)
                    end
                end
            end

            # Return only unique elements, since dictionary keys must be unique
            result["elements"] = unique(elements)
            return result
        elseif importas==:Tuple || importas==:tuple || importas==:NamedTuple
            # Import as NamedTuple (more efficient future default)
            t = Bool[(skipnameless && e !== "") for e in elements]
            elements = sanitizevarname.(elements[t])
            i₀ = firstindex(data) + skipstart
            values = (columnformat(data[i₀:end, j], standardize, floattype) for j in findall(vec(t)))
            return NamedTuple{symboltuple(elements)}(values)
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
    Convert a Dict or Named Tuple of vectors into a 2-D array with variables as columns

    ### Examples
    ```julia
    julia> D
    NamedTuple with 3 elements:
      La  = Vector{Float64}(2,) [1.5 ... 3.7]
      Ce  = Vector{Float64}(2,) [1.1 ... 2.9]
      Pr  = Vector{Float64}(2,) [1.0 ... 2.5]

    julia> unelementify(D)
    3×3 Matrix{Any}:
      "La"   "Ce"   "Pr"
     1.5    1.1    1.0
     3.7    2.9    2.5
    ```
    """
    function unelementify(dataset::Dict, elements=sort(collect(keys(dataset)));
            floatout::Bool=false,
            floattype=Float64,
            findnumeric::Bool=false,
            skipnan::Bool=false,
            rows=:
        )

        # Find the elements in the input dict if they exist and aren't otherwise specified
        if any(elements .== "elements")
            elements = stringarray(dataset["elements"])
        end

        # Figure out how many are numeric (if necessary), so we can export only
        # those if `findnumeric` is set
        if findnumeric
            is_numeric_element = Array{Bool}(undef,length(elements))
            for i ∈ eachindex(elements)
                is_numeric_element[i] = sum(isnumeric.(dataset[elements[i]])) > sum(nonnumeric.(dataset[elements[i]]))
            end
            elements = elements[is_numeric_element]
        end

        # Generate output array
        if floatout
            # Allocate output Array{Float64}
            result = Array{Float64}(undef, length(dataset[first(elements)][rows]), length(elements))

            # Parse the input dict. No column names if `floatout` is set
            for i ∈ eachindex(elements)
                result[:,i] = floatify.(dataset[elements[i]][rows], floattype)
            end
        else
            # Allocate output Array{Any}
            result = Array{Any}(undef, length(dataset[first(elements)][rows])+1, length(elements))

            # Parse the input dict
            for i ∈ eachindex(elements)
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
        elements = symboltuple(elements)
        if findnumeric
            elements = filter(x -> sum(isnumeric.(dataset[x])) > sum(nonnumeric.(dataset[x])), elements)
        end

        # Generate output array
        if floatout
            # Allocate output Array{Float64}
            result = Array{floattype}(undef,length(dataset[first(elements)][rows]),length(elements))

            # Parse the input dict. No column names if `floatout` is set
            for i ∈ eachindex(elements)
                result[:,i] = floatify.(dataset[elements[i]][rows], floattype)
            end
        else
            # Allocate output Array{Any}
            result = Array{Any}(undef,length(dataset[first(elements)][rows])+1,length(elements))

            # Parse the input dict
            for i ∈ eachindex(elements)
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

    # Fill an array with the designated empty type
    emptys(::Type, s...) = fill(missing, s...)
    emptys(::Type{T}, s...) where T <: AbstractString = fill("", s...)
    emptys(::Type{T}, s...) where T <: Number = fill(NaN, s...)
    emptys(::Type{T}, s...) where T <: AbstractFloat = fill(T(NaN), s...)

    """
    ```julia
    concatenatedatasets(d1::NamedTuple, d2::NamedTuple,... ;[elements::Vector{Symbol}])
    concatenatedatasets(d1::AbstractDict, d2::AbstractDict,... ;[elements::Vector{String}])
    ```
    Vertically concatenate two or more Dict- or Tuple-based datasets, variable-by-variable.
    Optionally, a list of variables to include may be specified in `elements`

    ### Examples
    ```julia
    julia> d1 = Dict("La" => rand(5), "Yb" => rand(5))
    Dict{String, Vector{Float64}} with 2 entries:
      "Yb" => [0.221085, 0.203369, 0.0657271, 0.124606, 0.0975556]
      "La" => [0.298578, 0.481674, 0.888624, 0.632234, 0.564491]

    julia> d2 = Dict("La" => rand(5), "Ce" => rand(5))
    Dict{String, Vector{Float64}} with 2 entries:
      "Ce" => [0.0979752, 0.108585, 0.718315, 0.771128, 0.698499]
      "La" => [0.538215, 0.633298, 0.981322, 0.908532, 0.77754]

    julia> concatenatedatasets(d1,d2)
    Dict{String, Vector{Float64}} with 3 entries:
      "Ce" => [NaN, NaN, NaN, NaN, NaN, 0.0979752, 0.108585, 0.718315, 0.771128, 0.698499]
      "Yb" => [0.221085, 0.203369, 0.0657271, 0.124606, 0.0975556, NaN, NaN, NaN, NaN, NaN]
      "La" => [0.298578, 0.481674, 0.888624, 0.632234, 0.564491, 0.538215, 0.633298, 0.981322, 0.908532, 0.77754]
    ```
    """
    concatenatedatasets(args...; kwargs...) = concatenatedatasets((args...,); kwargs...)
    function concatenatedatasets(dst::Tuple; kwargs...)
        if length(dst) == 1
            only(dst)
        elseif length(dst) == 2
            concatenatedatasets(dst[1], dst[2]; kwargs...)
        else
            c = concatenatedatasets(dst[1], dst[2]; kwargs...)
            concatenatedatasets((c, dst[3:end]...); kwargs...)
        end
    end
    function concatenatedatasets(d1::AbstractDict, d2::AbstractDict; elements=String[])
        # Return early if either is empty
        isempty(d1) && return d2
        isempty(d2) && return d1

        # Determine keys to include. Use "elements" field if it exists
        d1ₑ = haskey(d1,"elements") ? d1["elements"] : sort(collect(keys(d1)))
        d2ₑ = haskey(d2,"elements") ? d2["elements"] : sort(collect(keys(d2)))
        available = d1ₑ ∪ d2ₑ
        if isempty(elements)
            elementsᵢ = available
        else
            elementsᵢ = elements ∩ available
        end

        # Combine datasets
        s1, s2 = size(d1[first(d1ₑ)]),  size(d2[first(d2ₑ)])
        result = typeof(d1)(e => vcombine(d1,d2,e,s1,s2) for e in elementsᵢ)
        haskey(d1,"elements") && (result["elements"] = elementsᵢ)
        return result
    end
    function concatenatedatasets(d1::NamedTuple, d2::NamedTuple; elements=Symbol[])
        # Return early if either is empty
        isempty(d1) && return d2
        isempty(d2) && return d1

        # Determine keys to include
        available = keys(d1) ∪ keys(d2)
        if isempty(elements)
            elementsᵢ = available
        else
            elementsᵢ = elements ∩ available
        end

        # Combine datasets
        s1, s2 = size(d1[first(keys(d1))]), size(d2[first(keys(d2))])
        return NamedTuple{(elementsᵢ...,)}(vcombine(d1,d2,e,s1,s2) for e in elementsᵢ)
    end
    # Vertically concatenate the fields `e` (if present) of two named tuples
    function vcombine(d1, d2, e, s1=size(d1[first(keys(d1))]), s2=size(d2[first(keys(d2))]))
        if haskey(d1,e) && ~haskey(d2,e)
            T = eltype(d1[e])
            vcat(d1[e], emptys(T, s2))
        elseif ~haskey(d1,e) && haskey(d2,e)
            T = eltype(d2[e])
            vcat(emptys(T, s1), d2[e])
        else
            vcat(d1[e], d2[e])
        end
    end
    export concatenatedatasets


## --- Renormalization of imported datasets

    """
    ```julia
    renormalize!(A::AbstractArray; dim, total=1.0)
    ```
    Normalize an array `A` in place such that it sums to `total`. Optionally may
    specify a dimension `dim` along which to normalize.
    """
    function renormalize!(A::AbstractArray; dim=:, total=1.0)
        current_sum = NaNStatistics._nansum(A, dim)
        A .*= total ./ current_sum
    end
    """
    ```julia
    renormalize!(dataset, [elements]; total=1.0)
    ```
    Normalize in-place a (i.e., compositional) `dataset` defined by a `Dict` or
    `NamedTuple` of one-dimensional numerical arrays, such that all the `elements`
    (i.e., variables -- by default all keys in the datset) sum to a given `total`
    (by default, `1.0`).

    Note that the arrays representing each element or variable are assumed to be
    of uniform length
    """
    function renormalize!(dataset::Union{Dict,NamedTuple}, elements=keys(dataset); total=1.0)
        # Note that this assumes all variables in the dataset are the same length!
        current_sum = zeros(size(dataset[first(keys(dataset))]))
        for e in elements
            current_sum .+= dataset[e] .|> x -> isnan(x) ? 0 : x
        end
        current_sum[current_sum .== 0] .= NaN

        for e in elements
            dataset[e] .*= total ./ current_sum
        end
        return dataset
    end
    export renormalize!


## --- High-level import/export functions

    function guessdelimiter(s::AbstractString)
        if length(s)>3 
            if s[end-3:end] == ".csv"
                ','
            elseif s[end-3:end] == ".tsv"
                '\t'
            elseif s[end-3:end] == ".psv"
                '|'
            else
                '\t'
            end
        else
            '\t'
        end
    end

    """
    ```julia
    function importdataset(filepath, [delim];
        \timportas=:Dict,
        \telements=nothing,
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

        \timportas
    Specify the format of the imported dataset. Options include `:Dict` and `:Tuple`

        \telements
    Specify the names to be used for each element (i.e., column) of the dataset.
    Default value (`nothing`) will cause `elements` to be read from the first row of the file

        \tstandardize
    Convert columns to uniform type wherever possible. Boolean; `true` by default.

        \tfloattype
    Preferred floating-point type for numerical data. `Float64` by default.

        \tskipstart
    Ignore this many rows at the start of the input file (useful if input file has
    a header or other text before the column names). `0` by default.

        \tskipnameless
    Skip columns with no column name. Boolean; `true` by default

        \tmindefinedcolumns
    Skip rows with fewer than this number of delimiters. `0` by default.
    """
    function importdataset(filepath::AbstractString, delim::AbstractChar=guessdelimiter(filepath);
            importas=:Dict,
            elements=nothing,
            standardize::Bool=true,
            floattype=Float64,
            skipstart::Integer=0,
            skipnameless::Bool=true,
            mindefinedcolumns::Integer=0
        )

        # Read file
        io = open(filepath, "r")
        if read(io, Char) == '\ufeff'
            @warn """Skipping hidden \'\\ufeff\' (U+FEFF) character at start of input file.

            This character is often added to CSV files by Microsoft Excel (and some other
            Microsoft products) as what appears to be what we might call an "extension",
            which would would cause file parsing to fail if we didn't manually remove it.

            Try using open software like LibreOffice instead of Excel to make this warning go away.
            """
        else
            seekstart(io)
        end
        data = readdlm(io, delim, skipstart=skipstart)
        close(io)
        # Exclude rows with fewer than `mindefinedcolumns` columns
        if mindefinedcolumns > 0
            definedcolumns = vec(sum(.~ isempty.(data), dims=2))
            t = definedcolumns .>= mindefinedcolumns
            data = data[t,:]
        end

        if isnothing(elements)
            return elementify(data,
                importas=importas,
                standardize=standardize,
                floattype=floattype,
                skipnameless=skipnameless
            )
        else
            return elementify(data, elements,
                importas=importas,
                standardize=standardize,
                floattype=floattype,
                skipnameless=skipnameless
            )
        end
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
    Specify a number of absolute or significant digits to which to round the printed output.
    Default is no rounding.

        \tskipnan
    Leave `NaN`s as empty cells in the delimited output file. Boolean; `true` by default.

        \tfloatout
    Force all output to be represented as a floating-point number, or else `NaN`.
    Boolean; `false` by default.

        \tfindnumeric
    Export only numeric columns. Boolean; `false` by default.

        \trows
    specify which rows of the dataset to export. Default `:` exports all rows.
    """
    function exportdataset(dataset::Union{Dict,NamedTuple}, filepath::AbstractString, delim::AbstractChar=guessdelimiter(filepath);
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
        if digits > 0
            map!(x -> isa(x, Number) ? round(x, digits=digits) : x, data, data)
        end
        if sigdigits > 0
            map!(x -> isa(x, Number) ? round(x, sigdigits=sigdigits) : x, data, data)
        end

        # Write to file
        return writedlm(filepath, data, delim)
    end
    # As above, but specifying which elements to export
    function exportdataset(dataset::Union{Dict,NamedTuple}, elements::Array, filepath::AbstractString, delim::AbstractChar=guessdelimiter(filepath);
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
        if digits > 0
            map!(x -> isa(x, Number) ? round(x, digits=digits) : x, data, data)
        end
        if sigdigits > 0
            map!(x -> isa(x, Number) ? round(x, sigdigits=sigdigits) : x, data, data)
        end

        # Write to file
        return writedlm(filepath, data, delim)
    end
    export exportdataset

## --- End of File
