## --- Parse a delimited string

    # Parse a delimited string, return the results in a pre-allocated array provided as input
    function parseDelimString!(parsed::Array, str::String, delim::Char, parseType::Type; offset=0, merge=false)

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
    export parseDelimString!

    # Parse a delimited string, return an array as output
    function parseDelimString(str::String, delim::Char, parseType::Type; merge=false)
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
    export parseDelimString

## --- End of File
