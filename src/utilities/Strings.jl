## --- Parse a delimited string

    function parseDelimString!(parsed::Array, str::String, delim::Char, parseType::Type; offset=0)

        # Ignore initial delimiter
        last_delim_pos=0;
        if str[1]==delim
            last_delim_pos=1;
        end

        # Cycle through string parsing text betweeen delims
        delim_pos=0;
        n = offset;
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

        # Check for final value after last delim
        if length(str)>last_delim_pos
            n += 1;
            parsed[n] = parse(parseType, str[(last_delim_pos+1):length(str)])
        end

        # Return the number of parsed values
        return n-offset
    end
    export parseDelimString


## --- End of File
