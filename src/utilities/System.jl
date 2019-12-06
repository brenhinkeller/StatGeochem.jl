## --- Direct access to real libc system() function

    # Direct system() access without stripping special characters like Julia's
    # rather protective "run()". Learning new syntax for pipelining, etc. is
    # annoying. Is it really that common to be running shell commands from an
    # untrusted source?
    function system(cmdstr)
        return ccall((:system,), Int, (Cstring,), cmdstr)
    end
    export system

## --- Retain deprecated functions with matlab-like syntax,
    # to avoid breakages in user scripts that may depend on them

    if VERSION>=v"1.0"
        function linspace(l::Number,u::Number,n::Number)
            return range(l,stop=u,length=n)
        end
        export linspace

        linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
        export linreg

        function repmat(A::AbstractArray, vert::Integer)
            return repeat(A, outer=vert)
        end
        function repmat(A::AbstractArray, vert::Integer, horiz::Integer)
            return repeat(A, outer=(vert, horiz))
        end
        export repmat

        """
        ```julia
        contains(haystack::AbstractString, needle::Union{AbstractString,Regex,AbstractChar})
        ```

        Tests whether the string `haystack` contains the string or char `needle`,
        matching case. Identical to `occursin`, but with opposite argument order.
        """
        function contains(haystack::AbstractString, needle::Union{AbstractString,Regex,AbstractChar})
            return occursin(needle, haystack)
        end
        """
        ```julia
        contains(haystack, needle)
        ```

        Converts both `haystack` and `needle` to strings and checks whether
        `string(haystack)` contains `string(needle)`.
        """
        function contains(haystack, needle)
            return occursin(string(needle), string(haystack))
        end
        export contains

        """
        ```julia
        containsi(haystack::AbstractString, needle::Union{AbstractString,AbstractChar})
        ```

        Tests whether the string `haystack` contains the string or char `needle`,
        ignoring case
        """
        function containsi(haystack::AbstractString, needle::Union{AbstractString,AbstractChar})
            return occursin(lowercase(needle), lowercase(haystack))
        end
        """
        ```julia
        containsi(haystack, needle)
        ```

        Converts both `haystack` and `needle` to strings and checks whether
        `string(haystack)` contains `string(needle)`, ignoring case.
        """
        function containsi(haystack, needle)
            return occursin(lowercase(string(needle)), lowercase(string(haystack)))
        end
        export containsi
    end

## --- End of File
