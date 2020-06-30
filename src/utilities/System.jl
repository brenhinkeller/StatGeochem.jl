## --- Direct access to real libc system() function

    """
    ```julia
    system(cmdstr::AbstractString)
    ```

    Direct access to the command line through C's `system` function -- without
    stripping/sanitizing special characters,  in contrast to Julia's safer
    `run()` function. This allows pipelining, etc. in shell commands. Returns
    0 on success.
    """
    function system(cmdstr::AbstractString)
        return ccall((:system,), Int, (Cstring,), cmdstr)
    end
    export system

## --- Retain deprecated functions with matlab-like syntax,
    # to avoid breakages in user scripts that may depend on them

    if VERSION>=v"1.0"
        """
        ```julia
        linspace(l::Number,u::Number,n::Number)
        ```

        Returns a linearly spaced array with `n` points between the starting
        bound `l` and ending bound `u`
        """
        function linspace(l::Number,u::Number,n::Number)
            return range(l,stop=u,length=n)
        end
        export linspace

        """
        ```julia
        (a,b) = linreg(x::AbstractVector, y::AbstractVector)
        ```

        Returns the coefficients for a simple linear least-squares regression of
        the form `y = a + bx`
        """
        linreg(x::AbstractVector, y::AbstractVector) = hcat(fill!(similar(x), 1), x) \ y
        export linreg

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
