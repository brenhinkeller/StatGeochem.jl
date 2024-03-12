## --- Direct access to real libc system() function

    """
    ```julia
    system(cmdstr::AbstractString)
    ```

    Direct access to the command line through C's `system` function -- without
    stripping/sanitizing special characters,  in contrast to Julia's safer
    `run()` function. This allows pipelining, etc. in shell commands. Returns
    0 on success.

    ### Examples
    ```julia
    julia> system("pwd")
    /Users/foo/code/StatGeochem.jl
    0
    ```
    """
    function system(cmdstr::AbstractString)
        return ccall((:system,), Int, (Cstring,), cmdstr)
    end
    export system

## --- Some utilities for manipulating data types

    """
    ```julia
    materialize(x)
    ```
    Convert an array-like object to an materialized, actual allocated `Array`,
    and leaving other types unchanged.

    Unlike `collect`, will merely pass through an already-allocated `Array`
    without change, rather than allocating new memory and making a copy.

    ### Examples
    ```julia
    julia> StatGeochem.materialize(1:100)
    100-element Vector{Int64}:
       1
       2
       3
       ⋮
      99
     100

    julia> StatGeochem.materialize(5)
    5
    ```
    """
    materialize(x) = x
    materialize(x::Array) = x
    materialize(x::Number) = x
    materialize(x::Collection) = collect(x)

## --- End of File
