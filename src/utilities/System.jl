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

## --- Some utilities for manipulating data types

    materialize(x) = x
    materialize(x::Array) = x
    materialize(x::Number) = x
    materialize(x::AbstractArray) = collect(x)

## --- End of File
