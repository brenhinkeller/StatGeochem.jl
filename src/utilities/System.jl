## --- Direct access to real libc system() function

    # Direct system() access without stripping special characters like Julia's
    # rather protective "run()". Learning new syntax for pipelining, etc. is
    # annoying. Is it really that common to be running shell commands from an
    # untrusted source?
    function system(cmdstr)
        return ccall((:system,), Int, (Cstring,), cmdstr)
    end
    export system

## --- Retain deprecated functions with matlab-like syntax, to avoid breakages

    if VERSION>=v"1.0"
        function linspace(l::Number,u::Number,n::Number)
            sp = (u-l)/(n-1)
            return l:sp:u
        end
        export linspace

        function repmat(A::AbstractArray, vert::Integer)
            return repeat(A, outer=vert)
        end
        function repmat(A::AbstractArray, vert::Integer, horiz::Integer)
            return repeat(A, outer=(vert, horiz))
        end
        export repmat

        # function contains(haystack::AbstractString, needle::Union{AbstractString,Regex,AbstractChar})
        #     return occursin(needle::Union{AbstractString,Regex,AbstractChar}, haystack::AbstractString)
        # end
        # export contains
    end

## --- End of File
