## --- Direct access to real libc system() function

    # Direct system() access without stripping special characters like Julia's
    # rather protective "run()". Learning new syntax for pipelining, etc. is
    # annoying. Is it really that common to be running shell commands from an
    # untrusted source?
    function system(cmdstr)
        return ccall((:system, "libc"), Int, (Cstring,), cmdstr)
    end
    export system

## --- End of File
