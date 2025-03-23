
## -- Pretty printing

    # Custom pretty printing for named tuples as datasets
    # TODO: possibly avoid blatant type piracy 🏴‍☠️ 🏴‍☠️ 🏴‍☠️
    function Base.show(io::IO, ::MIME"text/plain", x::NamedTuple)
        println(io, "NamedTuple with $(length(keys(x))) elements:")
        showcollection(io, x)
    end

    # Compact show methods for custom types
    function Base.show(io::IO, x::T) where {T<:AbstractComposition}
        k = first(keys(x))
        l = last(keys(x))
        print(io, "$T($k = $(x[k]), ... $l = $(x[l]))")
    end
    function Base.show(io::IO, x::CompositionArray{T,N}) where {T,N}
        print(io, "CompositionArray{$T,$N}$(size(x)) with $(length(keys(x))) elements $(first(keys(x))) ... $(last(keys(x)))")
    end

    # Verbose show methods for custom types
    function Base.show(io::IO, ::MIME"text/plain", x::C) where {C<:AbstractComposition}
        println(io, "$C with $(length(keys(x))) elements:")
        showcollection(io, x)
    end
    function Base.show(io::IO, ::MIME"text/plain", x::CompositionArray{T,N}) where {T,N}
        println(io, "CompositionArray{$T,$N}$(size(x)) with $(length(keys(x))) elements:")
        showcollection(io, x)
    end

    # Generic pretty printing for any collection that can be indexed by `keys`
    function showcollection(io::IO, x)
        i = 1
        l = max(length.(string.(keys(x)))...)
        for s in keys(x)
            t = typeof(x[s])
            sp = " "^(l-length(string(s)))
            print(io, "  $s$sp  = $t")
            if t<:Number
                print(io, "\t$(x[s])")
            elseif t<:AbstractRange
                print(io, "\t$(x[s])")
            elseif t<:AbstractArray
                print(io, size(x[s]))
                if length(x[s]) < 2
                    print(io, "\t[$(x[s])]")
                else
                    print(io, "\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:NTuple
                if length(x[s]) < 2
                    print(io, "\t[$(x[s])]")
                else
                    print(io, "\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:AbstractString
                if length(x[s]) < 50
                    print(io, "\t\"$(x[s])\"")
                else
                    print(io, "\t\"$(x[s][firstindex(x[s])+(1:50)])...")
                end
            end
            print(io, "\n")
            i += 1
            if i > 222
                print(io, ".\n.\n.\n")
                break
            end
        end
    end
