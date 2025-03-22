
## -- Import and extend Base.display for various custom types

    import Base.display

    # Custom pretty printing for named tuples as datasets
    # TODO: possibly avoid blatant type piracy 🏴‍☠️ 🏴‍☠️ 🏴‍☠️
    function display(x::NamedTuple)
        i = 1
        println("NamedTuple with $(length(keys(x))) elements:")
        l = max(length.(string.(keys(x)))...)
        for s in keys(x)
            t = typeof(x[s])
            sp = " "^(l-length(string(s)))
            print("  $s$sp  = $t")
            if t<:Number
                print("\t$(x[s])")
            elseif t<:AbstractRange
                print("\t$(x[s])")
            elseif t<:AbstractArray
                print(size(x[s]))
                if length(x[s]) < 2
                    print("\t[$(x[s])]")
                else
                    print("\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:NTuple
                if length(x[s]) < 2
                    print("\t[$(x[s])]")
                else
                    print("\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:AbstractString
                if length(x[s]) < 50
                    print("\t\"$(x[s])\"")
                else
                    print("\t\"$(x[s][firstindex(x[s])+(1:50)])...")
                end
            end
            print("\n")
            i += 1
            if i > 222
                print(".\n.\n.\n")
                break
            end
        end
    end

    # Generic pretty printing for any AbstractComposition
    function Base.display(x::C) where {C<:Union{AbstractComposition, CompositionArray}}
        i = 1
        println("$C with $(length(keys(x))) elements:")
        l = max(length.(string.(keys(x)))...)
        for s in keys(x)
            t = typeof(x[s])
            sp = " "^(l-length(string(s)))
            print("  $s$sp  = $t")
            if t<:Number
                print("\t$(x[s])")
            elseif t<:AbstractRange
                print("\t$(x[s])")
            elseif t<:AbstractArray
                print(size(x[s]))
                if length(x[s]) < 2
                    print("\t[$(x[s])]")
                else
                    print("\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:NTuple
                if length(x[s]) < 2
                    print("\t[$(x[s])]")
                else
                    print("\t[$(first(x[s])) ... $(last(x[s]))]")
                end
            elseif t<:AbstractString
                if length(x[s]) < 50
                    print("\t\"$(x[s])\"")
                else
                    print("\t\"$(x[s][firstindex(x[s])+(1:50)])...")
                end
            end
            print("\n")
            i += 1
            if i > 222
                print(".\n.\n.\n")
                break
            end
        end
    end
