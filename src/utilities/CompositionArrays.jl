# Arrays of compositions, using StructArrays.jl
struct CompositionArray{T,N,C,I} <: AbstractArray{T,N}
    data::StructArray{T,N,C,I}
end
CompositionArray(args...) = CompositionArray(StructArray(args...))
CompositionArray{T}(args...) where {T<:AbstractComposition} = CompositionArray(StructArray{T}(args...))
export CompositionArray

# Forward properties from wrapped StructArray
@inline Base.propertynames(x::CompositionArray) = propertynames(getfield(x, :data))
@inline Base.hasproperty(x::CompositionArray, key::Symbol) = hasproperty(getfield(x, :data), key)
@inline Base.getproperty(x::CompositionArray, key::Symbol) = getproperty(getfield(x, :data), key)

# Forward array interface to wrapped StructArray
@inline Base.setindex!(x::CompositionArray, args...) = setindex!(getfield(x, :data), args...)
@inline Base.getindex(x::CompositionArray, args...) = CompositionArray(getindex(getfield(x, :data), args...))
@inline Base.getindex(x::CompositionArray, inds::Vararg{Int}) = getindex(getfield(x, :data), inds...)
@inline Base.size(x::CompositionArray, args...) = size(getfield(x, :data), args...)
@inline Base.axes(x::CompositionArray, args...) = axes(getfield(x, :data), args...)

# Key-based indexing
@inline Base.keys(x::CompositionArray) = propertynames(x)
@inline Base.haskey(x::CompositionArray, key::Symbol) = hasproperty(x, key)
@inline Base.getindex(x::CompositionArray, key::Symbol) = getproperty(x, key)

# Major and trace elements: general fallback method
majors(x::CompositionArray) = filter(k->contains(String(k),"O"), keys(x))
traces(x::CompositionArray) = filter(k->!contains(String(k),"O"), keys(x))
# Major and trace elements: forward to Composition type if known
majors(x::CompositionArray{T}) where {T<:AbstractComposition} = majors(T)
traces(x::CompositionArray{T}) where {T<:AbstractComposition} = traces(T)


function normalize!(x::CompositionArray; anhydrous::Bool=false)
    for i in eachindex(x)
        normconst = zero(eltype(eltype(x)))
        for e in majors(x)
            if !isnan(x[e][i]) && (!anhydrous || !(e === :H2O || e === :CO2))
                normconst += x[e] / 100
            end
        end
        for e in traces(x)
            if !isnan(x[e][i])
                normconst += x[e] / 1_000_000
            end
        end
        for e in keys(x)
            x[e][i] /= normconst
        end
    end
    return x
end
function normalize!(x::CompositionArray{C}; anhydrous::Bool=false) where {T, C<:LogTraceComposition{T}}
    for i in eachindex(x)
        normconst = zero(T)
        for e in majors(x)
            if !isnan(x[e]) && (!anhydrous || !(e === :H2O || e === :CO2))
                normconst += x[e] / 100
            end
        end
        for e in traces(x)
            if !isnan(x[e])
                normconst += exp(x[e]) / 1_000_000
            end
        end
        lognormconst = log(normconst)
        for e in majors(x)
            x[e][i] /= normconst
        end
        for e in majors(x)
            x[e][i] -= lognormconst
        end
    end
    return x
end