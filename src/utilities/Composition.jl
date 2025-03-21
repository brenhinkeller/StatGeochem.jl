# Abstract type for any Composition-type object and implementations thereof
abstract type AbstractComposition{T} end

# Variants where trace elements may or may not be stored in log concentration space
abstract type LinearTraceComposition{T} <: AbstractComposition{T} end
abstract type LogTraceComposition{T} <: AbstractComposition{T} end

# Generic pretty printing for any AbstractComposition
function Base.display(x::C) where {C<:AbstractComposition}
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

# Default methods which assume fields are elements, which will be used
# if a concrete type does not override with something more specific
Base.keys(x::C) where {C<:AbstractComposition} = fieldnames(C)
Base.haskey(x::C, key::Symbol) where {C<:AbstractComposition} = hasfield(C, key)
Base.getindex(x::AbstractComposition, key::Symbol) = getfield(x, key)
Base.zero(x::C) where {T, C<:AbstractComposition{T}} = C((zero(T) for _ in fieldnames(C))...,)
Random.rand(rng::AbstractRNG, ::Random.SamplerType{C}) where {T, C<:LinearTraceComposition{T}} = normalize(C((rand(rng, T)*20 for _ in majorelements(C))..., (exp(randn(rng, T)) for _ in traceelements(C))...,))
Random.rand(rng::AbstractRNG, ::Random.SamplerType{C}) where {T, C<:LogTraceComposition{T}} = normalize(C((rand(rng, T)*20 for _ in majorelements(C))..., (randn(rng, T) for _ in traceelements(C))...,))

# Major elements are assumed to be oxides, if concrete type does not override
majorelements(::Type{T}) where {T<:AbstractComposition} = filter(k->contains(String(k),"O"), fieldnames(T))
majorelements(::T) where {T<:AbstractComposition} = majorelements(T)
export majorelements

# Trace elements are assumed to be everything but oxides, if concrete type does not override
traceelements(::Type{T}) where {T<:AbstractComposition} = filter(k->!contains(String(k),"O"), fieldnames(T))
traceelements(::T) where {T<:AbstractComposition} = traceelements(T)
export traceelements

# Normalization
function normalize(x::C; anhydrous::Bool=false) where {T, C<:LinearTraceComposition{T}}
    normconst = zero(T)
    for e in majorelements(x)
        if !isnan(x[e]) && (!anhydrous || !(e === :H2O || e === :CO2))
            normconst += x[e] / 100
        end
    end
    for e in traceelements(x)
        if !isnan(x[e])
            normconst += x[e] / 1_000_000
        end
    end
    return C((x[e]/normconst for e in fieldnames(C))...,)
end
function normalize(x::C; anhydrous::Bool=false) where {T, C<:LogTraceComposition{T}}
    normconst = zero(T)
    for e in majorelements(x)
        if !isnan(x[e]) && (!anhydrous || !(e === :H2O || e === :CO2))
            normconst += x[e] / 100
        end
    end
    for e in traceelements(x)
        if !isnan(x[e])
            normconst += exp(x[e]) / 1_000_000
        end
    end
    lognormconst = log(normconst)
    return C((x[e]/normconst for e in majorelements(x))..., (x[e]-lognormconst for e in traceelements(x))...,)
end
export normalize

struct NCKFMASHTOtrace{T} <: LinearTraceComposition{T}
    SiO2::T
    TiO2::T
    Al2O3::T
    FeO::T
    MgO::T
    CaO::T
    Na2O::T
    K2O::T
    O2::T
    H2O::T
    P::T
    Rb::T
    Cs::T
    Sr::T
    Ba::T
    Sc::T
    V::T
    Cr::T
    Mn::T
    Co::T
    Ni::T
    La::T
    Ce::T
    Nd::T
    Sm::T
    Eu::T
    Gd::T
    Tb::T
    Dy::T
    Yb::T
    Lu::T
    Y::T
    Zr::T
    Hf::T
    Nb::T
    Ta::T
    Mo::T
    W::T
    Th::T
    U::T
end
export NCKFMASHTOtrace

struct NCKFMASHTOlogtrace{T} <: LogTraceComposition{T}
    SiO2::T
    TiO2::T
    Al2O3::T
    FeO::T
    MgO::T
    CaO::T
    Na2O::T
    K2O::T
    O2::T
    H2O::T
    P::T
    Rb::T
    Cs::T
    Sr::T
    Ba::T
    Sc::T
    V::T
    Cr::T
    Mn::T
    Co::T
    Ni::T
    La::T
    Ce::T
    Nd::T
    Sm::T
    Eu::T
    Gd::T
    Tb::T
    Dy::T
    Yb::T
    Lu::T
    Y::T
    Zr::T
    Hf::T
    Nb::T
    Ta::T
    Mo::T
    W::T
    Th::T
    U::T
end
export NCKFMASHTOlogtrace

# majorelements, traceelements and conversions
majorelements(::Type{<:Union{NCKFMASHTOtrace, NCKFMASHTOlogtrace}}) = (:SiO2, :TiO2, :Al2O3, :FeO, :MgO, :CaO, :Na2O, :K2O, :O2, :H2O)
traceelements(::Type{<:Union{NCKFMASHTOtrace, NCKFMASHTOlogtrace}}) = (:P, :Rb, :Cs, :Sr, :Ba, :Sc, :V, :Cr, :Mn, :Co, :Ni, :La, :Ce, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Yb, :Lu, :Y, :Zr, :Hf, :Nb, :Ta, :Mo, :W, :Th, :U)
NCKFMASHTOlogtrace(x::NCKFMASHTOtrace) = NCKFMASHTOlogtrace((x[e] for e in majorelements(x))..., (log(x[e]) for e in traceelements(x))...,)
NCKFMASHTOtrace(x::NCKFMASHTOlogtrace) = NCKFMASHTOtrace((x[e] for e in majorelements(x))..., (exp(x[e]) for e in traceelements(x))...,)

struct NCKFMASHTOCrtrace{T} <: LinearTraceComposition{T}
    SiO2::T
    TiO2::T
    Al2O3::T
    Cr2O3::T
    FeO::T
    MgO::T
    CaO::T
    Na2O::T
    K2O::T
    O2::T
    H2O::T
    P::T
    Rb::T
    Cs::T
    Sr::T
    Ba::T
    Sc::T
    V::T
    Mn::T
    Co::T
    Ni::T
    La::T
    Ce::T
    Nd::T
    Sm::T
    Eu::T
    Gd::T
    Tb::T
    Dy::T
    Yb::T
    Lu::T
    Y::T
    Zr::T
    Hf::T
    Nb::T
    Ta::T
    Mo::T
    W::T
    Th::T
    U::T
end
export NCKFMASHTOCrtrace

struct NCKFMASHTOCrlogtrace{T} <: LogTraceComposition{T}
    SiO2::T
    TiO2::T
    Al2O3::T
    Cr2O3::T
    FeO::T
    MgO::T
    CaO::T
    Na2O::T
    K2O::T
    O2::T
    H2O::T
    P::T
    Rb::T
    Cs::T
    Sr::T
    Ba::T
    Sc::T
    V::T
    Mn::T
    Co::T
    Ni::T
    La::T
    Ce::T
    Nd::T
    Sm::T
    Eu::T
    Gd::T
    Tb::T
    Dy::T
    Yb::T
    Lu::T
    Y::T
    Zr::T
    Hf::T
    Nb::T
    Ta::T
    Mo::T
    W::T
    Th::T
    U::T
end
export NCKFMASHTOCrlogtrace

# majorelements, traceelements and conversions
majorelements(::Type{<:Union{NCKFMASHTOCrtrace, NCKFMASHTOCrlogtrace}}) = (:SiO2, :TiO2, :Al2O3, :Cr2O3, :FeO, :MgO, :CaO, :Na2O, :K2O, :O2, :H2O)
traceelements(::Type{<:Union{NCKFMASHTOCrtrace, NCKFMASHTOCrlogtrace}}) = (:P, :Rb, :Cs, :Sr, :Ba, :Sc, :V, :Mn, :Co, :Ni, :La, :Ce, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Yb, :Lu, :Y, :Zr, :Hf, :Nb, :Ta, :Mo, :W, :Th, :U)
NCKFMASHTOCrlogtrace(x::NCKFMASHTOCrtrace) = NCKFMASHTOCrlogtrace((x[e] for e in majorelements(x))..., (log(x[e]) for e in traceelements(x))...,)
NCKFMASHTOCrtrace(x::NCKFMASHTOCrlogtrace) = NCKFMASHTOCrtrace((x[e] for e in majorelements(x))..., (exp(x[e]) for e in traceelements(x))...,)
