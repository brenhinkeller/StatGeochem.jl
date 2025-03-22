# Abstract type for any Composition-type object and implementations thereof
abstract type AbstractComposition{T} end

# Variants where trace elements may or may not be stored in log concentration space
abstract type LinearTraceComposition{T} <: AbstractComposition{T} end
abstract type LogTraceComposition{T} <: AbstractComposition{T} end

# Default methods which assume fields are elements, which will be used
# if a concrete type does not override with something more specific
Base.keys(x::C) where {C<:AbstractComposition} = fieldnames(C)
Base.haskey(x::C, key::Symbol) where {C<:AbstractComposition} = hasfield(C, key)
Base.getindex(x::AbstractComposition, key::Symbol) = getfield(x, key)

# Partial math interface, using generated functions so that we don't have to manually
# write out all the field names for every concrete subtype of AbstractComposition
@generated function Base.:*(x::C, n::Number) where {C<:LinearTraceComposition}
    result = :($C())
    for e in fieldnames(C)
        push!(result.args, :(x.$e * n))
    end
    return result
end
@generated function Base.:*(x::C, n::Number) where {C<:LogTraceComposition}
    result = :($C())
    for e in majorelements(C)
        push!(result.args, :(x.$e * n))
    end
    for e in traceelements(C)
        push!(result.args, :(x.$e + logn))
    end
    return Expr(
        :block,
        :(logn = log(n)),
        :(return $result),
    )
end
Base.:*(n::Number, x::AbstractComposition) = x * n          # Scalar multiplication is commutative
Base.:/(x::AbstractComposition, n::Number) = x * inv(n)     # Division by a scalar is multiplciation by multiplicative inverse
@generated function Base.:+(x::C, y::C) where {C<:LinearTraceComposition}
    result = :($C())
    for e in fieldnames(C)
        push!(result.args, :(x.$e + y.$e))
    end
    return result
end
@generated function Base.:+(x::C, y::C) where {C<:LogTraceComposition}
    result = :($C())
    for e in majorelements(C)
        push!(result.args, :(x.$e + y.$e))
    end
    for e in traceelements(C)
        push!(result.args, :(logaddexp(x.$e, y.$e)))
    end
    return result
end

# Generating zero and random compositions
@generated function Base.zero(::Type{C}) where {T, C<:AbstractComposition{T}}
    result = :($C())
    for e in fieldnames(C)
        push!(result.args, :(zero(T)))
    end
    return result
end
@generated function Random.rand(rng::AbstractRNG, ::Random.SamplerType{C}) where {T, C<:LinearTraceComposition{T}}
    result = :($C())
    for e in majorelements(C)
        push!(result.args, :(rand(rng, T)*20))
    end
    for e in traceelements(C)
        push!(result.args, :(exp(randn(rng, T))))
    end
    return :(normalize($result))
end
@generated function Random.rand(rng::AbstractRNG, ::Random.SamplerType{C}) where {T, C<:LogTraceComposition{T}}
    result = :($C())
    for e in majorelements(C)
        push!(result.args, :(rand(rng, T)*20))
    end
    for e in traceelements(C)
        push!(result.args, :(randn(rng, T)))
    end
    return :(normalize($result))
end

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
    return x/normconst
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
    return x/normconst
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
