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
@inline Base.view(x::CompositionArray, args...) = CompositionArray(view(getfield(x, :data), args...))
@inline Base.copy(x::CompositionArray) = CompositionArray(copy(getfield(x, :data)))

# Key-based indexing
@inline Base.keys(x::CompositionArray) = propertynames(x)
@inline Base.haskey(x::CompositionArray, key::Symbol) = hasproperty(x, key)
@inline Base.getindex(x::CompositionArray, key::Symbol) = getproperty(x, key)

# Major and trace elements: general fallback method
majorelements(x::CompositionArray) = filter(k->contains(String(k),"O"), keys(x))
traceelements(x::CompositionArray) = filter(k->!contains(String(k),"O"), keys(x))
# Major and trace elements: forward to Composition type if known
majorelements(x::CompositionArray{T}) where {T<:AbstractComposition} = majorelements(T)
traceelements(x::CompositionArray{T}) where {T<:AbstractComposition} = traceelements(T)

# Extend StatGeochemBase.renormalize! for CompositionArray{<:AbstractComposition}
function StatGeochemBase.renormalize!(x::CompositionArray{<:AbstractComposition}; anhydrous::Bool=false)
    for i in eachindex(x)
        xᵢ = x[i]
        c = anhydrous ? normconstanhydrous(xᵢ) : normconst(xᵢ)
        if !(c ≈ 1)
            x[i] = xᵢ/c
        end
    end
    return x
end

function partiallymix!(x::CompositionArray, mixingfraction::Number, ifirst=findfirst(!isnan, getfield(x,:data)), ilast=findlast(!isnan, getfield(x,:data)))
    @assert 0 <= mixingfraction <= 1 "Mixing fraction must be between 0 and 1"
    unmixingfraction = 1 - mixingfraction
    xfirst, xlast = x[ifirst], x[ilast]
    Δi = ilast-ifirst
    for i in ifirst:ilast
        fᵢ = (i - ifirst)/Δi
        mixᵢ = (fᵢ*xlast + (1-fᵢ)*xfirst)
        x[i] = unmixingfraction * x[i] + mixingfraction * mixᵢ
    end
    return x
end
export partiallymix!

# Extend NaNStatistics for CompositionArrays
function NaNStatistics.nanmean(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    μ = ntuple(i->T(nanmean(x[e[i]])), fieldcount(C))
    return C(μ)
end
function NaNStatistics.nanvar(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    σ² = ntuple(i->T(nanvar(x[e[i]])), fieldcount(C))
    return C(σ²)
end
function NaNStatistics.nanstd(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    σ = ntuple(i->T(nanstd(x[e[i]])), fieldcount(C))
    return C(σ)
end
function NaNStatistics.nansem(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    σ = ntuple(i->T(nansem(x[e[i]])), fieldcount(C))
    return C(σ)
end
function NaNStatistics.nancov(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    Σ = zeros(T, fieldcount(C), fieldcount(C))
    @inbounds for i in eachindex(e)
        for j in 1:i
            Σ[i,j] = Σ[j,i] = nancov(x[e[i]], x[e[j]])
        end
    end
    return Σ
end
function NaNStatistics.nancovem(x::CompositionArray{C}) where {T, C<:AbstractComposition{T}}
    e = fieldnames(C)
    Σ = zeros(T, fieldcount(C), fieldcount(C))
    @inbounds for i in eachindex(e)
        for j in 1:i
            Σ[i,j] = Σ[j,i] = nancovem(x[e[i]], x[e[j]])
        end
    end
    return Σ
end

## -- Distributions of compositions

# Abstract type for composition distributions
# Concrete subtypes are expected to have fields including 
#   `dist` - a Distributions.jl-compatible distribution of length equal to numbe of fields in the composition `C`
#   `buffer` - a Vector of length equal to numbe of fields in the composition `C`
abstract type CompositionDistribution{C} end


# Distributions.jl-style interface: forward to underlying distribution
Distributions.pdf(d::CompositionDistribution, x::AbstractVector) = pdf(d.dist, x)
Distributions.logpdf(d::CompositionDistribution, x::AbstractVector) = logpdf(d.dist, x)
function Distributions.pdf(d::CompositionDistribution{C}, x::C) where {C<:AbstractComposition}
    d.buffer .= ntuple(x)
    return pdf(d, d.buffer)
end
function Distributions.logpdf(d::CompositionDistribution{C}, x::C) where {C<:AbstractComposition}
    d.buffer .= ntuple(x)
    return logpdf(d, d.buffer)
end
Distributions.mean(d::CompositionDistribution{C}) where C = C(mean(d.dist))
Distributions.var(d::CompositionDistribution{C}) where C = C(var(d.dist))
Distributions.std(d::CompositionDistribution{C}) where C = C(sqrt.(var(d.dist)))
Distributions.cov(d::CompositionDistribution) = cov(d.dist)
Base.:(==)(a::CompositionDistribution, b::CompositionDistribution) = false
Base.:(==)(a::C, b::C) where {C<:CompositionDistribution} = isequal(a.dist, b.dist)

# The default concret type, based on multivariate normal
struct CompositionNormal{T, C<:AbstractComposition{T}, D<:MvNormal{T}} <: CompositionDistribution{C}
    dist::D
    buffer::Vector{T}
end
function CompositionNormal(::Type{C}, μ::AbstractVector, Σ::AbstractMatrix) where {T, C<:AbstractComposition{T}}
    @assert length(μ) == fieldcount(C)
    @assert size(Σ,1) == size(Σ,1) == fieldcount(C)
    d = MvNormal(μ, Σ)
    return CompositionNormal{T,C,typeof(d)}(d, zeros(T, fieldcount(C)))
end
function CompositionNormal(μ::C, Σ::AbstractMatrix) where {T, C<:AbstractComposition{T}}
    @assert size(Σ,1) == size(Σ,1) == fieldcount(C)
    d = MvNormal(collect(ntuple(μ)), Σ)
    return CompositionNormal{T,C,typeof(d)}(d, zeros(T, fieldcount(C)))
end
export CompositionNormal


# Interface for drawing Compositions from MvNormal distribution
# (using Random interface, rather than Distributions.jl interface)
Random.gentype(::Type{<:CompositionDistribution{C}}) where {C<:AbstractComposition} = C
function Random.rand(rng::AbstractRNG, d::Random.SamplerTrivial{<:CompositionNormal{T,C}}) where {T,C<:AbstractComposition{T}}
    rand!(rng, d.self.dist, d.self.buffer) 
    return renormalize(C(d.self.buffer))
end
