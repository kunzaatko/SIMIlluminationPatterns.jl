using InterfaceFunctions
using Unitful: 𝐋

const DimT{T,D} = Union{Unitful.Quantity{T,D,U},Unitful.Level{L,S,Unitful.Quantity{T,D,U}} where {L,S}} where {U}
const Freq{T} = DimT{T,𝐋^-1}

"""
    IlluminationPattern{N}
Super type for `N`-dimensional illumination patterns. 
"""
abstract type IlluminationPattern{N} end

Broadcast.broadcastable(ip::IlluminationPattern) = Ref(ip)

"""
    @interface intensity(ip::IlluminationPattern{N}, r::SVector{N})
    @interface intensity(ip::IlluminationPattern{N}, R::HyperRectangle{<:SVector{N}})
Return the intensity of the illumination pattern `ip` at the location `r` ([`SVector{N}`](@extref SVector)) or under the
cell integration hyper-rectangle `R` ([`HyperRectangle{<:SVector{N}}`](@extref DomainSets.HyperRectangle)).

Every `IlluminationPattern` must implement these methods.

Method `intensity` when `r isa AbstractVector{<:Length}` has a default implementation converting `r` to `SVector{N}` and
calling the child implementation.
"""
@interface intensity(ip::IlluminationPattern, loc)
@inline function intensity(ip::IlluminationPattern{N}, r::AbstractVector{<:Length}) where {N}
    length(r) == N || throw(DimensionMismatch("Location must have a length of $N for getting the intensity of a an illumination pattern `<:IllumiatnionPattern{$N}`. Got `$r` of length $(length(r))"))
    intensity(ip, SVector{N}(r))
end
@interface Base.eltype(ip::IlluminationPattern) = Float64

# NOTE: Taken from Distributions.jl <kunzaatko> 
for func in (:(==), :isequal, :isapprox)
    @eval function Base.$func(ip1::A, ip2::B; kwargs...) where {A<:IlluminationPattern,B<:IlluminationPattern}
        nameof(A) === nameof(B) || return false
        fields = fieldnames(A)
        fields === fieldnames(B) || return false

        for f in fields
            isdefined(ip1, f) && isdefined(ip2, f) || return false
            # perform equivalence check to support types that have no defined equality, such
            # as `missing`
            getfield(ip1, f) === getfield(ip2, f) || $func(getfield(ip1, f), getfield(ip2, f); kwargs...) || return false
        end

        return true
    end
end

# NOTE: Taken from Distributions.jl <kunzaatko> 
function Base.hash(ip::IP, h::UInt) where {IP<:IlluminationPattern}
    hashed = hash(IlluminationPattern, h)
    hashed = hash(nameof(IP), hashed)

    for f in fieldnames(IP)
        hashed = hash(getfield(ip, f), hashed)
    end

    return hashed
end
