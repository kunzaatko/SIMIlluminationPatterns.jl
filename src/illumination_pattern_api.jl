using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

@doc raw"""
Super type for `N`-dimensional illumation patterns. Any particular illumination pattern is a subtype of this. See the 
documentation for [`IlluminationPattern` interface](@ref "Unified `IlluminationPattern` Interface") to see how to implement your own illumination pattern subtype.

See also: [`Harmonic2D`](@ref)
"""
abstract type IlluminationPattern{N} end
const IP{N} = IlluminationPattern{N}

(ip::IP{N})(::Vararg{Length,N}) where {N} = error("Not Implemented!")
# TODO: Add to interface docs <30-10-23> 
prefered_type(::IP{N}) where {N} = Float64

@doc raw"""
Super type for `N`-dimensional illumination pattern realization. You can fix the pixels sizes (`Δxy`) and sample the 
pattern on your senosor and optical system setup.
"""
struct IlluminationPatternRealization{T<:Real,N}
    "physical illumination pattern"
    pattern::IP{N}
    "pixel size"
    Δxy::NTuple{N,Length}
    function IlluminationPatternRealization{T,N}(ip::IP{N}, Δxy::NTuple{N,Length}) where {N,T}
        new{T,N}(ip, Δxy)
    end
end
const IPR{T,N} = IlluminationPatternRealization{T,N}

# TODO: Add docs and examples <30-10-23> 
function (ip::IP{N})(T::Type; Δxy::NTuple{N,Length}) where {N}
    IPR{T,N}(ip, Δxy)
end
(ip::IP{N})(; Δxy::NTuple{N,Length}) where {N} = (ip)(prefered_type(ip); Δxy)

# FIX: Why is this so slow? Compare with python implementation... The slow part is evaluation of the harmonic <30-10-23> 
# PERF:  
#  @benchmark [1 + m / 2 * cos(2π * sum(sincos(θ) .* ν .* (y, x)) + ϕ) for y in 1:1024, x in 1:1024]
#  Time  (mean ± σ):   799.778 ms ±   6.472 ms  ┊ GC (mean ± σ):  3.58% ± 0.24%
#
# @benchmark Harmonic(0.5, π/4, 2/61u"nm", 0.0)(;Δxy = (30.5u"nm", 30.5u"nm"))(1:1024, 1:1024)
#  Time  (mean ± σ):   845.213 ms ±  10.466 ms  ┊ GC (mean ± σ):  4.44% ± 0.25%

# FIX: T should be used for generating of the given type this should be generalized in the pattern type <30-10-23> 
# FIX: Generate grid <30-10-23> 
(ipr::IPR{T,2})(x::Real, y::Real) where {T} = ipr.pattern(x * ipr.Δxy[1], y * ipr.Δxy[2])
(ipr::IPR{T,3})(x::Real, y::Real, z::Real) where {T} = ipr.pattern(x .* ipr.Δxy[1], y .* ipr.Δxy[2], z .* ipr.Δxy[3])
(ipr::IPR{T,N})(r::Vararg{Real,N}) where {T,N} = splat(ipr.pattern)(r .* ipr.Δxy)

# FIX: Generalize for any dimensions <30-10-23> 
(ipr::IPR{T,2})(xs::AbstractVector, ys::AbstractVector) where {T} = [ipr(x, y) for x in xs, y in ys]
(ipr::IPR{T,3})(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector) where {T} = [ipr(x, y, z) for x in xs, y in ys, z in zs]
