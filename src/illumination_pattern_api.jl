@doc raw"""
`N`-dimensional illumination pattern realization. You can fix the pixels sizes (`Δxy`) and sample the 
pattern on your sensor and optical system setup.
"""
struct IlluminationPatternRealization{T<:Real,N,IP<:IlluminationPattern{N}}
    pattern::IP
    "Pixel dimensions in the object space"
    Δxy::NTuple{N,Length}
end
const IPR = IlluminationPatternRealization

# TODO: Add docs and examples <30-10-23> 
@doc raw"""
    (ip::IP{N})(T::Type{<:Real}; Δxy)::IPR{T,N}
    (ip::IP{N})(; Δxy)::IPR{T,N}

Create an illumination pattern realization with pixel dimensions `Δxy`.

# Parameters
+ `Δxy`: can be a `NTuple{Length}` or a `Length`

```jldoctest
julia> h = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0)
Harmonic2D(m=1.0, θ=0.7853981633974483, ν=0.03278688524590164 nm^-1, ϕ=0.0)

julia> h(;Δxy=30.5u"nm")
Harmonic2D(1.0, 0.7853981633974483, 0.03278688524590164 nm^-1, 0.0){2}(Δxy = 30.5 nm) with eltype Float64

julia> h(Float16;Δxy=(33.5u"nm", 30.5u"nm"))
Harmonic2D(1.0, 0.7853981633974483, 0.03278688524590164 nm^-1, 0.0){2}(Δxy = (33.5 nm, 30.5 nm)) with eltype Float16
```
"""
function (ip::IP{N})(T::Type{<:Real}; Δxy::Union{NTuple{N,Length},Length}) where {N}
    Δxy = Δxy isa Length ? tuple(fill(Δxy, N)...) : Δxy
    IPR{T,N,typeof(ip)}(ip, Δxy)
end
(ip::IP{N})(; Δxy) where {N} = (ip)(Float64; Δxy)

# FIX: Why is this so slow? The slow part is evaluation of the harmonic <30-10-23> 
# NOTE: This would be a terrible problem for the optimization reconstruction...
# PERF:  
# JULIA:
#  @benchmark [1 + m / 2 * cos(2π * sum(sincos(θ) .* ν .* (y, x)) + ϕ) for y in 1:1024, x in 1:1024]
#  Time  (mean ± σ):   799.778 ms ±   6.472 ms  ┊ GC (mean ± σ):  3.58% ± 0.24%
#
# @benchmark Harmonic(0.5, π/4, 2/61u"nm", 0.0)(;Δxy = (30.5u"nm", 30.5u"nm"))(1:1024, 1:1024)
#  Time  (mean ± σ):   845.213 ms ±  10.466 ms  ┊ GC (mean ± σ):  4.44% ± 0.25%
#
# PYTHON:
# In [2]: %time A = illumination_pattern(0.6, 0.12, 0.5, 1, 1024)
# CPU times: user 18.5 ms, sys: 10.2 ms, total: 28.7 ms
# Wall time: 35.6 ms

# FIX: T should be used for generating of the given type this should be generalized in the pattern type <30-10-23> 
# FIX: Generate grid <30-10-23> 
# (ipr::IPR{T,2})(x::Real, y::Real) where {T} = ipr.pattern(x * ipr.Δxy[1], y * ipr.Δxy[2])
# (ipr::IPR{T,3})(x::Real, y::Real, z::Real) where {T} = ipr.pattern(x .* ipr.Δxy[1], y .* ipr.Δxy[2], z .* ipr.Δxy[3])
(ipr::IPR{T,N})(r::Vararg{T,N}) where {T,N} = ipr.pattern((r .* ipr.Δxy)...)
(ipr::IPR{T,N})(r::Vararg{Real,N}) where {T,N} = splat(ipr.pattern)(r .* ipr.Δxy)

# FIX: Generalize for any dimensions <30-10-23> 
(ipr::IPR{T,2})(xs::AbstractVector, ys::AbstractVector) where {T} = [ipr(x, y) for x in xs, y in ys]
(ipr::IPR{T,3})(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector) where {T} = [ipr(x, y, z) for x in xs, y in ys, z in zs]
(ipr::IPR{T,N})(size::NTuple{N,Integer}) where {T,N} = ipr(map(s -> 0:(s-1), size)...)

function Base.show(io::IO, ::MIME"text/plain", ipr::IPR{T,N}) where {T,N}
    show(io, MIME("text/plain"), ipr.pattern)
    print(io, "(Δxy = $(allequal(ipr.Δxy) ? ipr.Δxy[1] : ipr.Δxy)) with eltype $T")
end

# TODO: Add interface `map_frequencies` to general illumination pattern... <06-11-23> 
# It needs to be decided:
# - how to abstract over the image size
# - not every combination of image acquisitions with illumination patterns are possible to separate. Runtime errors or
# some use of the dispatch... Can a dispatch be made only for the types that are possible to separate and map. Perhaps
# a trait which would make it clear, whether a backward model is determined... Otherwise only possible reconstruction
# methods will be inversion methods.
