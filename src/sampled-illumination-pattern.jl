raw"""
    SampledIlluminationPattern{T<:Real,N,IP<:IlluminationPattern{N}}
    SampledIlluminationPattern([T=Float64], ip::IlluminationPattern, Δxy)
`N`-dimensional sampled illumination pattern. You can fix the pixels sizes (`Δxy`) and sample the 
pattern on your sensor and optical system setup.

# Examples
```jldoctest;  filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
julia> h = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0)
Harmonic2D(m=1.0, θ=0.785, ν=0.0328 nm^-1, ϕ=0.0)

julia> sampled_ip = SampledIlluminationPattern(h, (30u"nm", 30u"nm"))
Harmonic2D(m=1.0, θ=0.785, ν=0.0328 nm^-1, ϕ=0.0)(Δxy = 30 nm) with eltype Float64

julia> sampled_ip(3, 4.5)
1.1048934112868387

julia> sampled_ip(CartesianIndex(1,1))
0.6126893879715403

julia> img = rand(10,10)
10×10 Matrix{Float64}:
[...]

julia> sampled_ip((10, 10))
10×10 Matrix{Float64}:
[...]

julia> sampled_ip(img)
10×10 Matrix{Float64}:
[...]

julia> sampled_ip(img) == sampled_ip(img)
true

julia> sampled_ip(CartesianIndices(img)) != sampled_ip(img) # BEWARE! This is not equivalent
true

julia> sampled_ip(CartesianIndices((0:9, 0:9))) == sampled_ip(img) # This is equivalent
true
```
"""
struct SampledIlluminationPattern{T<:Real,N,IP<:IlluminationPattern{N}}
    pattern::IP
    "Pixel dimensions in the object space"
    Δxy::NTuple{N,Length}
end
function SampledIlluminationPattern(
    T::Type{<:Real},
    ip::IP,
    Δxy::Length
) where {N,IP<:IlluminationPattern{N}}
    Δxy = ntuple(_ -> Δxy, Val(N))
    return SampledIlluminationPattern(T, ip, Δxy)
end
function SampledIlluminationPattern(
    T::Type{<:Real},
    ip::IP,
    Δxy
) where {N,IP<:IlluminationPattern{N}}
    return SampledIlluminationPattern{T,N,IP}(ip, Δxy)
end
SampledIlluminationPattern(ip, Δxy) = SampledIlluminationPattern(Float64, ip, Δxy)
const SampledIP = SampledIlluminationPattern

(ipr::SampledIP{T,N})(r::Vararg{Real,N}) where {T,N} = convert(T, ipr.pattern((r .* ipr.Δxy)...))
(ipr::SampledIP{T,N})(ind::CartesianIndex{N}) where {T,N} = ipr(Tuple(ind)...)

(ipr::SampledIP{T,N})(inds::CartesianIndices) where {T,N} = ipr.(inds)

# TODO: There should be a warning in the docs about method being different than the one when calling with `size` <18-12-24> 
(ipr::SampledIP{T,N})(size::NTuple{N,Integer}) where {T,N} = ipr(CartesianIndices(size) .- CartesianIndex(1, 1))
(ipr::SampledIP{T,N})(img::AbstractArray{TA,N}) where {TA,T,N} = ipr(size(img))

function Base.show(io::IO, ::MIME"text/plain", ipr::SampledIP{T,N}) where {T,N}
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

const GenericGrayImage{T<:Real,N} = AbstractArray{<:Union{T,AbstractGray{T}},N}

# TODO: This should probably be in a different package?! Perhaps, this could be a part of `SIMReconstruction` <30-11-23> 
# PERF: This could potentially be very bad on memory, since we are holding a potentially `MeasuredTransferFunction` with its data for every image...
# FIX: This could be solved by using `Ref` <30-11-23> 
# NOTE: https://discourse.julialang.org/t/how-to-create-struct-where-type-parameter-is-a-parametric-type-itself/101723/7
struct IlluminatedImage{T,N,IP<:IlluminationPattern{N}}
    img::GenericGrayImage{T,N}
    "Illumination pattern used in the acquisition"
    illumination_pattern::SampledIlluminationPattern{T,N}
    "Pixel dimensions in object space"
    Δxy::NTuple{N,Length}
end
IlluminatedImage(img::GenericGrayImage{T,N}, ip::IP{N}, Δxy::NTuple{N,Length}) where {T,N} = IlluminatedImage{T,N,typeof(ip)}(img, ip(T; Δxy), Δxy)
IlluminatedImage(img::GenericGrayImage{T,N}, ip::IP{N}, Δxy::Length) where {T,N} = IlluminatedImage(img, ip, Tuple(fill(Δxy, N)))
IlluminatedImage(img::GenericGrayImage{T,N}, ipr::SampledIP{N}) where {T,N} = IlluminatedImage(img, ipr.pattern, ipr.Δxy)

function Base.show(io::IO, ::MIME"text/plain", iimg::IlluminatedImage{T}) where {T}
    print(join(size(iimg.img), "×"))
    print(" ")
    show(io, MIME("text/plain"), iimg.illumination_pattern)
    print(" with eltype $(T) and Δxy = $(allequal(iimg.Δxy) ? iimg.Δxy[1] : iimg.Δxy)\n")
    Base.print_array(io, iimg.img)
end
