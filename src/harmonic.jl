using Unitful: Length
using Unitful: 𝐋
using StaticArrays: @SVector, SVector
using DomainSets: volume, components, DomainSets, Rectangle
# TODO: Instead of testing every one of the constructors of `Harmonic` separately, I should define a set of the
# parameters and compute the corresponding other parameters from them in a previous example. It would be more
# explanatory for the user <19-12-24> 

"""
    Harmonic{T} <: IlluminationPattern{2}
2D Harmonic (sinusoidal) illumination pattern. 

The intensity is given by
```math
    I(r) = 1 + (m/2) * cos(2π * (k ⋅r + ϕ))
```
"""
struct Harmonic{T} <: IlluminationPattern{2}
    m::T
    Pθ::SVector{2,<:Freq{T}}
    ϕ::T
    I0::DimT{T,𝐋^-2}

    function Harmonic{T}(m::T, Pθ::SVector{2,<:Frequency{T}}, ϕ::T, I0::DimT{T,𝐋^-2}) where {T}
        m >= zero(T) || throw(DomainError(m, "amplitude modulation ∈(-1,0) is equivalent to shifting the phase offset by π (180°) and is not allowed"))
        zero(T) <= ϕ < 2π || throw(DomainError(ϕ, "use a phase offset ∈[0, 2π) (between 0° and 360°) (perhaps you should use: `mod(ϕ, 2π)`)"))

        return new{T}(m, Pθ, ϕ, I0)
    end
end
function Harmonic(m::Real, Pθ::AbstractVector{<:Frequency{<:Real}}, ϕ::Real, I0::DimT{<:Real,𝐋^-2})
    Pθ_unit = unit.(Pθ)
    I0_unit = unit(I0)
    @assert length(Pθ) == 2 "Wave vector must have length 2"
    m, Pθ1, Pθ2, ϕ, I0_val = promote(m, ustrip.(Pθ)..., ϕ, ustrip(I0))
    Pθ = SVector{2}([Pθ1, Pθ2] .* Pθ_unit)
    I0 = I0_val * I0_unit
    Harmonic{typeof(m)}(m, Pθ, ϕ, I0)
end

Base.eltype(h::Harmonic{T}) where {T} = T

cossin(θ) = reverse(sincos(θ))

δ(Pθ::SVector{2,<:Frequency}, A::SpatialMatrix) = Pθ .* sampling(A) .* size(A)
δ(h::Harmonic, A::SpatialMatrix) = δ(h.Pθ, A)
Pθ(θ::Real, δ::Real, Δ::PixelSize{2}) = SVector{2}(cossin(θ) ./ (Δ ./ δ))
θ(Pθ::SVector{2,<:Frequency}) = atan(Pθ[2], Pθ[1])
θ(h::Harmonic) = θ(h.Pθ)
ν(Pθ::SVector{2,<:Frequency}) = hypot(Pθ...)
ν(h::Harmonic) = ν(h.Pθ)
λ(x) = 1 / ν(x)

@inline function intensity(h::Harmonic, r::SVector{2,<:Length})
    return (1 - (h.m / 2) * cos(2π * dot(h.Pθ, r) + h.ϕ)) * h.I0
end

@inline function intensity(h::Harmonic, R::Rectangle)
    Dx, Dy = components(R)
    Δx, Δy = volume(Dx), volume(Dy)
    xc, yc = minimum(Dx) + (Δx / 2), minimum(Dy) + (Δy / 2)
    rc = SVector{2}([xc, yc])
    kx, ky = h.Pθ
    # NOTE: Julia definition of `sinc` is `sinc(x) == sin(π * x)/(π * x)`
    return h.I0 * Δx * Δy - (h.I0 * Δx * Δy) * sinc(kx * Δx) * sinc(ky * Δy) * (h.m / 2) * cos(2π * dot(h.Pθ, rc) + h.ϕ)
end

# δ(h::Harmonic) = cossin(h.θ) .* h.ν
# FIX: What does `size` represent? If it is the size of the image in pixels, it should only take Integer values. <25-11-24> 
# δ(hr::SampledIP{T,N,Harmonic{N}}, size::NTuple{N,Real}) where {T,N} = hr.Δxy .* δ(hr.pattern) .* size

# CHECK: This gives the shift δ that is determined by the size of the image... Is this correct? <05-12-23> 
# δ(hr::SampledIP{T,N,Harmonic{N}}, size::Real) where {T,N} = δ(hr, ntuple(_ -> size, Val(N)))

function Base.show(io::IO, ::MIME"text/plain", h::Harmonic)
    Base.showarg(io, h, true)
    print(io,"(", rounded(("m", h.m), ("Pθ", h.Pθ), ("ϕ", h.ϕ), ("I0", h.I0); sigdigits=3), ")")
end

export Harmonic
