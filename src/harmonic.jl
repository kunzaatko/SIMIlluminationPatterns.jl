# TODO: 3D Harmonic has the same parameters as 2D harmonic, but can move z-sections. Implement. <30-10-23> 
# TODO: Generalize to multiple dimensions <30-10-23> 
# TODO: Implement constructors <11-10-23> 
# TODO: Check if all the methods that I want to include are included <11-10-23> 
# TODO: Add `@ref` links to `Frequency` and `Length` <12-10-23> 
@doc raw"""
Harmonic (sinusoidal) illumination pattern in the form
```math
    I(\vec{r})=1+{\frac{m}{2}}\cos\left(2π⋅ (kₓ⋅(\vec{r})ₓ + k_y ⋅ (\vec{r})_y) + \phi\right)
```

Parameters have types `Real`, `Frequency` or `Length` and denote:
* `m`: modulation factor
* `ν`: frequency (`\nu`)
* `λ`: wavelength (`\lambda`)
* `θ`: orientation angle (from the ``x``-axis) (`\theta`)
* `(kx, ky)`: wave vector (``(k_x, k_y) = (\sin(θ) ⋅ ν , \, \cos(θ) ⋅ ν)``)
* `(δx, δy)`: wavepeak vector (``(δ_x, δ_y) = (\sin(θ) ⋅ λ , \, \cos(θ) ⋅ λ)``)
* `ϕ`: phase offset (`\phi`)
* `Δxy` or `(Δx, Δy)`: ``x``-axis and ``y``-axis pixel sizes
"""
struct Harmonic{d} <: IP{d}
    "amplitude modulation"
    m::Real
    "orientation angle"
    θ::Real
    "frequency"
    ν::Frequency
    "phase offset"
    ϕ::Real

    function Harmonic{N}(m::Real, θ::Real, ν::Frequency, ϕ::Real) where {N}
        m <= one(m) || throw(DomainError(m, "amplitude modulation >1 would produce negative illumination intensities, which does not make sense"))
        m >= zero(m) || throw(DomainError(m, "amplitude modulation ∈(-1,0) is equivalent to shifting the phase offset by π (180°) and is not allowed"))
        zero(θ) <= θ < π || throw(DomainError(θ, "use an orientation ∈[0, π) (between 0° and 180°) (perhaps you should use: `mod(θ, π)`)"))
        zero(ϕ) <= ϕ < 2π || throw(DomainError(ϕ, "use a phase offset ∈[0, 2π) (between 0° and 360°) (perhaps you should use: `mod(ϕ, 2π)`)"))
        return new{N}(m, θ, ν, ϕ)
    end
end
const Harmonic2D = Harmonic{2}

# FIX: Is this the best way to set the default N? <30-10-23> 
Harmonic(a...) = Harmonic{2}(a...)

@doc raw"""
    Harmonic2D(m::Real, θ::Real, ν::Frequency, ϕ::Real)
    Harmonic2D(m::Real, (kx, ky)::Tuple{Frequency,Frequency}, ϕ::Real)
    Harmonic2D(m::Real, θ::Real, λ::Length, ϕ::Real)
    Harmonic2D(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
    Harmonic2D(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Length)
    Harmonic2D(m::Real, (kx, ky)::Tuple{Real,Real}, ϕ::Real, Δxy::Length)
    Harmonic2D(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
    Harmonic2D(m::Real, θ::Real, λ::Real, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})
    Harmonic2D(m::Real, (δx, δy)::Tuple{Real,Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

```@repl
Harmonic2D <: IlluminationPattern{2}
```
"""
Harmonic2D
Harmonic{N}(m::Real, θ::Real, λ::Length, ϕ::Real) where {N} = Harmonic{N}(m, θ, 1 / λ, ϕ)
# FIX: Check if hypot(k[1], k[2]) is correct <12-10-23> 
Harmonic{N}(m::Real, k::Tuple{Frequency,Frequency}, ϕ::Real) where {N} = Harmonic{N}(m, atan(ustrip(k[2]), ustrip(k[1])), hypot(k[1], k[2]), ϕ)
Harmonic{N}(m::Real, δ::Tuple{Length,Length}, ϕ::Real) where {N} = Harmonic{N}(m, atan(ustrip(δ[2]), ustrip(δ[1])), hypot(ustrip(δ[1]), ustrip(δ[2])), ϕ)

# Harmonic{N}(m::Real, θ::Real, λ::Length, ϕ::Real)
# Harmonic{N}(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
# Harmonic{N}(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Length)
# Harmonic{N}(m::Real, (kx, ky)::Tuple{Real,Real}, ϕ::Real, Δxy::Length)
# Harmonic{N}(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
# Harmonic{N}(m::Real, θ::Real, λ::Real, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})
# Harmonic{N}(m::Real, (δx, δy)::Tuple{Real,Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

function (h::Harmonic2D)(x::Length, y::Length)
    # TODO: Monomorphize the Length and Frequency <12-10-23> 
    # m, θ, ϕ = promote(h.m, h.θ, h.ϕ)
    return 1 + h.m / 2 * cos(2π * sum(sincos(h.θ) .* h.ν .* (y, x)) + h.ϕ)
end

# TODO: Abstract to multiple dimensions <30-10-23> 
# function (h::Harmonic{N})(v::Vararg{Length,N}) where {N}
#     m, θ, ϕ = promote(h.m, h.θ, h.ϕ)
#     return 1 + m / 2 * cos(2π * sum(sincos(θ) .* h.ν .* (y, x)) + ϕ)
# end

# TODO: Realizations of the Harmonic can be separate types created by supplying a pixel-size 
# (h::Harmonic)(Δxy)[1:400, 5:10], this type could implement abstract array interface. This could be done similarly with
# TODO: This should not be subtyped but should store the original type and should be monomorphized when sampled. This 
# would ensure that when there is a change in the attributes that it is still the most accurate representation <12-10-23> 
# the transfer functions <12-10-23> 

# TODO: Print the parameters that the Harmonic was created with with show. This would be done by storing them in the
# type itself <12-10-23> 
function Base.show(io::IO, ::MIME"text/plain", h::Harmonic{2})
    print(io, "Harmonic2D(", "m=", h.m, ", θ=", h.θ, ", ν=", h.ν, ", φ=", h.ϕ, ")")
end

const Harmonic3D = Harmonic{3}
@doc raw"""
```@example
Harmonic3D <: IlluminationPattern{3}
```
"""
Harmonic3D

function Base.show(io::IO, ::MIME"text/plain", h::Harmonic{3})
    print(io, "Harmonic3D(", "m=", h.m, ", θ=", h.θ, ", ν=", h.ν, ", φ=", h.ϕ, ")")
end

export Harmonic, Harmonic2D, Harmonic3D

# TODO: Implement Base.isequal <12-10-23> 
