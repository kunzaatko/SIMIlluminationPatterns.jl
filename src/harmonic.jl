# TODO: Implement constructors <11-10-23> 
# TODO: Check if all the methods that I want to include are included <11-10-23> 
@doc raw"""
    Harmonic(m::Real, θ::Real, ν::Frequency, ϕ::Real)
    Harmonic(m::Real, (kx, ky)::Tuple{Frequency,Frequency}, ϕ::Real)
    Harmonic(m::Real, θ::Real, λ::Length, ϕ::Real)
    Harmonic(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
    Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Length)
    Harmonic(m::Real, (kx, ky)::Tuple{Real,Real}, ϕ::Real, Δxy::Length)
    Harmonic(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
    Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})
    Harmonic(m::Real, (δx, δy)::Tuple{Real,Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

Harmonic (sinusoidal) illumination pattern in the form
```math
    I(\vec{r})=1+{\frac{m}{2}}\cos\left(2π⋅ (kₓ⋅\vec{r}ₓ + k_y ⋅ \vec{r}_y) + \phi\right)
```

Parameters have types `Real`, [`Frequency`](@ref Frequency) or [`Length`](@ref Length) and denote:
* `m`: modulation factor
* `ν`: frequency (`\nu`)
* `λ`: wavelength (`\lambda`)
* `θ`: orientation angle (from the ``x``-axis) (`\theta`)
* `(kx, ky)`: wave vector (``(k_x, k_y) = (\sin(θ) ⋅ ν , \, \cos(θ) ⋅ ν)``)
* `(δx, δy)`: wavepeak vector (``(δ_x, δ_y) = (\sin(θ) ⋅ λ , \, \cos(θ) ⋅ λ)``)
* `ϕ`: phase offset (`\phi`)
* `Δxy` or `(Δx, Δy)`: ``x``-axis and ``y``-axis pixel sizes
"""
struct Harmonic <: IP
    "amplitude modulation"
    m::Real
    "orientation angle"
    θ::Real
    "frequency"
    ν::Frequency
    "phase offset"
    ϕ::Real

    function Harmonic(m::Real, θ::Real, ν::Frequency, ϕ::Real)
        m <= one(m) || throw(DomainError(m, "amplitude modulation >1 would produce negative illumination intensities, which does not make sense"))
        m >= zero(m) || throw(DomainError(m, "amplitude modulation ∈(-1,0) is equivalent to shifting the phase offset by π (180°) and is not allowed"))
        zero(θ) <= θ < π || throw(DomainError(θ, "use an orientation ∈[0, π) (between 0° and 180°) (perhaps you should use: `mod(θ, π)`)"))
        zero(ϕ) <= ϕ < 2π || throw(DomainError(ϕ, "use a phase offset ∈[0, 2π) (between 0° and 360°) (perhaps you should use: `mod(ϕ, 2π)`)"))
        return new(m, θ, ν, ϕ)
    end
end

Harmonic(m::Real, θ::Real, λ::Length, ϕ::Real) = Harmonic(m, θ, 1 / λ, ϕ)
# FIX: Check if hypot(k[1], k[2]) is correct <12-10-23> 
Harmonic(m::Real, k::Tuple{Frequency,Frequency}, ϕ::Real) = Harmonic(m, atan(ustrip(k[2]), ustrip(k[1])), hypot(k[1], k[2]), ϕ)
Harmonic(m::Real, δ::Tuple{Length,Length}, ϕ::Real) = Harmonic(m, atan(ustrip(δ[2]), ustrip(δ[1])), hypot(ustrip(δ[1]), ustrip(δ[2])), ϕ)

# Harmonic(m::Real, θ::Real, λ::Length, ϕ::Real)
# Harmonic(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
# Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Length)
# Harmonic(m::Real, (kx, ky)::Tuple{Real,Real}, ϕ::Real, Δxy::Length)
# Harmonic(m::Real, (δx, δy)::Tuple{Length,Length}, ϕ::Real)
# Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})
# Harmonic(m::Real, (δx, δy)::Tuple{Real,Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

function (h::Harmonic)(x::Length, y::Length)
    # TODO: Monomorphize the Length and Frequency <12-10-23> 
    m, θ, ϕ = promote(h.m, h.θ, h.ϕ)
    return 1 + m / 2 * cos(2π * sum(sincos(θ) .* h.ν .* (y, x)) + ϕ)
end

# TODO: Realizations of the Harmonic can be separate types created by supplying a pixel-size 
# (h::Harmonic)(Δxy)[1:400, 5:10], this type could implement abstract array interface. This could be done similarly with
# TODO: This should not be subtyped but should store the original type and should be monomorphized when sampled. This 
# would ensure that when there is a change in the attributes that it is still the most accurate representation <12-10-23> 
# the transfer functions <12-10-23> 

# TODO: Print the parameters that the Harmonic was created with with show. This would be done by storing them in the
# type itself <12-10-23> 
function Base.show(io::IO, ::MIME"text/plain", h::Harmonic)
    print(io, "Harmonic(", "m=", h.m, ", θ=", h.θ, ", ν=", h.ν, ", φ=", h.ϕ, ")")
end

# TODO: Implement Base.isequal <12-10-23> 
