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

    Harmonic(m::Real, θ::Real, ν::Frequency, ϕ::Real)
    Harmonic(m::Real, (ν_x, ν_y)::Tuple{Frequency,Frequency}, ϕ::Real)
    Harmonic(m::Real, (ν_x, ν_y)::Tuple{Real,Real}, ϕ::Real, Δxy::Length)
    Harmonic(m::Real, (ν_x, ν_y)::Tuple{Real,Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

    Harmonic(m::Real, θ::Real, λ::Length, ϕ::Real)
    Harmonic(m::Real, (λ_x, λ_y)::Tuple{Length,Length}, ϕ::Real)
    Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Length)
    Harmonic(m::Real, θ::Real, λ::Real, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})

    Harmonic(m::Real, (δ_x, δ_y)::Tuple{Real, Real}, size::Union{Tuple{Real, Real}, Real}, ϕ::Real, Δxy::Length)
    Harmonic(m::Real, (δ_x, δ_y)::Tuple{Real, Real}, size::Union{Tuple{Real, Real}, Real}, ϕ::Real, (Δx, Δy)::Tuple{Length,Length})


Parameters have types `Real`, `Frequency` or `Length` and denote:
+ `m`: modulation factor
+ `θ`: orientation angle (from the ``x``-axis) (`\theta`)
+ `ν`: frequency (`\nu`)
+ `(ν_x, ν_y)`: wave vector ('shift' in physical units) (``(ν_x, ν_y) = (\sin(θ) ⋅ ν , \, \cos(θ) ⋅ ν)``)
+ `(δ_x, δ_y)`: shift in pixels
+ `size`: size of acquisition image, from which the shift `(δ_x, δ_y)` is determined
+ `λ`: wavelength (`\lambda`)
+ `(λ_x, λ_y)`: wavepeak vector (``(λ_x, λ_y) = (\sin(θ) ⋅ λ , \, \cos(θ) ⋅ λ)``)
+ `ϕ`: phase offset (`\phi`)
+ `Δxy` or `(Δx, Δy)`: ``x``-axis and ``y``-axis pixel sizes

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
        -π < θ <= π || throw(DomainError(θ, "use an orientation ∈(-π, π) (between -180° and 180°) (perhaps you should use: `mod(θ + π, 2π) - π`)"))
        zero(ϕ) <= ϕ < 2π || throw(DomainError(ϕ, "use a phase offset ∈[0, 2π) (between 0° and 360°) (perhaps you should use: `mod(ϕ, 2π)`)"))
        return new{N}(m, θ, ν, ϕ)
    end
end

const Harmonic2D = Harmonic{2}
@doc raw"""
```@repl
Harmonic2D === Harmonic{2}
Harmonic{2} <: IlluminationPattern{2}
```
"""
Harmonic2D

const Harmonic3D = Harmonic{3}
@doc """
```@repl
Harmonic3D === Harmonic{3}
Harmonic{3} <: IlluminationPattern{3}
```
"""
Harmonic3D

@inline cossin(θ::Real) = reverse(sincos(θ))
θν(ν::Tuple{Frequency,Frequency}) = (atan(ν[2], ν[1]), hypot(ν...))
θν(ν::Tuple{Real,Real}, Δxy::Union{Tuple{Length,Length},Length}) = θν(ν ./ Δxy)
θν(λ::Tuple{Length,Length}) = (atan(λ[2], λ[1]), ν(prod(λ) / hypot(λ...)))
ν(λ::Length) = 1 / λ
θν(θ::Real, λ::Real, Δxy::Union{Tuple{Length,Length},Length}) = θν(λ .* cossin(θ) .* Δxy)
θν(δ::Tuple{Real,Real}, size::Union{Tuple{Real,Real},Real}, Δxy::Union{Tuple{Length,Length},Length}) = θν(δ ./ (size .* Δxy))

Harmonic(a...) = Harmonic{2}(a...)
Harmonic{N}(m::Real, ν::Tuple{Frequency,Frequency}, ϕ::Real) where {N} = Harmonic{N}(m, θν(ν)..., ϕ)
Harmonic{N}(m::Real, ν::Tuple{Real,Real}, ϕ::Real, Δxy::Union{Tuple{Length,Length},Length}) where {N} = Harmonic{N}(m, θν(ν, Δxy)..., ϕ)
Harmonic{N}(m::Real, λ::Tuple{Length,Length}, ϕ::Real) where {N} = Harmonic{N}(m, θν(λ)..., ϕ)
Harmonic{N}(m::Real, θ::Real, λ::Length, ϕ::Real) where {N} = Harmonic{N}(m, θ, ν(λ), ϕ)
Harmonic{N}(m::Real, θ::Real, λ::Real, ϕ::Real, Δxy::Union{Tuple{Length,Length},Length}) where {N} = Harmonic{N}(m, θ, θν(θ, λ, Δxy), ϕ)
Harmonic{N}(m::Real, δ::Tuple{Real,Real}, size::Union{Tuple{Real,Real},Real}, ϕ::Real, Δxy::Union{Tuple{Length,Length},Length}) where {N} = Harmonic{N}(m, θν(δ, size, Δxy)..., ϕ)

function (h::Harmonic{2})(x::Length, y::Length)
    # TODO: Monomorphize the Length and Frequency <12-10-23> 
    # m, θ, ϕ = promote(h.m, h.θ, h.ϕ)
    return 1 + h.m / 2 * cos(2π * sum(cossin(h.θ) .* h.ν .* (x, y)) + h.ϕ)
end

δ(h::Harmonic{N}) where {N} = cossin(h.θ) .* h.ν
δ(hr::IPR{T,N,Harmonic{N}}, size::NTuple{N,Real}) where {T,N} = hr.Δxy .* δ(hr.pattern) .* size
δ(hr::IPR{T,N,Harmonic{N}}, size::Real) where {T,N} = δ(hr, tuple(fill(size, N)...))

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
function Base.show(io::IO, ::MIME"text/plain", h::Harmonic{N}) where {N}
    print(io, "Harmonic$(N)D(", "m=", h.m, ", θ=", h.θ, ", ν=", h.ν, ", ϕ=", h.ϕ, ")")
end

export Harmonic, Harmonic2D, Harmonic3D
