using StaticArrays: SVector, SMatrix, @SMatrix, MArray

function mixin_matrix(m::SVector{N}, ϕ::SVector{N}) where {N}
    M = ones(Complex, N, 3)
    M[:, 2] = @. m / 2 * exp(im * ϕ)
    M[:, 3] = @. m / 2 * exp(-im * ϕ)

    return SMatrix{3,3}(M)
end

function mixin_matrix(m::AbstractVector, ϕ::AbstractVector)
    @assert length(m) == length(ϕ)
    N = length(m)
    mixin_matrix(SVector{N}(m), SVector{N}(ϕ))
end
mixin_matrix(m::AbstractVector, ϕ0::Real) = mixin_matrix(m, SVector{3}(ϕ0 .+ [2π * (n - 1) / 3 for n in 1:3]))
mixin_matrix(m::Real, ϕ) = mixin_matrix(@SVector(fill(m, 3)), ϕ)
# mixin_matrix(ϕ0::Real, N::Int) = mixin_matrix(ϕ0, ntuple(_ -> 1, Val(N)))

"""
    separation_matrix(ϕ, m)
    separation_matrix(ϕ0, m) # Assume equidistant phases
Construct a separation matrix `M` with `eltype(M) = Complex{T}` for given phase shifts `ϕ` and modulations `μ`.
"""
function separation_matrix(
    m::SVector{3}, ϕ::SVector{3}
)
    ϕ₁, ϕ₂, ϕ₃ = ϕ
    m₁, m₂, m₃ = m
    Δ = exp(im * (ϕ₂ - ϕ₁)) - exp(im * (ϕ₁ - ϕ₂)) - exp(im * (ϕ₃ - ϕ₁)) + exp(im * (ϕ₁ - ϕ₃)) + exp(im * (ϕ₃ - ϕ₂)) - exp(im * (ϕ₂ - ϕ₃))
    SMatrix{3,3}([
        (exp(im * (ϕ₂ - ϕ₃)) - exp(im * (ϕ₃ - ϕ₂)))
        (exp(im * (ϕ₃ - ϕ₁)) - exp(im * (ϕ₁ - ϕ₃)))
        (exp(im * (ϕ₂ - ϕ₁)) - exp(im * (ϕ₁ - ϕ₂)));
        (2(exp(im * ϕ₃) - exp(im * ϕ₂)) / m₁)
        (2(exp(im * ϕ₁) - exp(im * ϕ₃)) / m₂)
        (2(exp(im * ϕ₂) - exp(im * ϕ₁)) / m₃);
        (2(exp(-im * ϕ₂) - exp(-im * ϕ₃)) / m₁)
        (2(exp(-im * ϕ₃) - exp(-im * ϕ₁)) / m₂)
        (2(exp(-im * ϕ₁) - exp(-im * ϕ₂)) / m₃)
    ] ./ Δ)
end
separation_matrix(m::Real, ϕ) = separation_matrix(@SVector(fill(m, 3)), ϕ)
separation_matrix(m::SVector{3}, ϕ0::Real) = separation_matrix(m, SVector{3}(ϕ0 .+ [2π * (n - 1) / 3 for n in 1:3]))
