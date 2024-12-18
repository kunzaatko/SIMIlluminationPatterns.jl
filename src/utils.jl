# TODO: This should be in the SIMIlluminations package <18-11-24> 
using Interpolations: InterpolationType, BoundaryCondition, interpolate, extrapolate
using Tullio

@doc """
    TranslationAlgorithm

The abstract type for algorimthms for translating an array

# subtypes:
    - [`InterpolateExtrapolate`](@ref)
    - [`Fourier`](@ref)
"""
abstract type TranslationAlgorithm end
struct InterpolateExtrapolate <: TranslationAlgorithm
    intp::InterpolationType
    extp::Union{BoundaryCondition,Number}
    # Default values? Or defaults only for specific types?
    # intp=BSpline(Cubic(Flat(OnGrid()))), 
    # extp=zero(eltype(A))
end

# FIX: Why do we need `shifted` and what does it mean for the setup? (It is needed for component shifting, but not for
# the OTF, shifting) <10-12-23> 
struct Fourier <: TranslationAlgorithm
    domain::Union{Val{:fourier},Val{:spatial}}
    shifted::Bool
end
Fourier(dom::Symbol, args...) = dom ∈ (:fourier, :spatial) ? Fourier(Val(dom), args...) : throw(ArgumentError("domain must be ∈ (:fourier, :spatial). Got \'$dom\'"))

# TODO: Test <10-12-23> 
# FIX: This function should be possible to use for MeasuredPSF shifting <10-12-23> 
# IDEA!: Shift can function differently based on the type it gets to shift. e.g. "new" OTFEvalution (planned) should 
# handle the ifftshift and shift back and could have different defaults for the boundary conditions of the interpolation
# <24-10-23> 
# FIX: Should be dimension invariant <10-12-23> 
function shift(
    alg::InterpolateExtrapolate,
    A::AbstractMatrix,
    Δ::NTuple{2,<:Real}
    # IDEA: Add scale argument (as in interpolations... That is add argument for axes of A) <22-09-23> 
)
    # TODO: Check if bounds make sense for interpolation (if they are in the bounds of the data) <22-09-23> 

    A_intp = interpolate(A, alg.intp)
    A_extp = extrapolate(A_intp, alg.extp)

    out_x, out_y = (axes(A, i) .- Δ[i] for i in 1:2)
    return A_extp(out_x, out_y)
end

function shift(
    alg::Fourier,
    A::AbstractMatrix,
    Δ::NTuple{2,<:Real}
)
    # FIX: Is this correct for multiple dimensions? <10-12-23> 
    # PERF: Using `fftfreq` is 2× slower <10-12-23> 
    # xs, ys = fftshift.(fftfreq.(size(A)))
    xs, ys = fftfreq.(size(A))
    xs, ys = alg.shifted ? (fftshift(xs), fftshift(ys)) : (xs, ys)
    # xs, ys = map(size(A)) do s
    #     (Base.OneTo(s) .- 1 .- s ÷ 2) ./ s
    # end

    xexp = cispi.(-2Δ[1] * xs) # cis(x) === exp(im*x)
    yexp = cispi.(-2Δ[2] * ys)

    f_A = alg.domain == Val(:fourier) ? ifft(A) : fft(A)

    local fΔ_A, i, j
    @tullio fΔ_A[i, j] := f_A[i, j] * xexp[i] * yexp[j]

    Δ_A = alg.domain == Val(:fourier) ? fft(fΔ_A) : ifft(fΔ_A)
    return Δ_A
end
