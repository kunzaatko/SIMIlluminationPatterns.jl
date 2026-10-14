using Unitful

rounded(params::Vararg{Tuple{String, Any}}; kwargs...) = join([rounded(name, val; kwargs...) for (name, val) in params], ", ")
rounded(name::String, val; kwargs...) = name * "=" * rounded(val; kwargs...)
rounded(val::Quantity; kwargs...) = string(round(unit(val), val; kwargs...))
rounded(val; sigdigits=3, kwargs...) = string(round(val; sigdigits=sigdigits, kwargs...))
rounded(vals::AbstractVector; kwargs...)  = "[" * join([rounded(val; kwargs...) for val in vals], ", ") * "]"

#
# FIX: Ultimately, this should be removed as it should only be a utility function and the other functions should decide
# whether this to use the Fourier, i.e. `fftshift` version or the ordinary spatial version... <08-09-25> 
# TODO: This should be in the SIMIlluminations package? <18-11-24> 
# using Interpolations: InterpolationType, BoundaryCondition, interpolate, extrapolate
# using Tullio
#
# """
#     TranslationAlgorithm
#
# Algorithms for translating an array
# """
# abstract type TranslationAlgorithm end
#
# # TODO: Finish documentation <19-12-24> 
# """
#     InterpolateExtrapolate <: TranslationAlgorithm
# """
# struct InterpolateExtrapolate <: TranslationAlgorithm
#     intp::InterpolationType
#     extp::Union{BoundaryCondition,Number}
#     # Default values? Or defaults only for specific types?
#     # intp=BSpline(Cubic(Flat(OnGrid()))), 
#     # extp=zero(eltype(A))
# end
#
# # FIX: Why do we need `shifted` and what does it mean for the setup? (It is needed for component shifting, but not for
# # the OTF, shifting) <10-12-23> 
# """
#     Fourier <: TranslationAlgorithm
# """
# struct Fourier <: TranslationAlgorithm
#     domain::Union{Val{:fourier},Val{:spatial}}
#     shifted::Bool
# end
# Fourier(dom::Symbol, args...) = dom ∈ (:fourier, :spatial) ? Fourier(Val(dom), args...) : throw(ArgumentError("domain must be ∈ (:fourier, :spatial). Got \'$dom\'"))
#
# function shift(
#     alg::InterpolateExtrapolate,
#     A::AbstractMatrix,
#     Δ::NTuple{2,<:Real}
#     # IDEA: Add scale argument (as in interpolations... That is add argument for axes of A) <22-09-23> 
# )
#     # TODO: Check if bounds make sense for interpolation (if they are in the bounds of the data) <22-09-23> 
#
#     A_intp = interpolate(A, alg.intp)
#     A_extp = extrapolate(A_intp, alg.extp)
#
#     out_x, out_y = (axes(A, i) .- Δ[i] for i in 1:2)
#     return A_extp(out_x, out_y)
# end
#
# function shift(
#     alg::Fourier,
#     A::AbstractMatrix,
#     Δ::NTuple{2,<:Real}
# )
#     # FIX: Is this correct for multiple dimensions? <10-12-23> 
#     # PERF: Using `fftfreq` is 2× slower <10-12-23> 
#     # xs, ys = fftshift.(fftfreq.(size(A)))
#     xs, ys = fftfreq.(size(A))
#     xs, ys = alg.shifted ? (fftshift(xs), fftshift(ys)) : (xs, ys)
#     # xs, ys = map(size(A)) do s
#     #     (Base.OneTo(s) .- 1 .- s ÷ 2) ./ s
#     # end
#
#     xexp = cispi.(-2Δ[1] * xs) # cis(x) === exp(im*x)
#     yexp = cispi.(-2Δ[2] * ys)
#
#     f_A = alg.domain == Val(:fourier) ? ifft(A) : fft(A)
#
#     local fΔ_A, i, j
#     @tullio fΔ_A[i, j] := f_A[i, j] * xexp[i] * yexp[j]
#
#     Δ_A = alg.domain == Val(:fourier) ? fft(fΔ_A) : ifft(fΔ_A)
#     return Δ_A
# end
