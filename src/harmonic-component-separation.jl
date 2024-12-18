using ImageFiltering
# TODO: Should I overload `stack` to get the full orientation acquisition stack and check the correctness of the matrix
# separation? <25-11-24> 
# PERF: Could probably be memory efficient if it held only a reference to the array itself and was a view. But in the
# context of component separation, the matrices are created anyway. Could however be a good choice for
# `separate_components!` which needs to hold the  <25-11-24> 
# NOTE: There could also be a 3D component and the components under non-linear SIM have a higher indices <25-11-24> 
struct SeparatedComponent{T,N} <: AbstractArray{T,N}
    "Components in the Fourier domain"
    component::AbstractArray{T,N}
    "Direction of the shift. For ordinary SIM, this is -1, 0 and 1."
    shift_direction::Integer
end
Base.size(sc::SeparatedComponent) = size(sc.component)
Base.IndexStyle(::Type{<:SeparatedComponent}) = IndexLinear()
Base.getindex(sc::SeparatedComponent, i::Integer) = sc.component[i]
Base.setindex!(sc::SeparatedComponent, v, i::Integer) = (sc.component[i] = v)

function Base.show(io::IO, ::MIME"text/plain", c::SeparatedComponent{T}) where {T}
    print(join(size(c.component), "×"))
    print(" ShiftedComponent{$(typeof(c.component))}")
    print(" i=$(c.shift_direction) with eltype $(T)\n")
    Base.print_array(io, c.component)
end

# NOTE: https://discourse.julialang.org/t/how-to-create-struct-where-type-parameter-is-a-parametric-type-itself/101723/7
struct ShiftedComponent{T,N} <: AbstractArray{T,N}
    "Components in the Fourier domain"
    component::AbstractArray{T,N}
    "Frequency shift of the component in pixels"
    shift::NTuple{N,Real}
end
function Base.show(io::IO, ::MIME"text/plain", c::ShiftedComponent{T}) where {T}
    print(join(size(c.component), "×"))
    print(" ShiftedComponent{$(typeof(c.component))}")
    print(" shifted by δ=$(c.shift) with eltype $(T)\n")
    Base.print_array(io, c.component)
end
Base.size(sc::ShiftedComponent) = size(sc.component)
Base.IndexStyle(::Type{<:ShiftedComponent}) = IndexLinear()
Base.getindex(sc::ShiftedComponent, i::Integer) = sc.component[i]
Base.setindex!(sc::ShiftedComponent, v, i::Integer) = (sc.component[i] = v)

# TODO: Documentation <29-08-24> 
function mixin_matrix(
    T::Type{<:Real},
    ϕ::Tuple{TP,Vararg{TP,N}},
    μ::Tuple{TM,Vararg{TM,N}}=ntuple(_ -> 1, N + 1)
) where {N,TP<:Real,TM<:Real}
    M = ones(Complex{T}, N + 1, 3)
    M[:, 2] = collect(@. μ / 2 * exp(im * ϕ))
    M[:, 3] = collect(@. μ / 2 * exp(-im * ϕ))

    return M
end
# NOTE: There are unbound type parameters here, but it does not matter really <29-08-24> 
function mixin_matrix(ϕ::Tuple{TP,Vararg{TP,N}}, μ::Tuple{TM,Vararg{TM,N}}=ntuple(_ -> 1, N + 1)) where {N,TP<:Real,TM<:Real}
    mixin_matrix(promote_type(TP, TM), ϕ, μ)
end
function mixin_matrix(ϕ_0::Real, μ::Tuple{T,Vararg{T,N}}) where {N,T<:Real}
    mixin_matrix(ϕ_0 .+ Tuple(LinRange(0, 2π, N + 2)[begin:(end-1)]), μ)
end
mixin_matrix(ϕ_0::Real, N::Int) = mixin_matrix(ϕ_0, ntuple(_ -> 1, Val(N)))

# FIX: The math should not appear in the documentation <25-11-24> 
@doc raw"""
    separation_matrix(ϕ::NTuple{3,<:Real}, μ::NTuple{3,<:Real}=(1, 1, 1))
    separation_matrix(ϕ_0::Real,...) # Assume equidistant phases

Construct a separation matrix `M` with `eltype(M) = Complex{T}` for given phase shifts `ϕ` and modulations `μ`.

The matrix `M` is such that for acquisitions ``D_1(\tilde{k})``, ``D_2(\tilde{k})`` and ``D_3(\tilde{k})`` the
components are separated as
```math
    \begin{pmatrix}
    C_{0}(\tilde{k}) \\
    C_{-1}(\tilde{k}) \\
    C_{+1}(\tilde{k})
    \end{pmatrix} =
    \begin{pmatrix}
    H(\tilde{k})S(\tilde{k}) \\
    H(\tilde{k})S(\tilde{k} - \tilde{\nu}) \\
    H(\tilde{k})S(\tilde{k} + \tilde{\nu})
    \end{pmatrix} =
    \begin{pmatrix}
    1 & \tfrac{\mu_1}{2} e^{-i\phi_1} & \tfrac{\mu_1}{2} e^{i \phi_1} \\
    1 & \tfrac{\mu_2}{2} e^{-i\phi_2} & \tfrac{\mu_2}{2} e^{i \phi_2} \\
    1 & \tfrac{\mu_3}{2} e^{-i\phi_3} & \tfrac{\mu_3}{2} e^{i \phi_3}
    \end{pmatrix}^{-1}
    \begin{pmatrix}
    D_{1}(\tilde{k}) \\
    D_{2}(\tilde{k}) \\
    D_{3}(\tilde{k})
    \end{pmatrix}
```
"""
function separation_matrix(
    args...
)
    # PERF: Could be faster and more precise, if the matrix was created from the analytical inversion
    # FIX: This does not work for a non-square matrix. It needs to be handled differently for N != 3. Perhaps it will
    # lead to more variants of components -1 and 1 and then they will be averaged?? This could instead be
    # a pseudoinverse similar to linear regression and least squares <29-08-24> 
    return inv(mixin_matrix(args...))
end

# TODO: There should be a mutating version. Some times we only want the components and do not need to keep them. Mixin
# should be symmetric to this. <29-08-24> 
# TODO: Docs. Should include the fact that it is a 4 dim array by default <29-08-24> 
function mixin_components(f_imgs::AbstractArray{<:Number,3}, M::AbstractMatrix{<:Complex})
    nphases = size(M, 1)
    @assert mod(size(f_imgs, 3), nphases) == 0 """The number of images supplied must be a multiple of the number of phases.
    The complete stack of images should be of `size(f_imgs, 3)` == `norientations`×`nphases`"""
    comps = stack(Iterators.partition(eachslice(f_imgs, dims=3), nphases)) do single_orientation
        stack(row -> sum(single_orientation .* row), eachrow(M))
    end
    return comps
end
function mixin_components(f_imgs::AbstractArray{<:Number,4}, M::AbstractMatrix{<:Complex})
    mixin_components(reshape(f_imgs, size(f_imgs, 1), size(f_imgs, 2), :), M)
end

# TODO: Test these methods <25-11-24> 
mixin_components(f_imgs::AbstractArray{<:Number,3}, mixin_args...) = mixin_components(f_imgs, mixin_matrix(mixin_args...))
mixin_components(f_imgs::AbstractArray{<:Number,4}, mixin_args...) = mixin_components(f_imgs, mixin_matrix(mixin_args...))

# FIX: There is an issue with the separation of components if an OffsetArray is supplied.... <12-09-24> 
# TODO: There should be a mutating version. Some times we only want the components and do not need to keep them. Mixin
# should be symmetric to this. <29-08-24> 
# TODO: Docs. Should include the fact that it is a 4 dim array by default <29-08-24> 
function separate_components(f_imgs::AbstractArray{<:Number,3}, M_inv::AbstractMatrix{<:Complex})
    # NOTE: This currently only works if there are 3 phases... <29-08-24> 
    nphases = size(M_inv, 1)
    @assert mod(size(f_imgs, 3), nphases) == 0 """The number of images supplied must be a multiple of the number of phases.
    The complete stack of images should be of `size(f_imgs, 3)` == `norientations`×`nphases`"""
    comps = stack(Iterators.partition(eachslice(f_imgs, dims=3), nphases)) do single_orientation
        stack(row -> sum(single_orientation .* row), eachrow(M_inv))
    end
    comps = map(Iterators.cycle([0, -1, 1]), eachslice(comps, dims=(3, 4))) do d, comp
        SeparatedComponent(comp, d)
    end
    return comps
end
separate_components(f_imgs::AbstractArray{<:Number,4}, M_inv::AbstractMatrix{<:Complex}) = separate_components(reshape(f_imgs, size(f_imgs, 1), size(f_imgs, 2), :), M_inv)

# TODO: Test these methods <25-11-24> 
separate_components(f_imgs::AbstractArray{<:Number,3}, separation_args...) = separate_components(f_imgs, separation_matrix(separation_args...))
separate_components(f_imgs::AbstractArray{<:Number,4}, separation_args...) = separate_components(f_imgs, separation_matrix(separation_args...))


function shift_component(sc::SeparatedComponent{TC,N}, s_ip::SampledIlluminationPattern{TP,N,IP}) where {TC,TP,N,IP}
    Δ = sc.shift_direction .* δ(s_ip, size(sc))
    if all(iszero.(Δ)) # Zero index component
        return ShiftedComponent(sc.component, Δ)
    end
    padding_lower = Tuple(abs.(min.(zeros(Int, length(Δ)), Int.(sign.(Δ)) .* ceil.(Int, abs.(Δ)))))
    padding_upper = Tuple(abs.(max.(zeros(Int, length(Δ)), Int.(sign.(Δ)) .* ceil.(Int, abs.(Δ)))))
    comp_padded = padarray(sc.component, Fill(zero(eltype(sc.component)), padding_lower, padding_upper))
    # FIX: A problem with the offset getting lost. It should be possible to restore the offset.. This should be done by
    # specializing the shift function for the OffsetArray method. <25-11-24> 
    @show typeof(comp_padded), comp_padded.offsets
    shifted_comp = shift(Fourier(:fourier, true), fftshift(comp_padded), Δ) |> ifftshift
    return ShiftedComponent(shifted_comp, Δ)
end

# FIX: Is this correct for 3D SIM? <30-11-23> 
# function separate_components(siis::NTuple{3,IlluminatedImage{T,N,Harmonic{N}}}) where {T,N}
#     θs = map(sii -> sii.illumination_pattern.pattern.θ, siis)
#     allequal(θs) || ArgumentError("All images must have the same orientation. Given orientations $(θs).")
#
#     ϕs = map(sii -> sii.illumination_pattern.pattern.ϕ, siis)
#     ms = map(sii -> sii.illumination_pattern.pattern.m, siis)
#     M_inv = separation_matrix(T, ϕs, ms)
#
#     # NOTE: FFTW doesn't know how to work with Gray types <30-11-23> 
#     imgs = map(sii -> eltype(sii.img) isa AbstractGray ? gray.(sii.img) : sii.img, siis)
#     fft_lr_imgs = stack(map(img -> fft(img), imgs))
#
#     local C, x, y, i, j
#     @tullio C[x, y, i] := M_inv[i, j] * fft_lr_imgs[y, x, j]
#
#     return map(zip(eachslice(C, dims=ndims(C)), siis, (0, +1, -1))) do (component, sii, shift_ind)
#         # FIX: This is not right, we want to give the shift for a particular size not for making the size 2×... <05-12-23> 
#         ShiftedComponent(component, shift_ind .* δ(sii.illumination_pattern, size(sii.img) .* 2), sii.Δxy)
#     end |> splat(tuple)
# end

# TODO: Test <12-12-23> 
# NOTE: Taken from Distributions.jl <kunzaatko> 
for func in (:(==), :isequal, :isapprox)
    @eval function Base.$func(sc1::A, sc2::B; kwargs...) where {A<:ShiftedComponent,B<:ShiftedComponent}
        nameof(A) === nameof(B) || return false
        fields = fieldnames(A)
        fields === fieldnames(B) || return false

        for f in fields
            isdefined(sc1, f) && isdefined(sc2, f) || return false
            # perform equivalence check to support types that have no defined equality, such
            # as `missing`
            getfield(sc1, f) === getfield(sc2, f) || $func(getfield(sc1, f), getfield(sc2, f); kwargs...) || return false
        end

        return true
    end
end

# TODO: Test <12-12-23> 
# NOTE: Taken from Distributions.jl <kunzaatko> 
function Base.hash(sc::ShiftedComponent, h::UInt)
    hashed = hash(ShiftedComponent, h)
    hashed = hash(nameof(ShiftedComponent), hashed)

    for f in fieldnames(ShiftedComponent)
        hashed = hash(getfield(sc, f), hashed)
    end

    return hashed
end
