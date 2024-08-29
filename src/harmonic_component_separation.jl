# FIX: This and a lot more can be made significantly better with StructArrays.jl <10-12-23> 

# TODO: This should probably be in a different package?! Perhaps, this could be a part of `SIMReconstruction` <30-11-23> 
# PERF: This could potentially be very bad on memory, since we are holding a potentially `MeasuredTransferFunction` with its data for every image... 
# FIX: This could be solved by using `Ref` <30-11-23> 
# NOTE: https://discourse.julialang.org/t/how-to-create-struct-where-type-parameter-is-a-parametric-type-itself/101723/7
# POLICY: This type can be constructed only internally with separation of components from `StructuredIlluminationImage`s
struct ShiftedComponent{T,N}
    "Components in the Fourier domain"
    component::AbstractArray{Complex{T},N}
    "Frequency shift of the component in pixels"
    shift::NTuple{N,Real}
    "Pixel dimensions in object space"
    Δxy::NTuple{N,Length} # FIX: This should be separated from the ShiftedComponent it is not necessary to determine it <01-12-23> 
end
function Base.show(io::IO, ::MIME"text/plain", c::ShiftedComponent{T}) where {T}
    print(join(size(c.component), "×"))
    print(" ShiftedComponent{$(typeof(c.component))}")
    print(" shifted by δ=$(c.shift) with eltype $(T) and Δxy = $(allequal(c.Δxy) ? c.Δxy[1] : c.Δxy)\n")
    Base.print_array(io, c.component)
end

"""
    separation_matrix(T::Type{<:Real}, ϕs::NTuple{3,<:Real}, m::NTuple{3,<:Real}=(1, 1, 1))

Construct a separation matrix `M` with `eltype(M) = Complex{T}` for given phase shifts `ϕs` and modulations `m`.
"""
function separation_matrix(T::Type{<:Real}, ϕs::NTuple{3,<:Real}, m::NTuple{3,<:Real}=(1, 1, 1))
    M = ones(Complex{T}, 3, 3)
    M[:, 2] = collect(@. m / 2 * exp(im * ϕs))
    M[:, 3] = collect(@. m / 2 * exp(-im * ϕs))
    # PERF: Could be faster, if the matrix was created as inv at the beginning
    return inv(M)
end
separation_matrix(ϕs::NTuple{3,<:Real}, m::NTuple{3,<:Real}=(1, 1, 1)) = separation_matrix(Float64, ϕs, m)

# FIX: Is this correct for 3D SIM? <30-11-23> 
function separate_components(siis::NTuple{3,IlluminatedImage{T,N,Harmonic{N}}}) where {T,N}
    θs = map(sii -> sii.illumination_pattern.pattern.θ, siis)
    allequal(θs) || ArgumentError("All images must have the same orientation. Given orientations $(θs).")

    ϕs = map(sii -> sii.illumination_pattern.pattern.ϕ, siis)
    ms = map(sii -> sii.illumination_pattern.pattern.m, siis)
    M_inv = separation_matrix(T, ϕs, ms)

    # NOTE: FFTW doesn't know how to work with Gray types <30-11-23> 
    imgs = map(sii -> eltype(sii.img) isa AbstractGray ? gray.(sii.img) : sii.img, siis)
    fft_lr_imgs = stack(map(img -> fft(img), imgs))

    local C, x, y, i, j
    @tullio C[x, y, i] := M_inv[i, j] * fft_lr_imgs[y, x, j]

    return map(zip(eachslice(C, dims=ndims(C)), siis, (0, +1, -1))) do (component, sii, shift_ind)
        # FIX: This is not right, we want to give the shift for a particular size not for making the size 2×... <05-12-23> 
        ShiftedComponent(component, shift_ind .* δ(sii.illumination_pattern, size(sii.img) .* 2), sii.Δxy)
    end |> splat(tuple)
end

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
