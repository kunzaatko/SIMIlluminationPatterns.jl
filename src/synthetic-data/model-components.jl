# TODO: Could be a generic method resample <18-12-24> 
@doc raw"""
    Synthetic.DownSampling <: Synthetic.ModelComponent
Reduce the sampling of the data by a factor of `ratio` with the reduce function `reduce`

Fields: `ratio::Int`, `reduce::Function`

# Examples
```jldoctest; setup = :(using Statistics: mean)
julia> ds = DownSampling(3)
DownSampling(3) with reduce `mean`

julia> ds = DownSampling()
DownSampling(2) with reduce `mean`

julia> ds = DownSampling(reduce=maximum)
DownSampling(2) with reduce `maximum`
```
"""
Base.@kwdef struct DownSampling <: ModelComponent
    ratio::Int = 2
    reduce::Function = mean
end
DownSampling(ratio::Int) = DownSampling(ratio=ratio)
Base.show(io::IO, ::MIME"text/plain", ds::DownSampling) = print(io, "DownSampling(", ds.ratio, ") with reduce `", nameof(ds.reduce), "`")

@doc raw"""
    apply(ds::DownSampling, data)

# Examples
```jldoctest
julia> img = testimage("moonsurface.tiff");

julia> ds = DownSampling(2);

julia> size(img)
(256, 256)

julia> img_ds = apply(ds, img);

julia> size(img_ds)
(128, 128)
```
"""
function apply(ds::DownSampling, data::AbstractArray)
    indices = NTuple(map(axes(data)) do ax
        ax[1]:ds.ratio:ax[end]
    end
    )
    data = mapwindow(ds.reduce, data, fill(0:(ds.ratio-1), ndims(data)), border=Inner(), indices=indices)
    return data
end
apply(ds::DownSampling, datas::Vector{<:AbstractArray}) = map(data -> apply(ds, data), datas)

# TODO: Document the mathematical formula of the noise <19-12-24> 
"""
    PhotonShotNoise <: ModelComponent
Add Poisson noise simulating the photon shot noise (inherent to the quantum nature of light) to the data

See also [`AdditiveNoise`](@ref)

# Examples
```jldoctest
julia> noise_ps = PhotonShotNoise(0.1)
PhotonShotNoise(0.1)
```
"""
struct PhotonShotNoise <: ModelComponent
    α::Real
end
@doc raw"""
    apply(psn::PhotonShotNoise, data, ground_truth)
Add Poisson noise to `data` where the ``λ`` parameter of the Poisson distribution is proportional to the `ground_truth` (see [`PhotonShotNoise`](@ref))

# Examples
```jldoctest
julia> data = ground_truth = testimage("moonsurface.tiff");

julia> noise_ps = PhotonShotNoise(0.1);

julia> apply(noise_ps, data, ground_truth);
```
"""
function apply(psn::PhotonShotNoise, data, ground_truth)
    data = map(data, ground_truth) do d, gt
        d + rand(Poisson(gray(gt) * psn.α))
    end
    return data
end

"""
    AdditiveNoise{D<:Distribution} <: ModelComponent
Additive noise component of the synthetic data model.

See also [`PhotonShotNoise`](@ref)
# Examples
```jldoctest
julia> noise = AdditiveNoise(Normal(0, 0.1))
AdditiveNoise(Normal{Float64}(μ=0.0, σ=0.1))
```
"""
struct AdditiveNoise{D<:Distribution} <: ModelComponent
    dist::D
end
function Base.show(io::IO, ::MIME"text/plain", an::AdditiveNoise)
    print(io, "AdditiveNoise(")
    show(io, MIME("text/plain"), an.dist)
    print(io, ")")
end

@doc raw"""
    apply(noise::AdditiveNoise, data)
Add noise to data from `noise.dist`

# Examples
```jldoctest
julia> img = testimage("moonsurface.tiff");

julia> noise = AdditiveNoise(Normal(0, 0.1));

julia> img_δ = apply(noise, img);
```
"""
function apply(an::AdditiveNoise, data)
    data = data .+ rand(an.dist, size(data))
    return data
end

"""
    Illumination <: ModelComponent
Illuminate the image with a sampled illumination pattern.

# Examples
```jldoctest
julia> ip = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0);

julia> sampled_ip = SampledIlluminationPattern(ip, 61u"nm");

julia> ill = Illumination(sampled_ip)
Illumination(Harmonic2D(m=1.0, θ=0.785, ν=0.0328 nm^-1, ϕ=0.0)(Δxy = 61 nm) with eltype Float64)
```
"""
struct Illumination <: ModelComponent
    pattern::SampledIlluminationPattern
end
function Base.show(io::IO, ::MIME"text/plain", ill::Illumination)
    print(io, "Illumination(")
    show(io, MIME("text/plain"), ill.pattern)
    print(io, ")")
end
@doc raw"""
    apply(ill::Illumination, data::AbstractArray)
Illuminate the image `data` with the sampled illumination pattern `ill.pattern`.

# Examples
```jldoctest
julia> img = testimage("moonsurface.tiff");

julia> ip = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0);

julia> sampled_ip = SampledIlluminationPattern(ip, 61u"nm");

julia> ill = Illumination(sampled_ip);

julia> img_sim = apply(ill, img);
```
"""
apply(ill::Illumination, data::AbstractArray) = data .* ill.pattern(data)

"""
    OpticalTransfer <: `ModelComponent`
Simulate light transfer through the optical system by via a transfer function ([`SampledTransferFuction`](@extref TransferFunctions `TransferFunctions.SampledTransferFunction`))

# Examples
```jldoctest
julia> using TransferFunctions: BornWolf, SampledPSF

julia> bw = BornWolf(444u"nm", 1.4, 1.2);

julia> sampled_bw = SampledPSF(bw, 61u"nm");

julia> ot = OpticalTransfer(sampled_bw)
OpticalTransfer(SampledPSF{2, BornWolf{Float64}}(BornWolf{Float64}(444.0 nm, 1.4, 1.2), (61 nm, 61 nm), (0, 0)))
```
"""
struct OpticalTransfer <: ModelComponent
    transfer_function::SampledTransferFunction
end


"""
    apply(ot::OpticalTransfer, data::AbstractArray)
Simulate the light transfer by convolving with a transfer function `ot.transfer_function`.

# Examples
```jldoctest
julia> using TransferFunctions: BornWolf, SampledPSF

julia> bw = BornWolf(444u"nm", 1.4, 1.2);

julia> sampled_bw = SampledPSF(bw, 61u"nm");

julia> img = testimage("moonsurface.tiff");

julia> ot = OpticalTransfer(sampled_bw);

julia> apply(ot,img);
```
"""
apply(ot::OpticalTransfer, data::AbstractArray) = TF.apply(ot.transfer_function, data)
