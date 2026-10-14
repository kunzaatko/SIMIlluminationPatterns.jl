# TODO: Resampling with different ratios if the different dimensions <08-04-25>
# TODO: Could be a generic method resample <18-12-24>
"""
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
Base.show(io::IO, ::MIME"text/plain", ds::DownSampling) = print(io, "DownSampling(", ds.ratio, ") with reduce `", nameof(ds.reduce), "`") # FIX: Should be handled generically <08-04-25> 

"""
    apply(ds::DownSampling, data)

# Examples
```jldoctest
julia> img = SpatialArray(testimage("moonsurface.tiff"), 61u"nm");

julia> ds = DownSampling(2);

julia> size(img)
(256, 256)

julia> img_ds = apply(ds, img);

julia> size(img_ds)
(128, 128)
```
"""
# FIX: Has to change the `Δ` of the `SpatialArray`. Should use `similar` function and then set the sampling. <08-04-25> 
function apply(ds::DownSampling, data::SpatialArray)
    indices = NTuple(map(axes(data)) do ax
        ax[1]:ds.ratio:ax[end]
    end)
    data = mapwindow(ds.reduce, data, fill(0:(ds.ratio-1), ndims(data)), border=Inner(), indices=indices)
    return data
end
apply(ds::DownSampling, datas::Vector{<:AbstractArray}) = map(data -> apply(ds, data), datas)

# TODO: Document the mathematical formula of the noise. It should be in the documentation rather than in the docstring of the function <19-12-24> 
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
"""
    apply(psn::PhotonShotNoise, data, ground_truth)
Add Poisson noise to `data` where the ``λ`` parameter of the Poisson distribution is proportional to the `ground_truth` (see [`PhotonShotNoise`](@ref))

# Examples
```jldoctest
julia> data = ground_truth = SpatialArray(testimage("moonsurface.tiff"), 61u"nm");

julia> noise_ps = PhotonShotNoise(0.1);

julia> apply(noise_ps, data, ground_truth);
```
"""
function apply(psn::PhotonShotNoise, data::SpatialArray, ground_truth::SpatialArray)
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
function Base.show(io::IO, ::MIME"text/plain", an::AdditiveNoise) # TODO: Generic show for all the model components <08-04-25> 
    print(io, "AdditiveNoise(")
    show(io, MIME("text/plain"), an.dist)
    print(io, ")")
end

"""
    apply(noise::AdditiveNoise, data)
Add noise to data from `noise.dist`

# Examples
```jldoctest
julia> img = SpatialArray(testimage("moonsurface.tiff"), 61u"nm");

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

julia> ill = Illumination(ip)
Illumination(Harmonic(m=1.0, θ=0.785, ν=0.0328 nm^-1, ϕ=0.0)(Δxy = 61 nm) with eltype Float64)
```
"""
struct Illumination <: ModelComponent
    pattern::IlluminationPattern
end
function Base.show(io::IO, ::MIME"text/plain", ill::Illumination)
    print(io, "Illumination(")
    show(io, MIME("text/plain"), ill.pattern)
    print(io, ")")
end
"""
    apply(ill::Illumination, data::SpatialArray)
Illuminate the image `data` with the sampled illumination pattern `ill.pattern`.

# Examples
```jldoctest
julia> img = SpatialArray(testimage("moonsurface.tiff"), 61u"nm");

julia> ip = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0);

julia> sampled_ip = SampledIlluminationPattern(ip, 61u"nm");

julia> ill = Illumination(sampled_ip);

julia> img_sim = apply(ill, img);
```
"""
apply(ill::Illumination, data::SpatialArray) = data .* ill.pattern(data)

"""
    OpticalTransfer <: `ModelComponent`
Simulate light transfer through the optical system by via a transfer function ([`TransferFunction`](@extref
TransferFunctions `TransferFunctions.TransferFunction`))

# Examples
```jldoctest
julia> using TransferFunctions: BornWolf

julia> bw = BornWolf(444u"nm", 1.4, 1.2);

julia> ot = OpticalTransfer(bw, 61u"nm")
OpticalTransfer(SampledPSF{2, BornWolf{Float64}}(BornWolf{Float64}(444.0 nm, 1.4, 1.2), (61 nm, 61 nm), (0, 0)))
```
"""
struct OpticalTransfer <: ModelComponent
    tf::TransferFunction
end

"""
    apply(ot::OpticalTransfer, data::AbstractArray)
Simulate the light transfer by convolving with a transfer function `ot.transfer_function`.

# Examples
```jldoctest
julia> using TransferFunctions: BornWolf, SampledPSF

julia> bw = BornWolf(444u"nm", 1.4, 1.2);

julia> img = SpatialArray(testimage("moonsurface.tiff"), 61u"nm");

julia> ot = OpticalTransfer(bw);

julia> apply(ot,img);
```
"""
apply(ot::OpticalTransfer, data::SpatialArray) = transfer(ot.tf, data)
