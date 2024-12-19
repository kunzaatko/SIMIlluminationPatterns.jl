# TODO: Add links to the docstring to the relevant ModelComponents in sections <19-12-24> 
@doc """
Serves to generate synthetic SIM data from a ground truth image.

The model of synthetic data generation is characterised by
    a list of [`ModelComponent`](@ref)s which can be for instance:
- **Noise** -- black current noise, background noise, photon shot noise (additive
    / multiplicative / etc. data dependent / independent)
- **Downsampling** -- downsampling the data so that the reconstruction could have the same size
    as the ground truth ([`DownSampling`](@ref Synthetic.DownSampling))
- **Illumination** -- illumination of the focal plane / sample volume
- **Optical transfer** -- transfer of the light through the optical system
"""
module Synthetic
using OffsetArrays, TransferFunctions, ImageFiltering, FFTW, ColorTypes
using ImageFiltering: mapwindow
using TransferFunctions: TransferFunction, NotImplementedError, SampledTransferFunction
using SIMIlluminationPatterns: GenericGrayImage, Length, IlluminationPattern, SampledIlluminationPattern
using Unitful: Quantity, 𝐋, Length
using Distributions: Normal, Uniform, mean, Poisson, Distribution

const PerLength = Quantity{<:Any,inv(𝐋)}

"""
    Synthetic.GroundTruthGenerator
A struct that generates a ground truth synthetic image based on stored latent variables.

# Implementation
Any type that implements this interface must define:
- `generate(gtg::GroundTruthGenerator)` -- generate the ground truth image
"""
abstract type GroundTruthGenerator end

"""
    Synthetic.ModelComponent
A struct that represents a component of the model of synthetic data generation.

# Implementation
Any type `A <: ModelComponent` that implements this interface must define:
- `apply(c::A, data; <keyword arguments>)` or `apply(c:A, data, ground_truth; <keyword arguments>)` -- transform `data` optionally depending on the `ground_truth`.
"""
abstract type ModelComponent end
apply(mc::ModelComponent, data, ground_truth; kwargs...) = apply(mc, data; kwargs...)

struct SyntheticDataModel
  components::Vector{<:ModelComponent}

  function SyntheticDataModel(components::Vector)
    @assert !isempty(components) "The number of components in the chain must be atleast 1"
    return new(components)
  end
end

# `params` - either latent variables that are used to generate the ground truth image or a 
function (chain::SyntheticDataModel)(ground_truth)
  # NOTE: Constructor guarantees that chain is non empty <14-12-24> 
  data = apply(chain.components[1], copy(ground_truth), ground_truth)
  for comp in chain.components[2:end]
    data = apply(comp, data, ground_truth)
  end
  return data
end

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
Simulate light transfer through the optical system by via a transfer function ([`SampledTransferFuction`](@exref))

# Examples
```jldoctest
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


```
"""
apply(ot::OpticalTransfer, data::AbstractArray) = TransferFunctions.apply(ot.transfer_function, data)

"""
    bead([T=Float64], d::Length, α::PerLength, (Δxy::Length,Δxy::Length))
    bead(d, α, Δxy::Length)
Generate a model of a bead with diameter `d` and pixelsizes `Δxy`.

# Arguments
- `α::PerLength`: evanescent wave attenuation constant
- `pixel_grid_length::Int = 10`: length of each pixel in the grid
- `peak_intensity = 1.0`: peak intensity value
-  `subpixel_shift =(0, 0)`
"""
function bead(T::Type{<:Real}, d::Length, α::PerLength, Δxy::NTuple{2,Length}; pixel_grid_length=10, peak_intensity=one(T), subpixel_shift=(0, 0))::OffsetMatrix{T}
  half_wh_px = d ./ (2 .* Δxy) .+ 1 .|> ceil # ½ width-height in pixels +1 for sub-pixel shift
  wh_px_g = (2 .* half_wh_px .+ 1) .* pixel_grid_length .|> Int # ½ width-height in grid, +1 for center
  grid = Matrix{T}(undef, wh_px_g)
  grid_subpixel_shift = subpixel_shift .* pixel_grid_length
  grid_center = wh_px_g ./ 2 .+ 1 / 2 .+ grid_subpixel_shift
  r_grid = map(CartesianIndices(grid)) do xy
    (Tuple(xy) .- grid_center) ./ pixel_grid_length .* Δxy |> splat(hypot)
  end # radii from the center
  supp_grid = r_grid .< (d / 2) # support of the bead
  z_grid = similar(r_grid, Union{Missing,Length}) # z-axes offset from the focal plane
  z_grid[supp_grid.==0] .= missing
  z_grid[supp_grid] .= (d / 2) .- sqrt.((d / 2)^2 .- r_grid[supp_grid] .^ 2)
  intensity_grid = map(z_grid) do z
    ismissing(z) ? zero(T) : exp(-α * z)
  end
  buf = Matrix{T}(undef, Int.(2 .* half_wh_px .+ 1))
  map!(buf, CartesianIndices(buf)) do ind
    indsx, indsy = map(Tuple(ind)) do i
      ((i-1)*pixel_grid_length+1):(i*pixel_grid_length)
    end
    mean(intensity_grid[indsx, indsy])
  end
  buf = OffsetArray(buf, -1 .* (half_wh_px .+ 1) .|> Int)
  buf .*= peak_intensity / buf[0, 0]
  return buf
end
bead(d::Length, α::PerLength, Δxy::NTuple{2,Length}; vargs...) = bead(Float64, d, α, Δxy; vargs...)
bead(T::Type{<:Real}, d::Length, α::PerLength, Δxy::Length; vargs...) = bead(T, d, α, (Δxy, Δxy); vargs...)
bead(d::Length, α::PerLength, Δxy; vargs...) = bead(Float64, d, α, Δxy; vargs...)

"""
    beads(T=Float64; <keyword arguments>)
Generate a synthetic microscopy image of fluorescent beads with realistic optical properties.

# Arguments
- `image_size::Tuple{Int,Int}=(1024, 1024)`: Size of the output image in pixels
- `N::Int=1000`: Number of beads to generate
- `pxsize::Tuple{Quantity,Quantity}=(30.5u"nm", 30.5u"nm")`: Pixel size in physical units
- `bead_radius::Quantity=0.05u"μm"`: Radius of each bead
- `min_distance::Quantity=0.06u"μm"`: Minimum allowed distance between bead centers
- `α_evanescent::Quantity=1/200u"nm"`: Evanescent field decay constant
"""
function beads()
  # TODO:  <18-11-24> 
end

"""
    synthetic_beads_image(T=Float64; <keyword arguments>)
Generate a synthetic microscopy image `Matrix{T}` of fluorescent beads with realistic optical properties.

# Arguments
- `image_size::Tuple{Int,Int}=(1024, 1024)`: Size of the output image in pixels
- `N::Int=1000`: Number of beads to generate
- `pxsize::Tuple{Quantity,Quantity}=(30.5u"nm", 30.5u"nm")`: Pixel size in physical units
- `bead_radius::Quantity=0.05u"μm"`: Radius of each bead
- `min_distance::Quantity=0.06u"μm"`: Minimum allowed distance between bead centers
- `α_evanescent::Quantity=1/200u"nm"`: Evanescent field decay constant
- `noise_model::Distribution=Normal{Float32}(2.0503677f-12, 4*0.005893613f0)`: Statistical model for image noise
- `optical_transfer_function=IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.9)`: Optical transfer function for microscope simulation

# Extended help
Generates a synthetic microscopy image of fluorescent beads, simulating:
1. Random placement of beads with minimum separation distance (`min_distance`)
2. Evanescent field illumination
3. Realistic noise
4. Optical transfer function effects

The beads are positioned randomly but maintain a minimum separation distance. Each bead's 
intensity is modulated by the evanescent field decay. The image is then convolved with 
the specified optical transfer function and noise is added according to the provided 
noise model.

# Examples
```julia
# Generate a synthetic image with default parameters
img = synthetic_beads_image()

# Generate a 512x512 image with fewer beads
img = synthetic_beads_image(Float32; image_size=(512,512), N=500)
```
"""
function synthetic_beads_image(
  T::Type{<:Number}=Float64;
  image_size=(1024, 1024),
  N=1000,
  pxsize=(30.5u"nm", 30.5u"nm"),
  bead_radius=0.05u"μm",
  min_distance=0.06u"μm",
  α_evanescent=1 / 200u"nm",
  noise_model=Normal{Float32}(2.0503677f-12, 4 * 0.005893613f0),
  illumination_patterns=[nothing],
  optical_transfer_function=IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.9) # Adams
)
  # Determine the beads positions
  beads_positions = []
  local pdists = map(lims -> Uniform(lims...), map(d -> (1, d), image_size))
  while length(beads_positions) != N
    proposal = rand(pdists[1]), rand(pdists[2])
    valid = all(map(x -> hypot(((proposal .- x) .* pxsize)...), beads_positions) .>= min_distance)
    if valid
      push!(beads_positions, proposal)
    end
  end

  # Adding the beads model
  buf = zeros(T, image_size)
  # TODO: padding should be determined based on the size of the bead (how far can it overstep the border?)
  buf_padded = padarray(buf, Fill(zero(T), (3, 3), (3, 3)))
  for p in beads_positions
    whole_px = round.(Int, p)
    subpixel_shift = p .- whole_px
    b = bead(bead_radius * 2, α_evanescent, pxsize; subpixel_shift)
    window = map(extrema(axes(b)), whole_px) do ex, px
      ex .+ px
    end
    buf_padded[window...] .+= b
  end
  buf_gt = buf_padded[1:image_size[1], 1:image_size[2]]

  # Adding noise
  if !isnothing(noise_model)
    buf_padded += OffsetArray(rand(noise_model, size(buf_padded)), -3, -3) # FIX: Better do this generically
  end

  # TODO: Change to IlluminationPattern instead of IlluminationPatternRealization <18-11-24> 
  # Illuminating the image
  @assert all(typeof(ip) <: Union{Nothing,SampledIlluminationPattern} for ip in illumination_patterns) "`illumination_patterns` must be `<:IlluminationPattern` or `nothing`"
  # TODO: Generalize for the size <18-11-24> 
  sampled_patterns = map(p -> p isa Nothing ? nothing : p(-2:1027, -2:1027), illumination_patterns)
  # TODO: This should be done in one go to avoid checking that it is a matrix... <18-11-24> 
  bufs_illuminated = map(sampled_patterns) do pat
    pat isa Matrix ? buf_padded .* OffsetMatrix(pat, -3, -3) : buf_padded
  end

  # Optical transfer
  # FIX: This malforms the padding. It should be done differently
  if !isnothing(optical_transfer_function)
    sampled_otf = SampledOTF(optical_transfer_function, pxsize)
    # FIX: This should really be done generically. It will fail if I do not have the predefined arguments
    bufs_padded = map(bufs_illuminated) do buf
      buf = ifft(otf(sampled_otf, (1030, 1030)) .* fft(buf))[4:end-3, 4:end-3]
      return real(buf)
    end
  else
    bufs_padded = bufs_illuminated
  end

  return [clamp.(buf, zero(T), maximum(buf))[1:image_size[1], 1:image_size[2]] for buf in bufs_padded], buf_gt
end

export SyntheticDataModel, apply
export DownSampling, AdditiveNoise, PhotonShotNoise, Illumination, OpticalTransfer
end
