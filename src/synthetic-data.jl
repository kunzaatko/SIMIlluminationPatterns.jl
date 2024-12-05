module Synthetic
using OffsetArrays, TransferFunctions
using ImageFiltering
using TransferFunctions: TransferFunction
using SIMIlluminationPatterns: GenericGrayImage, Length, IlluminationPattern, IlluminationPatternRealization
using Unitful: Quantity, 𝐋, Length
using Distributions: Normal, Uniform, Sampleable, mean
using FFTW

const PerLength = Quantity{<:Any,inv(𝐋)}

# TODO: Never mind this, I will finish it later <18-11-24> 

struct GroundTruth{T<:Real,N}
    img::GenericGrayImage{T,N}
    # FIX: Is this necessary? When we would be creating a synthetic image from something like Lena, we do not need to define the pixel size <27-11-24> 
    Δxy::NTuple{N,Length}
end
GroundTruth(img::GenericGrayImage{T,N}, Δxy::Length) where {T,N} = GroundTruth(img::GenericGrayImage{T,N}, ntuple(_ -> Δxy, Val(N)))

# NOTE: Function should be able to generate a `Sampleable` from the intensity of a pixel or the position of the pixel
# within the image <18-11-24> 
const NoiseGenerator = Union{Function,Sampleable}

abstract type NoiseComponent end
struct NoiseSampleable{S<:Sampleable} <: NoiseComponent
    sampleable::S
end
function apply(i::GenericGrayImage, nc::NoiseSampleable{S}) where {S<:Sampleable}
    noisy = i.img .+ rand(nc.sampleable, size(i.img))
    return noisy
end

# apply(img::GenericGrayImage, nc::NoiseComponent) = apply!(copy(img), nc)
# struct NoiseComponent
#     # TODO:  <18-11-24> 
# end

apply

const ModelComponents{N} = Union{IlluminationPattern{N},TransferFunction{N},NoiseComponent}

function apply(img::GroundTruth, model::ModelComponents)
    apply(img, model)
end

function simulate(img::GroundTruth, model::Vector{ModelComponents})
    map(n -> apply(img, model[0:n]), 0:length(model))
end

"""
    bead([T=Float64], d::Length, α::PerLength, (Δxy::Length,Δxy::Length))
    bead(d, α, Δxy::Length)

Generate a model of a bead with diameter `d` and pixelsizes `Δxy`.

## Arguments
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
    beads(T::Type{<:Number}=Float64; kwargs...)

Generate a synthetic microscopy image of fluorescent beads with realistic optical properties.

# Arguments
- `T::Type{<:Number}=Float64`: The numeric type for the output image array

# Keyword Arguments
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
    synthetic_beads_image(T::Type{<:Number}=Float64; kwargs...)

Generate a synthetic microscopy image of fluorescent beads with realistic optical properties.

# Arguments
- `T::Type{<:Number}=Float64`: The numeric type for the output image array

# Keyword Arguments
- `image_size::Tuple{Int,Int}=(1024, 1024)`: Size of the output image in pixels
- `N::Int=1000`: Number of beads to generate
- `pxsize::Tuple{Quantity,Quantity}=(30.5u"nm", 30.5u"nm")`: Pixel size in physical units
- `bead_radius::Quantity=0.05u"μm"`: Radius of each bead
- `min_distance::Quantity=0.06u"μm"`: Minimum allowed distance between bead centers
- `α_evanescent::Quantity=1/200u"nm"`: Evanescent field decay constant
- `noise_model::Distribution=Normal{Float32}(2.0503677f-12, 4*0.005893613f0)`: Statistical model for image noise
- `optical_transfer_function=IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.9)`: Optical transfer function for microscope simulation

# Returns
A 2D array of type `T` containing the synthetic image with values clamped between 0 and the maximum intensity.

# Description
This function generates a synthetic microscopy image of fluorescent beads, simulating:
1. Random placement of beads with minimum separation distance
2. Evanescent field illumination
3. Realistic noise
4. Optical transfer function effects

The beads are positioned randomly but maintain a minimum separation distance. Each bead's 
intensity is modulated by the evanescent field decay. The image is then convolved with 
the specified optical transfer function and noise is added according to the provided 
noise model.

# Example
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
    @assert all(typeof(ip) <: Union{Nothing,IlluminationPatternRealization} for ip in illumination_patterns) "`illumination_patterns` must be `<:IlluminationPattern` or `nothing`"
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

end
