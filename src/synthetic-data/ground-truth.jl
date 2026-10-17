# TODO: Use the same syntax as `Augmentor.jl` <28-08-25> 

using TransferFunctions
using TransferFunctions.SyntheticData: bead

# TODO: Should be similar to the transfer functions model type that can be then sampled <06-03-25> 
abstract type SyntheticData end
# TODO: In the constructor, should check whether the positions are inside the frame <06-03-25> 

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
    optical_transfer_function=CircularPupilOTF(488u"nm", 1.4, 1.0, 0.9) # Adams
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
        b = bead(bead_radius * 2, pxsize; α=α_evanescent, position=subpixel_shift)
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
        # FIX: This should really be done generically. It will fail if I do not have the predefined arguments
        sampled_otf = otf(optical_transfer_function, pxsize, (1030, 1030))
        bufs_padded = map(bufs_illuminated) do buf
            buf = ifft(sampled_otf .* fft(buf))[4:end-3, 4:end-3]
            return real(buf)
        end
    else
        bufs_padded = bufs_illuminated
    end

    return [clamp.(buf, zero(T), maximum(buf))[1:image_size[1], 1:image_size[2]] for buf in bufs_padded], buf_gt
end
