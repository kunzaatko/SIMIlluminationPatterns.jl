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
using OffsetArrays, ImageFiltering, FFTW, ColorTypes, Unitful
using TransferFunctions: TransferFunctions as TF
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

struct ForwardModel
  components::Vector{<:ModelComponent}

  function ForwardModel(components::Vector)
    @assert !isempty(components) "The number of components in the model must be atleast 1"
    return new(components)
  end
end

# `params` - either latent variables that are used to generate the ground truth image or a 
function (chain::ForwardModel)(ground_truth)
  # NOTE: Constructor guarantees that chain is non empty <14-12-24> 
  data = apply(chain.components[1], copy(ground_truth), ground_truth)
  for comp in chain.components[2:end]
    data = apply(comp, data, ground_truth)
  end
  return data
end

include("model-components.jl")
include("ground-truth.jl")

export ForwardModel, apply
export DownSampling, AdditiveNoise, PhotonShotNoise, Illumination, OpticalTransfer
end
