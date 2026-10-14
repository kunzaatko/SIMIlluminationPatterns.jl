module SIMIlluminationPatterns
using LinearAlgebra
using FFTW, TransferFunctions, Unitful
using TransferFunctions: Frequency, Length, PixelSize, SpatialArray

include("utils.jl")

include("types.jl")

include("sampling.jl")

include("harmonic.jl")
include("harmonic-component-separation.jl")

# synthetic data
include("synthetic-data/synthetic-data.jl")

end
