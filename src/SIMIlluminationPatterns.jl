# FIX: Generic methods without argument types should be made concrete, because if not, the error when supplying a wrong
# type of argument is not a MethodError for that function but another or worse, it can be a different error type
# entirely <12-12-23> 
module SIMIlluminationPatterns

using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1 true

using FFTW
using ImageCore
using Tullio

### source files

include("utils.jl")

# type system
include("common.jl")

# Helper API
include("sampled-illumination-pattern.jl")

# specific illumination patterns
include("harmonic.jl")
include("harmonic-component-separation.jl")

# synthetic data
include("synthetic-data/synthetic-data.jl")

export SampledIlluminationPattern, Synthetic
end
