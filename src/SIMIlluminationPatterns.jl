# FIX: Generic methods without argument types should be made concrete, because if not, the error when supplying a wrong
# type of argument is not a MethodError for that function but another or worse, it can be a different error type
# entirely <12-12-23> 
module SIMIlluminationPatterns

using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

using TransferFunctions
using FFTW
using ImageCore
using Tullio

### source files

# type system
include("common.jl")

# Helper API
include("illumination_pattern_api.jl")

# specific illumination patterns
include("harmonic.jl")
include("harmonic_component_separation.jl")

end
