module SIMIlluminationPatterns

using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

### source files

# type system
include("common.jl")

# generic functions
include("illumination_pattern_api.jl")

# specific illumination patterns
include("harmonic.jl")


end
