module SIMIlluminationPatterns
using Reexport, Unitful

using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

abstract type IlluminationPattern end
const IP = IlluminationPattern

include("harmonic.jl")

@reexport using Unitful
export Harmonic

end
