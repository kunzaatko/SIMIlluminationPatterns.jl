```@meta
CurrentModule = SIMIlluminationPatterns
CollapsedDocStrings = true
```

# Interface
Any illumination pattern is a subtype of the [`IlluminationPattern`](@ref)
```@docs
SIMIlluminationPatterns.IlluminationPattern
```

Such an illumination pattern can be sampled at points within the object plane/space by
[`intensity(IP, r)`](@ref intensity).
```@docs
SIMIlluminationPatterns.intensity(::IlluminationPattern, ::Any)
```
or as an illumination of a object density described by a [`SpatialArray`](@extref TransferFunctions.SpatialArray) by the
`SIMIlluminationPatterns.sample(IP, A)` method
```@docs
SIMIlluminationPatterns.sample(::IlluminationPattern, ::SpatialArray)
```
 
# Harmonic Illumination Pattern

The harmonic illumination pattern is represented by the [`Harmonic`](@ref) type
```@docs
Harmonic
```
