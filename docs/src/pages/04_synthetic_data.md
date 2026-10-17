```@meta
CurrentModule = SIMIlluminationPatterns.Synthetic
```

# Synthetic Microscopy Data
```@docs
Synthetic
```
<!-- TODO:  <21-11-24> -->

## Components of the Model

```@docs
ModelComponent
```

Below are the implemented `ModelComponent`s in the natural order of application.
<!-- TODO:  <21-11-24> -->
### Illumination
```@docs
Illumination
apply(::Illumination, ::SpatialArray)
```

### Noise 
```@docs
AdditiveNoise
apply(::AdditiveNoise, ::Any)
```

```@docs
PhotonShotNoise
apply(::PhotonShotNoise, ::SpatialArray, ::SpatialArray)
```
<!-- TODO:  <21-11-24> -->
### Microscopy transfer

```@docs
OpticalTransfer
apply(::OpticalTransfer, ::SpatialArray)
```
<!-- TODO:  <21-11-24> -->

### Other Model Components
```@docs
DownSampling
apply(::DownSampling, ::SpatialArray)
```

## Ground Truth Models
### Fluorescent Microspheres (beads)
```@docs
GroundTruthGenerator
synthetic_beads_image
```
