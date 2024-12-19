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
<!-- TODO:  <21-11-24> -->
### Noise 
```@docs
GroundTruthGenerator

AdditiveNoise
apply(::AdditiveNoise, ::Any)
```
<!-- TODO:  <21-11-24> -->
### Microscopy transfer
<!-- TODO:  <21-11-24> -->

## Other Model Components
```@docs
DownSampling
apply(::DownSampling, ::Any)
```

## Ground Truth Models
### Fluorescent Microspheres (beads)
```@docs
bead
beads
synthetic_beads_image
```
