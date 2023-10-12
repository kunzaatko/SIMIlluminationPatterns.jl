```@meta
CurrentModule = SIMIlluminationPatterns
```
<!-- TODO: Use https://github.com/JuliaDocs/DocumenterCitations.jl for the citations <12-10-23> -->

# SIMIlluminationPatterns

Documentation for [SIMIlluminationPatterns](https://github.com/kunzaatko/SIMIlluminationPatterns.jl).

```@index
```

# Introduction

In the theory of image formation in a linear[^1] optical system, we model the data ``D(\vec{r})`` that is sensed by the "camera" (be
it a CCD, CMOS or a different type of sensor) as

```math
    D(\vec{r}) = E'(\vec{r}) + η(\vec{r}) = E(\vec{r}') \otimes H(\vec{r} - \vec{r}', \vec{r}) + η(\vec{r}),
```
where ``\otimes`` is the convolution operation, ``E(\vec{r})``[^2] denotes the intensity produced by the imaged object, disregarding noise,
``H(\vec{r} - \vec{r}', \vec{r})``[^3] is denotes the point spread function[^4] which is caused by the light passing through the 
objective and ``η(\vec{r})`` is the noise term. In a structured illumination microscopy (SIM) set-up, we further assume
that the emission is linearly dependent[^5] on the illumination intensity
```math
    E(\vec{r}) = S(\vec{r})I(\vec{r}).
```
This package is concerned with ``I(\vec{r})`` in this model. By far the most frequent type of illumination pattern used
in SIM is the [_harmonic_](harmonic.md) also called `sinusoidal`. Another type that is not yet implemented in this
package the random _speckle_ pattern used in, so called, _blind SIM_[^Mudry2012].

[^1]: A linear assumption is often necessary for any analysis, but it is approximately true for most optical systems
[^2]: "E" can stand for emission which is relevant for imaging "emitted" light, i.e. fluorescence microscopy
[^3]: For an _aplanic_ system, we can drop the ``\vec{r}'`` since every point in one focal plane is considered to have
      the same transfer function, i.e. ``H(\vec{r} - \vec{r}', \vec{r}) = H(\vec{r} - \vec{r}')``
[^4]: A package concerned with the transfer functions of optical systems is [`TransferFunctions.jl`](https://github.com/kunzaatko/TransferFunctions.jl)
[^5]: In fluorescence microscopy, the linearity of the relation is determined by the linearity of the relation between
      excitation and emission of the used fluorescent probes
[^Mudry2012]: Mudry, E., Belkebir, K., Girard, J., Savatier, J., Le Moal, E., Nicoletti, C., Allain, M., Sentenac, A., 2012. Structured illumination microscopy using unknown speckle patterns. Nature Photon 6, 312–315. https://doi.org/10.1038/nphoton.2012.83
