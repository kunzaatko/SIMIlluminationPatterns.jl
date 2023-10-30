```@meta
CurrentModule = SIMIlluminationPatterns
```

# SIMIlluminationPatterns

Documentation for [SIMIlluminationPatterns](https://github.com/kunzaatko/SIMIlluminationPatterns.jl).

```@index
```

# Introduction

In the theory of image formation in a linear[^1] optical system, we model the data ``D(\vec{x})`` that is sensed by the "camera" (be
it a CCD, CMOS or a different type of sensor) as

```math
    D(\vec{x}) = E'(\vec{x}) + η(\vec{x}) = E(\vec{y}) \otimes h(\vec{y}; \vec{x}) + η(\vec{x}),
```
where ``\otimes`` is the convolution operation i.e. , ``E(\vec{y})``[^2] denotes the intensity produced by the imaged object at the point ``\vec{y}``,
disregarding noise, ``h(\vec{y}; \vec{x})``[^3] denotes the point spread function (PSF) at the point ``\vec{x}``, evaluated at ``\vec{y}`` which is
caused by the light passing through the objective and ``η(\vec{x})`` is the noise term.

!!! info
    A package concerned with transfer functions of optical systems including the PSF is 
    [`TransferFunctions.jl`](https://github.com/kunzaatko/TransferFunctions.jl)


In the frequency domain ([Fourier Transform](https://en.wikipedia.org/wiki/Fourier_transform)), due to the [convolution
theorem of the Fourier Transform](https://en.wikipedia.org/wiki/Convolution_theorem), this transforms into
```math
    \tilde{D}(\vec{f}) = \tilde{E}(\vec{f}) \cdot H(\vec{f}) + N(\vec{f}),
```
where ``\vec{f}`` denotes the frequency and ``\tilde{D}``, ``\tilde{E}``, ``H``, called the Optical transfer function (OTF) and ``N`` are the Fourier 
transforms of ``D``, ``E``, ``h`` and ``η``, respectively. A particular property of the OTF of a diffraction limited
system is that it has a bounded support, which is a demonstration of the resolution limit. This can be intuited by the fact that frequencies outside of
the support `` \vec{\xi} \not\in \mathop{supp}(H)`` are not observed in the sensed data i.e. ``\tilde{D}(\vec{\xi}) = N(\vec{\xi})`` and thus any two points within distance of the reciprocal of those frequencies ``1/\|\vec{\xi}\|`` are observed as if produced by a single point.

In a regular structured illumination microscopy (SIM) set-up, we further assume
that the emission is linearly dependent[^4] on the illumination intensity
```math
    E(\vec{x}) = S(\vec{x})I(\vec{x}).
``` 
In the frequency domain, using the convolution theorem as above, we have 
```math
    \tilde{E}(\vec{f}) = \tilde{S}(\vec{f}) \otimes \tilde{I}(\vec{f}) = \int \tilde{S}(\vec{f} - \vec{\xi}) \cdot \tilde{I}(\vec{\xi}) \, d\vec{\xi},
```
which is useful since the higher frequencies that are outside of the OTF support ``\vec{\phi} \not\in \mathop{supp}(H)`` are now reflected in the sensed data 
``\tilde{D}`` as shifted into the support due to the convolution ``\mathop{supp}(H) \ni \vec{f} = \vec{\phi} - \vec{xi}``. This is the principle that enables to achieve higher resolution by using SIM (even overcoming the diffraction limit in resolution). The problem at hand is to recognize and separate the frequencies that are in-mixed in the sensed data ``\tilde{D}`` and move them to the correct places that they represent in the structure ``\tilde{S}``.

There are methods to induce a known non-linear relationship between the transfer function and the illumination intensity
such as saturated structured illumination microscopy, which can enable shifting higher frequencies into the sensed data
``D(\vec{x})``, than is possible to generate by a grating (due to diffraction effects), thus allowing to reach a higher
resolution.

This package is concerned with ``I(\vec{r})`` in this model. By far the most frequent type of illumination pattern used
in SIM is the [_harmonic_](harmonic.md) also called _sinusoidal_. Another type that is not yet implemented in this
package the random _speckle_ pattern used in, so called, _blind SIM_ [mudry2012](@cite).

[^1]: A linear assumption is often necessary for any analysis, but it is approximately true for most optical systems
[^2]: "E" can stand for emission which is relevant for imaging "emitted" light, i.e. fluorescence microscopy
[^3]: For an _aplanic_ system, we can drop the parameter of the PSF since every point in one focal plane is considered to have
      the same PSF, i.e. ``h(\vec{y}; \vec{x}) = h(\vec{y})``. This can be either a good or a bad approximation
      depending on the optical setup.
[^4]: In fluorescence microscopy, the linearity of the relation is determined by the linearity of the relation between
      excitation and emission of the used fluorescent probes
