using TransferFunctions
using TransferFunctions: sample_cells, SpatialArray

"""
    sample(ip::IlluminationPattern, sa::SpatialArray)
Sample the illumination pattern `ip` at the index locations of `sa`.

```jldoctest; setup=:(setup_params!()), filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> sa = TransferFunctions.SpatialArray(rand(10,10), 61u"nm");

julia> h = Harmonic(m, Pθ, 0, I0)
Harmonic{Float64}(m=0.5, Pθ=[0.00546 nm^-1, 0.0041 nm^-1], ϕ=0.0, I0=0.000269 nm^-2)

julia> SIMIlluminationPatterns.sample(h, sa)
10×10 SampledArray{Float64, Quantity{Int64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, 2, Matrix{Float64}, (61 nm, 61 nm)}:
 0.86838   0.86838   1.13162   1.13162   0.86838   …  1.13162   1.13162   0.86838   0.86838
 0.951824  1.1798    1.04818   0.820204  0.951824     1.04818   0.820204  0.951824  1.1798
 1.1798    0.951824  0.820204  1.04818   1.1798       0.820204  1.04818   1.1798    0.951824
 0.86838   0.86838   1.13162   1.13162   0.86838      1.13162   1.13162   0.86838   0.86838
 0.951824  1.1798    1.04818   0.820204  0.951824     1.04818   0.820204  0.951824  1.1798
 1.1798    0.951824  0.820204  1.04818   1.1798    …  0.820204  1.04818   1.1798    0.951824
 0.86838   0.86838   1.13162   1.13162   0.86838      1.13162   1.13162   0.86838   0.86838
 0.951824  1.1798    1.04818   0.820204  0.951824     1.04818   0.820204  0.951824  1.1798
 1.1798    0.951824  0.820204  1.04818   1.1798       0.820204  1.04818   1.1798    0.951824
 0.86838   0.86838   1.13162   1.13162   0.86838      1.13162   1.13162   0.86838   0.86838
```
"""
function sample(ip::IlluminationPattern, sa::SpatialArray)
    T = promote_type(eltype(ip), eltype(sa))
    out = similar(sa, T)
    copyto!(out, intensity.(ip, sample_cells(sa)))
    return out
end
