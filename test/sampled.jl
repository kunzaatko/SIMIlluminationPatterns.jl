using SIMIlluminationPatterns, Unitful, TransferFunctions

Δ = 61u"nm"
m, Pθ, ϕ, I0 = 0.5, [1 / 3, 1 / 4] / Δ, π / 2, 1 / Δ^2

h = Harmonic(m, Pθ, ϕ, I0)

sa = SpatialArray(rand(10, 10), 61u"nm")

@test SIMIlluminationPatterns.sample(h, sa) isa SpatialArray
