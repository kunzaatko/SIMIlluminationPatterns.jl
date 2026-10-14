using SIMIlluminationPatterns: SIMIlluminationPatterns as SIM_IP
using SIMIlluminationPatterns, Unitful

Δ = 61u"nm"
m, Pθ, ϕ, I0 = 0.5, [1 / 3, 1 / 4] / Δ, π / 2, 1 / Δ^2
@testset "Constructors" begin
    @test Harmonic(m, Pθ, ϕ, I0) isa Harmonic
    @testset "Constructor checks" begin
        @test_throws DomainError Harmonic(-1, Pθ, ϕ, I0)
        @test_throws AssertionError Harmonic(-1, [Pθ..., Pθ...], ϕ, I0)
        @test_throws DomainError Harmonic(m, Pθ, 3π, I0)
    end

    # @testset "Component Separation" begin
    #     # @testset "utils.jl" begin
    #     using SIMIlluminationPatterns: mixin_matrix, separation_matrix, mixin_components, separate_components
    #     @test mixin_matrix((0.5, 0.8, 1.0)) isa AbstractMatrix{<:Complex}
    #     @test mixin_matrix(0.0, 3) == mixin_matrix((0.0, 2π / 3, 4π / 3))
    #     @test mixin_matrix(0.5, 5) == mixin_matrix(0.5, (1, 1, 1, 1, 1))
    #
    #     @test separation_matrix(0.0, 3) isa AbstractMatrix{<:Complex}
    #     @test separation_matrix((0.5, 0.8, 1.0)) isa AbstractMatrix{<:Complex}
    #     @test separation_matrix((0.5, 0.8, 1.0), (0.3, 0.4, 0.2)) isa AbstractMatrix{<:Complex}
    #     @test_broken separation_matrix(0.5, (1, 1, 1, 1, 1)) isa AbstractMatrix{<:Complex}
    #
    #     M = mixin_matrix(0.0, 3)
    #     M_inv = separation_matrix(0.0, 3)
    #     comps = randn(10, 10, 3, 1)
    #     @test mixin_components(comps[:, :, :], M) isa AbstractArray{<:Complex,4}
    #     @test mixin_components(comps, M) == mixin_components(comps[:, :, :], M)
    #     raw = mixin_components(comps, M)
    #     @test separate_components(raw, M_inv) == separate_components(raw[:, :, :], M_inv)
    #     # FIX: BROKEN <18-12-24> 
    #     # @test comps ≈ separate_components(raw, M_inv)
    #
    #     @test separation_matrix((1.0, 2.0, 3.0), (2.0, 2.0, 2.0)) * transpose([1 1 1; exp(im) exp(2im) exp(3im); exp(-im) exp(-2im) exp(-3im)]) ≈ I(3)
    # end
end

@testset "intensity" begin
    h = Harmonic(m, Pθ, 0, I0)
    @test SIM_IP.intensity(h, [0u"nm", 0u"nm"]) * (Δ^2) == 0.75
end
