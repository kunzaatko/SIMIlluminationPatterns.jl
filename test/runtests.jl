using SIMIlluminationPatterns
using SIMIlluminationPatterns: SIMIlluminationPatterns as SIM
using Distributions, LinearAlgebra
using Test, Documenter, Aqua

macro no_error(ex)
    quote
        try
            $(esc(ex))
            true
        catch
            false
        end
    end
end

@testset "SIMIlluminationPatterns.jl" begin
    @testset "Code quality" begin
        @testset "Aqua.jl" begin
            if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
                Aqua.test_all(
                    SIMIlluminationPatterns;
                    ambiguities=false
                )
            else
                @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
            end
        end
        @testset "Ambiguities" begin
            @test length(Test.detect_ambiguities(SIMIlluminationPatterns)) == 0
        end
    end
    @testset "DocTests" begin
        # NOTE: Show for `Unitful.jl` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
        if !haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"
            # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
            # NOTE: When updating, must update also in `docs/make.jl` and  `test/fix_doctests.jl`<18-12-24> 
            DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(
                    using SIMIlluminationPatterns;
                    using SIMIlluminationPatterns.Synthetic;
                    using Distributions;
                    using TransferFunctions: TransferFunctions;
                    using TestImages;
                    filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
                    testimage.(filenames; download_only=false);
                    using Logging; # NOTE: This does not need to be in the `make.jl` of docs. We want `@warn ` to function there <19-12-24> 
                    Logging.disable_logging(Logging.Warn)
                ); recursive=true)
            doctest(SIMIlluminationPatterns)
        end
    end
    @testset "SampledIlluminationPattern" begin
        h = Harmonic(0.5, π / 4, 2 / 61u"nm", π)
        @test SampledIlluminationPattern(Float32, h, 61u"nm") isa SampledIlluminationPattern{Float32,2,Harmonic2D}
        @test SampledIlluminationPattern(h, (61u"nm", 61u"nm")) isa SampledIlluminationPattern{Float64,2,Harmonic2D}
        @test SampledIlluminationPattern(h, 61u"nm") isa SampledIlluminationPattern{Float64,2,Harmonic2D}

        sampledip = sampledip_f64 = SampledIlluminationPattern(h, 61u"nm")

        @test sampledip_f64(5.0, 3.1) isa Float64
        sampledip_f32 = SampledIlluminationPattern(Float32, h, 61u"nm")
        @test sampledip_f32(5.0, 3.1) isa Float32
        @test sampledip(CartesianIndex(1, 1)) isa Number

        @test sampledip((20, 20)) isa Matrix{Float64}
        @test size(sampledip((20, 20))) == (20, 20)

        img = zeros(10, 10)
        @test size(sampledip(img)) == (10, 10)
        @test sampledip(CartesianIndices(img))[1:9, 1:9] == sampledip(img)[2:10, 2:10]

    end
    @testset "Harmonic" begin
        Δxy = 61u"nm"
        m, θ, ν, ϕ = 0.5, π / 4, 2 / Δxy, π
        @testset "Constructors" begin
            @testset "Primary constructor checks" begin
                @test_throws DomainError Harmonic(-0.1, θ, ν, ϕ)
                # TODO: Make a test that checks that the function `@warn`s <11-12-23> 
                @test Harmonic(2.1, θ, ν, ϕ) isa SIM.IlluminationPattern # should only warn
                @test_throws DomainError Harmonic(m, -3π / 2, ν, ϕ)
                @test_throws DomainError Harmonic(m, 3π / 2, ν, ϕ)
                @test_throws DomainError Harmonic(m, θ, ν, -0.1)
                @test_throws DomainError Harmonic(m, θ, ν, 5π / 2)
            end
            @testset "Secondary constructors" begin
                # Harmonic(m::Real, θ::Real, ν::Frequency, ϕ::Real)
                @test @no_error Harmonic(m, θ, ν, ϕ)
                @test @no_error Harmonic(m, θ, Δxy / 2, ϕ)
                # Harmonic(m::Real, (kx, ky)::Tuple{Frequency,Frequency}, ϕ::Real)
                @test @no_error Harmonic(m, (2 / Δxy, 2 / Δxy), ϕ)

                @testset "Equivalences" begin
                    @test Harmonic(m, θ, Δxy / 2, ϕ) == Harmonic(m, θ, 2 / Δxy, ϕ)
                    @test Harmonic(m, (1 / (cos(θ) * Δxy), 1 / (sin(θ) * Δxy)), ϕ) ≈ Harmonic(m, θ, ν, ϕ)
                end
            end
            @testset "Component Separation" begin
                # @testset "utils.jl" begin
                using SIMIlluminationPatterns: mixin_matrix, separation_matrix, mixin_components, separate_components
                @test mixin_matrix((0.5, 0.8, 1.0)) isa Matrix{<:Complex}
                @test mixin_matrix(0.0, 3) == mixin_matrix((0.0, 2π / 3, 4π / 3))
                @test mixin_matrix(0.5, 5) == mixin_matrix(0.5, (1, 1, 1, 1, 1))

                @test separation_matrix(0.0, 3) isa Matrix{<:Complex}
                @test separation_matrix((0.5, 0.8, 1.0)) isa Matrix{<:Complex}
                @test separation_matrix((0.5, 0.8, 1.0), (0.3, 0.4, 0.2)) isa Matrix{<:Complex}
                @test_broken separation_matrix(0.5, (1, 1, 1, 1, 1)) isa Matrix{<:Complex}

                M = mixin_matrix(0.0, 3)
                M_inv = separation_matrix(0.0, 3)
                comps = randn(10, 10, 3, 1)
                @test mixin_components(comps[:, :, :], M) isa Array{<:Complex,4}
                @test mixin_components(comps, M) == mixin_components(comps[:, :, :], M)
                raw = mixin_components(comps, M)
                @test separate_components(raw, M_inv) == separate_components(raw[:, :, :], M_inv)
                # FIX: BROKEN <18-12-24> 
                # @test comps ≈ separate_components(raw, M_inv)

                @test separation_matrix((1.0, 2.0, 3.0), (2.0, 2.0, 2.0)) * transpose([1 1 1; exp(im) exp(2im) exp(3im); exp(-im) exp(-2im) exp(-3im)]) ≈ I(3)
            end
        end

        @testset "interfaces from Base" begin
            h1 = Harmonic(1, π / 4, 2 / 61u"nm", 0)
            h2 = Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0)

            # TODO: test `isequal` on missing values and equivalence operator on missing values <19-11-23> 
            # hmissing = Harmonic(missing, π / 4, 2 / 61u"nm", 0)

            @test h1 == h2
            @test hash(h1) == hash(h2)
            @test isequal(h1, h2)

        end
    end
    @testset "Synthetic data" begin
        using SIMIlluminationPatterns.Synthetic
        @test_throws AssertionError ForwardModel([])
        @testset "Ground Truth" begin
            using SIMIlluminationPatterns.Synthetic

            @test (Synthetic.bead(100u"nm", 30.5u"nm", peak_intensity=0.75) .<= 0.75) |> all
            @test_throws AssertionError Synthetic.bead(100u"nm", 30.5u"nm"; subpixel_shift=(-1.5, 0.7))
            @test_throws AssertionError Synthetic.bead(100u"nm", 30.5u"nm"; subpixel_shift=(1.5, 0.7))

            @test (Synthetic.bead(100u"nm", 30.5u"nm"; subpixel_shift=(0.5, 0.5)) .== reverse(Synthetic.bead(100u"nm", 30.5u"nm"; subpixel_shift=(-0.5, -0.5)))) |> all
            @test (2Synthetic.bead(100u"nm", 30.5u"nm"; peak_intensity=0.5) .== Synthetic.bead(100u"nm", 30.5u"nm")) |> all
        end
        @testset "Synthetic Model" begin
            @testset "Noise" begin
                data = ones(Float32, 100, 100)
                ground_truth = ones(Float32, 100, 100)

                poiss = PhotonShotNoise(0.1)
                @test poiss isa Synthetic.ModelComponent
                @test apply(poiss, data, ground_truth) != data

                gauss_additive = AdditiveNoise(Normal(0, 0.1))
                @test apply(gauss_additive, data) != data
            end

            data = ones(Float32, 100, 100)
            down_sampling = DownSampling(; ratio=2)
            @test apply(down_sampling, data) != size(data)
            datas = [data, data]
            @test size(apply(down_sampling, datas)[1]) != size(data)
        end
    end
end
