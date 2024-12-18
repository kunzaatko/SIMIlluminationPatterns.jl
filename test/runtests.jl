using SIMIlluminationPatterns
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
    end
    @testset "DocTests" begin
        # NOTE: Show for `Unitful.jl` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
        if !haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"
            # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
            DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(using SIMIlluminationPatterns); recursive=true)
            doctest(SIMIlluminationPatterns)
        end
    end
    @testset "Harmonic" begin
        Δxy = 61u"nm"
        m, θ, ν, ϕ = 0.5, π / 4, 2 / Δxy, π
        @testset "Constructors" begin
            @testset "Primary constructor checks" begin
                @test_throws DomainError Harmonic(-0.1, θ, ν, ϕ)
                # TODO: Make a test that checks that the function `@warn`s <11-12-23> 
                @test Harmonic(2.1, θ, ν, ϕ) isa IlluminationPattern # should only warn
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
                using LinearAlgebra
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
                @test comps ≈ separate_components(raw, M_inv)

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
end
