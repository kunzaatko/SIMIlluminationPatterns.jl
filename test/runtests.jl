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
    if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
        @testset "Code quality (Aqua.jl)" begin
            Aqua.test_all(
                SIMIlluminationPatterns;
                # ambiguities=VERSION >= v"1.1" ? (; broken=true) : false
            )
        end
    else
        @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
    end
    # NOTE: Show for `Unitful` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
    if !haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"
        @testset "DocTests" begin
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
                @test_throws DomainError Harmonic(1.1, θ, ν, ϕ)
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
