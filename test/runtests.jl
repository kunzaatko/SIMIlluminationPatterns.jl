using SIMIlluminationPatterns
using Test
using Aqua

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
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SIMIlluminationPatterns; ambiguities=VERSION >= v"1.1")
    end
    @testset "Harmonic" begin
        Δxy = 61u"nm"
        m, θ, ν, ϕ = 0.5, π / 4, 2 / Δxy, π
        @testset "Constructors" begin
            @testset "Primary constructor checks" begin
                @test_throws DomainError Harmonic(-0.1, θ, ν, ϕ)
                @test_throws DomainError Harmonic(1.1, θ, ν, ϕ)
                @test_throws DomainError Harmonic(m, -π / 2, ν, ϕ)
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
                    @test_broken Harmonic(m, (2 / (cos(θ) * Δxy), 2 / (sin(θ) * Δxy)), ϕ) == Harmonic(m, θ, ν, ϕ)
                end
            end
        end
        @testset "types with Base.equiv" begin
            @test_broken Harmonic(1, π / 4, 2 / 61u"nm", 0) == Harmonic(1.0, π / 4, 2 / 61u"nm", 0.0)
        end
    end
end
