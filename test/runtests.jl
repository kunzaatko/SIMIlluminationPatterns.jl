using SIMIlluminationPatterns
using SIMIlluminationPatterns: SIMIlluminationPatterns as SIM
using Distributions, LinearAlgebra
using Test, Documenter, Aqua

const run_all = isempty(ARGS) ? true : false

skip = Dict{String,Bool}(
    "aqua" => !haskey(ENV, "GITHUB_ACTIONS") && !haskey(ENV, "RUNTESTS_FULL"),
    "doctests" => !haskey(ENV, "RUNTESTS_FULL") && !(haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"),
    "ambiguities" => true # FIX: Fix the ambiguities <24-04-25> 
)

function should_test(arg::String)::Bool
    global run_all
    if run_all
        return !get(skip, arg, false)
    elseif arg in ARGS
        return true
    end
    return false
end

macro cond_testset(name, block)
    quote
        if should_test($name)
            @testset $name begin
                esc($block)
            end
        end
    end
end

@testset "SIMIlluminationPatterns.jl" begin

    @testset "Code quality" begin
        @cond_testset "aqua" begin
            Aqua.test_all(
                SIMIlluminationPatterns;
                ambiguities=false,
            )
        end

        @cond_testset "ambiguities" begin
            aqua_ambiguities = false
            if aqua_ambiguities
                Agua.test_ambiguities(SIMIlluminationPatterns)
            else
                @test length(Test.detect_ambiguities(SIMIlluminationPatterns)) == 0
            end
        end
    end

    @cond_testset "doctests" begin
        # FIX: When running locally, do not ask for SSH key password <10-12-23> 
        # NOTE: Show for `Unitful.jl` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
        # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
        # NOTE: When updating, must update also in `docs/make.jl` & `test/fix_doctests.jl` <18-12-24> 
        DocMeta.setdocmeta!(SIMIlluminationPatterns, :DocTestSetup, :(
                include(joinpath(@__DIR__, "doctestsetup.jl"));
                using Logging;
                # NOTE: Not necessary in `docs/make.jl`. `@warn` should work there <19-12-24>
                Logging.disable_logging(Logging.Warn)
            ); recursive=true)
        !haskey(ENV, "FIX_DOCTESTS") && @info "You can fix doctests by setting `ENV[\"FIX_DOCTESTS\"] = true`."
        doctest(SIMIlluminationPatterns; fix=ifelse(haskey(ENV, "FIX_DOCTESTS"), true, false))
    end

    @cond_testset "sampled" begin
        include("sampled.jl")
    end

    @cond_testset "harmonic" begin
        include("harmonic.jl")
    end

    # @cond_testset "synthetic" begin
    #     include("synthetic.jl")
    # end
end
