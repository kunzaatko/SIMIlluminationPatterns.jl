using SIMIlluminationPatterns
using Test
using Aqua

@testset "SIMIlluminationPatterns.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SIMIlluminationPatterns)
    end
    # Write your tests here.
end
