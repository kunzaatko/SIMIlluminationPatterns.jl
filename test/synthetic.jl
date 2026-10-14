using TransferFunctions: SpatialArray
using SIMIlluminationPatterns.Synthetic

@test_throws AssertionError ForwardModel([])
@testset "Ground Truth" begin
    using SIMIlluminationPatterns.Synthetic

end
@testset "Synthetic Model" begin
    @testset "Noise" begin
        data = SpatialArray(ones(Float32, 100, 100), 50u"nm")
        ground_truth = SpatialArray(ones(Float32, 100, 100), 50u"nm")

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
