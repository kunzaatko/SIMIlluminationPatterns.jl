if !isdefined(@__MODULE__, :LOADED_HASH) || LOADED_HASH != hash(readlines(@__FILE__))
    using SIMIlluminationPatterns
    using SIMIlluminationPatterns: SIMIlluminationPatterns as SIM_IP

    using Distributions, Unitful
    using TransferFunctions
    using TransferFunctions: TransferFunctions as TF

    using TestImages
    filenames = ["moonsurface", "mandril_gray"] # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
    testimage.(filenames; download_only=false)

    ENV["COLUMNS"] = 100
    ENV["LINES"] = 80

    setup_params!() = @eval begin
        Δ = 61u"nm"
        Δxy = 61u"nm"
        m, Pθ, ϕ, I0 = 0.5, [1 / 3, 1 / 4] / Δ, π / 2, 1 / Δ^2
        nothing
    end

    setup_makiemaestro!() = @eval begin
        using MakieMaestro
        using MakieMaestro.Recipes
        image = Recipes.image
        image! = Recipes.image!
        mosaic = Recipes.mosaic
        nothing
    end

    LOADED_HASH = hash(readlines(@__FILE__))

    return nothing

else
    @info "Skipping doctest setup... Same hash"
end
