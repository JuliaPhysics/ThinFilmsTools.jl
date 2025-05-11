using ThinFilmsTools
using Aqua
using Test

@testset "ThinFilmsTools" begin
    @testset "Aqua" begin
        Aqua.test_all(ThinFilmsTools)
    end
    include("test3omega.jl")

    include("testTMMO.jl")

    include("testFitTMMO_DBRExp.jl")
    include("testFitTMMO_LayerLLL.jl")
    include("testFitTMMO_LayerLLLSpaceSolution.jl")
    include("testFitTMMO_MC1dExp.jl")
    include("testFitBBO.jl")

    include("testFitCurves.jl")
end
