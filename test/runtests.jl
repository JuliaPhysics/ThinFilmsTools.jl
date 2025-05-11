using ThinFilmsTools
using Aqua
using Test

@testset "ThinFilmsTools" begin
    @testset "Aqua" begin
        Aqua.test_all(ThinFilmsTools)
    end
    @testset "test3omega" begin
        include("test3omega.jl")
    end

    @testset "testTMMO" begin
        include("testTMMO.jl")
    end

    @testset "testFitTMMO_DBRExp" begin
        include("testFitTMMO_DBRExp.jl")
    end
    @testset "testFitTMMO_LayerLLL" begin
        include("testFitTMMO_LayerLLL.jl")
    end
    @testset "testFitTMMO_LayerLLLSpaceSolution" begin
        include("testFitTMMO_LayerLLLSpaceSolution.jl")
    end
    @testset "testFitTMMO_MC1dExp" begin
        include("testFitTMMO_MC1dExp.jl")
    end
    @testset "testFitBBO" begin
        include("testFitBBO.jl")
    end

    @testset "testFitCurves" begin
        include("testFitCurves.jl")
    end
end
