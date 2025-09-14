include("../src/FixedPointSolver.jl")
include("test_utils.jl")
using InteractiveUtils, Test, ArgParse

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
    "--models"
		nargs = '*'
        arg_type=String
        default = string.(models_retrieval)
end

p_args = parse_args(parse_settings)
models_totest = p_args["models"]

X = LegendreQuadrature(1000)

@testset verbose = true "Models with parametric type for DiscreteExpectation" begin
    for TM in models_retrieval
        if string(TM) in models_totest
            @testset verbose = true "$(string(TM)) with parametric type for DiscreteExpectation" begin
                @testset verbose = true "$(string(TM)) with Expectation of type SingleVar" begin
                    conv, iter = test_model(TM(DummySingleVar()), X)
                    @test ( conv==1 && iter >1) 
                end
                @testset verbose = true "$(string(TM)) with Expectation of type SumVar" begin
                    conv, iter = test_model(TM(DummySumVar()), X)
                    @test ( conv==1 && iter >1) 
                end
                @testset verbose = true "$(string(TM)) with Expectation of type VectorVar" begin
                    conv, iter = test_model(TM(DummyVectorVar()), X)
                    @test ( conv==1 && iter >1) 
                end
            end
        end
    end
end
