include("../src/FixedPointSolver.jl")
include("test_utils.jl")
using InteractiveUtils, Test, ArgParse

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
    "--models"
		nargs = '*'
        arg_type=String
        default = ["BinaryPerceptronRS";"BinaryPerceptron1RSB";"CoupledBinaryPerceptrons"; "BinaryPerceptronYInfNoDelta";"BinaryPerceptronYInf"]
end

p_args = parse_args(parse_settings)
models_totest = p_args["models"]

X = LegendreQuadrature(200)

# dummy constructor for testing
BinaryPerceptronYInf(rule::PT) where PT <:Potential = BinaryPerceptronYInf(1.0,1.0,1.0,X.n,rule,ones(3)*1e-2)

@testset verbose = true "Perceptrons with Parametric Potentials" begin
    for TM in subtypes(Perceptron)
        if string(TM) in models_totest
            for PT in subtypes(Potential)
                @testset verbose = true "$(string(TM)) with Potential $(string(PT))" begin
                    conv, iter = test_model(TM(PT()), X)
                    @test ( conv==1 && iter >1) 
                end
            end
        end
    end
end