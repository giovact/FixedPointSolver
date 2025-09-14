include("../src/FixedPointSolver.jl")
include("test_utils.jl")
using InteractiveUtils, Test, ArgParse

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
    "--models"
		nargs = '*'
        arg_type=String
	"--history"
		arg_type = Bool
		default = true
end

p_args = parse_args(parse_settings)
models_totest = p_args["models"]
test_history = p_args["history"]

if length(models_totest)>0
	test_history=false
end

function test_model(Model::TM,X::TI) where {TM<: FPModel, TI <: IntegrationMethod}
    conv, iter, vareps = solve_SPequations!(Model, X; tol = 1e-5, niter = 10000,showinfo = false,printprogress = false)
	return conv, iter
end

X = UniformInterp(2048)
X = HermiteQuadrature(200)
X = LegendreQuadrature(200)
#X = HAdapt()


# missing constructors
CoupledHopfieldTSYInf() = CoupledHopfieldTSYInf(1.0,0.1,0.8,1.0,0,X.n, [0.1; 0.1; 0.1])


@testset verbose = true "Model convergence in simple settings" begin
	for TM in subtypes(FPModel)  
		if TM ∉ parametric_models ∪ [Perceptron]
			s = string(TM)
			if length(models_totest)==0 || s in models_totest
				@testset "Model $s" begin
					conv, iter = test_model(TM(), X)
					@test ( conv==1 && iter >1) 
				end
			end
		end
	end
end

if test_history
	@testset "History" begin
		Model = CurieWeiss()
		conv, vareps, history = solve_SPequations!(Model, NoInt(); showinfo = false, printprogress = false,savehistory = true)
		@test conv == 1

		Model = SK()
		conv, vareps, history = solve_SPequations!(Model,X; showinfo = false, printprogress = false,savehistory = true)
		@test conv == 1
	end
end


# add scanning test