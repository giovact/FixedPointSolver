include("../src/FixedPointSolver.jl")
using InteractiveUtils, Test, ArgParse, BenchmarkTools, Printf

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
    "--Model"
        arg_type = String
        required = true
    "--param"
        nargs = 2 # first argument name, second argument value, type inferred from the model type
        action = :append_arg
    "--parametric_type"
        nargs = '*'
        default = ["";""]
        
####### Model based Parameters #######
    "--IntMethod"
        arg_type = String
        default = "LegendreQuadrature"
    "--npointsInt"
        arg_type = Int
        default = 2000
    "--bound"
        arg_type = Float64
        default = 6.0
    "--tol"
        arg_type = Float64
        default = 1e-5
    "--maxiter"
        arg_type = Int
        default = 5000
    "--damp"
        arg_type = Float64
        default = 0.99
end

p_args = parse_args(parse_settings)

X = parse_integration_method(p_args)

ModelType, cparams_names, cparams_types = parse_Model_type(p_args["Model"])

Model = set_startingModel(ModelType,p_args["parametric_type"])

println("testing Model :", string(Model), " with Integration ", print_integration_params(X))
        
# just one iteration to compile
conv, iter, vareps = solve_SPequations!(Model, X; tol = p_args["tol"], K = FixedPoint(1.0), niter = 2, showinfo = false, printprogress = false)

teps = @elapsed conv, iter, vareps = solve_SPequations!(Model, X; tol = p_args["tol"], K = FixedPoint(p_args["damp"]), niter = p_args["maxiter"], showinfo = false, printprogress = false)
println("Converged in $iter")

str_teps = @sprintf "%.3f " teps
str_teps_perit = @sprintf "%.2e " teps / iter

println("Runtime in total ", str_teps, " seconds ***** # iter =  ",iter, "   time / iteration = ", str_teps_perit, "  seconds")
println(Model.O)


using BenchmarkTools
Model = set_startingModel(ModelType,p_args["parametric_type"]) 
@btime conv, iter, vareps = solve_SPequations!(Model, X; tol = p_args["tol"], K = FixedPoint(1.0), niter = 2, showinfo = false, printprogress = false)



println("Converged in $iter")