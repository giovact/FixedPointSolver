include("../src/FixedPointSolver.jl")
include("../src/visualization/Plotheaders.jl")
using DelimitedFiles, ArgParse, JLD, InteractiveUtils

parse_settings = ArgParseSettings()

@add_arg_table! parse_settings begin
####### Model and fixed parameters #######
    "--filename"
        arg_type = String
        required = true
    "--dir"
        arg_type = String
        default = "res/prove"
end

p_args = parse_args(parse_settings)



println("NOW PLOTTING")
println("Remember, give dir without the final /  and also give the filename without .jld extension")
plot_res_exploration(string(p_args["dir"], "/",p_args["filename"]))


