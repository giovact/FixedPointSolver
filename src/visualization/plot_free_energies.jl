using PlotlyJS
function plot_free_energy(Model::TM,X::TI;fix_op = Dict()) where {TM <: FPModel, TI<:IntegrationMethod}
    nOP = n_order_params(Model)
    nOP_reduced = nOP - length(keys(fix_op))
    
    if nOP_reduced == 0
        error("there is no order parameter to vary")
    elseif nOP_reduced == 1
        plot_free_energy_1d(Model,X; fix_op = fix_op)
    elseif nOP_reduced == 2
        plot_free_energy_2d(Model,X; fix_op = fix_op)
    else
        throw("still to do")
    end
end

function plot_free_energy_1d(Model::TM,X::TI;fix_op = Dict(), bounds = [0.0; 1.0]) where {TM <: FPModel, TI<:IntegrationMethod}
    list_OP = filter( x -> x ∉ collect(keys(fix_op)), names(Model.O,1) )
    println(list_OP)
    @assert length(list_OP) == 1
    for Os in keys(fix_op)
        println("here")
        Model.O[Os] = fix_op[Os]
    end
    Ovec = collect(range(bounds[1], stop = bounds[2], length = 1000))
    fvec = zero(Ovec)
    for (iO,O) in enumerate(Ovec)
        Model.O[list_OP[1]] = O
        fvec[iO] = free_energy(Model, X)
    end
    PyPlot.plot(Ovec, fvec, label = string(params_to_textstring(Model;separator = " ", equalsymbol = "=")) ) 
    PyPlot.xlabel(string(list_OP[1]))
    PyPlot.ylabel("free energy")
    PyPlot.legend()
    PyPlot.title(string(typeof(Model)))
end

# Helix equation

function plot_free_energy_2d(Model::TM,X::TI;fix_op = Dict()) where {TM <: FPModel, TI<:IntegrationMethod}
    list_OP = filter( x -> x ∉ collect(keys(fix_op)), names(Model.O,1) )
    @assert length(list_OP) == 2
    for Os in keys(fix_op)
        Model.O[Os] = fix_op[Os]
    end
    O1vec = collect(range(0., stop = 1.0, length = 100))
    O2vec = collect(range(0.0, stop = 1.0, length = 100))
    fvec = zeros(length(O1vec),length(O2vec))
    @showprogress for (iO1,O1) in enumerate(O1vec)
        for (iO2,O2) in enumerate(O2vec)
            Model.O[list_OP[1]] = O1
            Model.O[list_OP[2]] = O2
            fvec[iO1,iO2] = free_energy(Model, X)
        end
    end
    println("OK")
    trace = PlotlyJS.surface(x=O1vec, y = O2vec, z = log.(abs.(fvec)), contours_z=attr(show=true,usecolormap=true,project_z=true))
    PlotlyJS.plot(trace)
    
    #PyPlot.plot(Ovec, fvec, label = string(params_to_textstring(Model;separator = " ", equalsymbol = "=")) ) 
    #PyPlot.xlabel(string(list_OP[1]))
    #PyPlot.ylabel("free energy")
    #PyPlot.legend()
    #PyPlot.title(string(typeof(Model)))

    #return fvec
end

