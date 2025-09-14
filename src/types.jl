
const FT=Float64
const Vec = Vector{FT}
const NamedVec = NamedVector{FT, Vector{FT}, Tuple{OrderedDict{Symbol, Int64}}}
const NamedMat = NamedMatrix{FT, Matrix{FT}, Tuple{OrderedDict{Symbol, Int64}, OrderedDict{Symbol, Int64}}}

abstract type FPModel
end

# each model has at least the following fields (required for type-unspecific functions in solver.jl and other stuff)

# params : NamedTuples for control parameters
# O : Named Array for order parameters
# aux : dictionary for auxiliary stuff (among them, free energy)

# eventually also
# Oconj : Named Array for conjugate order parameters

function Base.string(Model::TM; separator::String="_") where {TM <: FPModel}
    out = string(typeof(Model))
    if occursin("{",out)
        out = replace(out, "{"=>separator)
    end
    if occursin("}",out)
        out = replace(out, "}"=>separator)
    end
    if hasfield(typeof(Model),:V)
        out = string(out, Model.V.storage ? "storage_" : "")
    end

    if hasfield(typeof(Model),:avg)
        out = string(out, "L",separator,Model.avg.L)
    end
    
    return out
end


dummy_model(T::Type{<:FPModel}) = T()

n_control_params(Model::TM) where {TM <: FPModel} = length(Model.params)
n_order_params(Model::TM) where {TM <: FPModel} = length(Model.O)
n_conj_order_params(Model::TM) where {TM <: FPModel} = hasfield(typeof(Model),:Oconj) ? length(Model.Oconj) : 0 
has_conj_order_params(Model::TM) where {TM <: FPModel} = hasfield(typeof(Model),:Oconj)

function printvalues(Model::TM,ε::FT,printaux::Vector{String}) where {TM <: FPModel} # TO DO, modify vector in place, declare in solve_SPequations!
    toprint = [ (n,Model.O[n]) for n in NamedArrays.names(Model.O,1) ]
    if has_conj_order_params(Model)
        append!(toprint, [(Symbol(n,:conj),Model.Oconj[n]) for n in NamedArrays.names(Model.Oconj,1)])
    end
    append!(toprint,[(:ε,ε)])
    if !isempty(printaux)
        append!(toprint, [(Symbol(aa),Model.aux[aa]) for aa in printaux])
    end
    return toprint
end

function print_op_withparam(Model::TM,param::Symbol) where {TM <: FPModel} 
    Onames = vcat([(param,Model.params[param])], [(n,Model.O[n]) for n in NamedArrays.names(Model.O,1) ])
    if has_conj_order_params(Model)
        append!(Onames, [(Symbol(n,:conj),Model.Oconj[n]) for n in NamedArrays.names(Model.Oconj,1)])
    end
    return Onames
end

function get_Oconj(Model::TM) where TM<:FPModel
    has_conj_order_params(Model) && return copy(Vector(Model.Oconj))
    return zeros(0)
end

function get_Oconjnames(Model::TM) where TM <: FPModel

    !has_conj_order_params(Model) && return nothing

    return [string(On, "conj") for On in names(Model.Oconj,1)]
end

function print_op(O::NamedVec; separator::String = " - ")
    out = ""; 
    for n in names(O,1)
        a = @sprintf "%.3f" O[n]
        out *= string(n,"=",a, separator)
    end
    out[1:end - length(separator)]
end
print_op(Model::TM; separator::String=" - ") where TM <: FPModel  = print_op(Model.O; separator=separator)


function Base.show(io::IO, Model::TM) where {TM <: FPModel}
    print(io, "Model : ", typeof(Model), " - params = (", params_to_textstring(Model; separator = " - ",equalsymbol="=")[1:end-3], ")  - OP --> ", print_op(Model) )
end