struct BAMHLT0 <: FPModel
    params::NamedTuple{(:α,:γ),Tuple{FT,FT}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BAMHLT0() = BAMHLT0(0.1, 1.0, ones(2)) # dummy constructor for testing

BAMHLT0(α,γ) = BAMHLT0((α=α,γ=γ),NamedArray(zeros(2),[:y;:yb]),NamedArray(zeros(2),[:χ;:χb]), Dict{String,FT}())
BAMHLT0(α,γ,O) = BAMHLT0((α=α,γ=γ),NamedArray(O,[:y;:yb]), NamedArray(zeros(2),[:χ;:χb]), Dict{String,FT}())

create_model(Model::BAMHLT0, params) = BAMHLT0(params, Model.O, Model.Oconj, copy(Model.aux)) # maybe 

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
function exprOconj!(Model::BAMHLT0)
    # compute chi, chi bar and determinant 
    @extract Model : params O Oconj
    Oconj[:χ] =  (2 * params[:γ] / sqrt(π)) * (O[:y]*exp(-O[:y]^2)) / erf(O[:yb]) 
    Oconj[:χb] = (2 / (params[:γ] * sqrt(π))) * (O[:yb]*exp(-O[:yb]^2)) / erf(O[:y])
    Model.aux["Δ"] = 1 - Oconj[:χ]*Oconj[:χb]
end

# order params SP
function expry(Model::BAMHLT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:yb]==0.0 && return 0.0
    return ( erf(O[:yb]) * Model.aux["Δ"] ) / sqrt( (2*params[:α]*params[:γ]) * (1 + Oconj[:χb]^2) )
end

function expryb(Model::BAMHLT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:y]==0.0 && return 0.0
    return ( erf(O[:y]) * Model.aux["Δ"] ) / sqrt( (2*params[:α]/params[:γ]) * (1 + Oconj[:χ]^2) )
end

function lhs!(Model::BAMHLT0,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    exprOconj!(Model)
    Onew[:y] = expry(Model,X)
    Onew[:yb] = expryb(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::BAMHLT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O 
    return -12345.0
end

function compute_aux!(Model::BAMHLT0,X::TI;setNaNs::Bool = false) where TI <: IntegrationMethod
    @extract Model : params O 
    Model.aux["M"] = erf(O[:y])
    Model.aux["Mb"] = erf(O[:yb])
end
