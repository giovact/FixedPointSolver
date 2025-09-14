struct RBMBinaryBinaryLL_onlyPure <: FPModel
    params::NamedTuple{(:β,:αh),Tuple{FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RBMBinaryBinaryLL_onlyPure() = RBMBinaryBinaryLL_onlyPure(.1,1.0,ones(1)) # dummy constructor for testing
RBMBinaryBinaryLL_onlyPure(β,αh) = RBMBinaryBinaryLL_onlyPure(β,αh,zeros(1))
RBMBinaryBinaryLL_onlyPure(β,αh,O) = RBMBinaryBinaryLL_onlyPure((β=β,αh=αh),NamedArray(O,[:M]), Dict("f"=>0.0))

create_model(Model::RBMBinaryBinaryLL_onlyPure, params) = RBMBinaryBinaryLL_onlyPure(params, Model.O,Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

# order params SP
function exprM(Model::RBMBinaryBinaryLL_onlyPure,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    O[:M]==0.0 && return 0.0
    return tanh( params[:β] * params[:αh] * tanh( params[:β] * O[:M]) ) 
end


function lhs!(Model::RBMBinaryBinaryLL_onlyPure,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew[:M] = exprM(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::RBMBinaryBinaryLL_onlyPure,X::TI) where TI <: IntegrationMethod
    @extract Model : params O

    return 1.0
end
