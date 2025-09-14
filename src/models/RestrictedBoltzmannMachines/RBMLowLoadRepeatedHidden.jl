struct RBMLowLoadRepeatedHidden <: FPModel
    params::NamedTuple{(:β,:αh),Tuple{FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RBMLowLoadRepeatedHidden() = RBMLowLoadRepeatedHidden(.1,1.0,ones(1)) # dummy constructor for testing
RBMLowLoadRepeatedHidden(β,αh) = RBMLowLoadRepeatedHidden(β,αh,zeros(1))
RBMLowLoadRepeatedHidden(β,αh,O) = RBMLowLoadRepeatedHidden((β=β,αh=αh),NamedArray(O,[:M]), Dict("f"=>0.0))

create_model(Model::RBMLowLoadRepeatedHidden, params) = RBMLowLoadRepeatedHidden(params, Model.O,Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function lhs!(Model::RBMLowLoadRepeatedHidden,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O params
    Onew[:M] = tanh( params[:β] * params[:αh] * tanh( params[:β] * O[:M]) )
end

########################################################## FREE ENERGY #########################################################
function free_energy(Model::RBMLowLoadRepeatedHidden,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    return 1.0
end
