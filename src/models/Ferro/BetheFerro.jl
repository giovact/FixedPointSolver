struct BetheFerro <: FPModel
    params::NamedTuple{(:β,:J,:h,:d),Tuple{FT,FT,FT,Int}}
    O::NamedVec #vector of orderparameters -> M #
    aux::Dict{String,FT} # other stuff, among them free energy
end

#constructors
BetheFerro() = BetheFerro(2.0,0.0,2,ones(1)) # dummy constructor for testing
BetheFerro(β,h,d,H0;J=1.0) = BetheFerro((β=β,J=J,h=h,d=d), NamedArray(H0,[:Heff]),Dict{String,FT}())
BetheFerro(β,h,d;J=1.0) = BetheFerro(β,h,d,zeros(1);J=J) 

create_model(Model::BetheFerro, params) = BetheFerro(params, Model.O, Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function exprHeff(Model::BetheFerro)
    @extract Model : params O
    return params[:β]*params[:h] + (2*params[:d] - 1) * atanh( tanh(params[:β]*params[:J])*tanh(O[:Heff])) 
end

function lhs!(Model::BetheFerro,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    Onew[:Heff] = exprHeff(Model)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::BetheFerro,X::TI) where TI <: IntegrationMethod
    return -1.0 # you should write the expression of the Bethe free energy
end

function compute_aux!(Model::BetheFerro,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    Model.aux["M"] = tanh(params[:β]*params[:h] + 2*params[:d]*atanh( tanh(params[:β]*params[:J])*tanh(O[:Heff])) )
end
    

