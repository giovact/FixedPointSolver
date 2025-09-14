struct CurieWeiss <: FPModel
    params::NamedTuple{(:β,:J,:h),Tuple{FT,FT,FT}}
    O::NamedVec
    aux::Dict{String,FT} # other stuff, among them free energy
end

#constructors
CurieWeiss() = CurieWeiss(2.0,1.0,0.0,ones(1)) # dummy constructor for testing
CurieWeiss(β,J,h) =  CurieWeiss(β,J,h,zeros(1))
CurieWeiss(β,J,h,O) = CurieWeiss((β=β,J=J,h=h), NamedArray(O,[:M]),Dict{String,FT}())
create_model(Model::CurieWeiss, params) = CurieWeiss(params, Model.O, Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function exprM(Model::CurieWeiss,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    return tanh(params[:β]* ( params[:h] + params[:J]*O[:M]))
end

function lhs!(Model::CurieWeiss,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew[:M] = exprM(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::CurieWeiss,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    energy = - params[:h] * O[:M] - 0.5 * params[:J] * O[:M]^2
    entropy =  0.5 * ( 1 + O[:M] ) * log( 0.5 * ( 1 + O[:M] ) ) + 0.5 * ( 1 - O[:M] ) * log( 0.5 * ( 1 - O[:M] ) ) 
    return energy + entropy/(params[:β])
end

function compute_aux!(Model::CurieWeiss,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    Model.aux["saddlepointMAX"] = - 0.5 * params[:J] * O[:M]^2 + log2cosh(params[:β]* ( params[:h] + params[:J]*O[:M]))/params[:β]
end


function free_energy2(Model::CurieWeiss,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    return 0.5 * params[:J] * O[:M]^2 - log2cosh(params[:β]* ( params[:h] + params[:J]*O[:M]))/params[:β]
end


_isCWmetastable(Model::CurieWeiss, X::TI, Oold::NamedVec)  where TI <: IntegrationMethod = sign(Model.params[:h] * Model.O[:M]) < 0