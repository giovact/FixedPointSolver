struct SK <: FPModel
    params::NamedTuple{(:β,:h,:J0,:J),Tuple{FT,FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
SK() = SK(0.1,0.0,0.0,1.0,ones(2)) # dummy constructor for testing
SK(β,h,J0,J) = SK(β,h,J0,J, zeros(2))
SK(β,h,J0,J,O) = SK((β=β,h=h,J0=J0,J=J),NamedArray(O,[:M;:Q]), Dict{String,FT}())

create_model(Model::SK, params) = SK(params, Model.O,Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function exprM(Model::SK,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    O[:M]==0.0 && return 0.0
    IM(z) = tanh( params[:β] * ( params[:J] * sqrt(O[:Q]) * z + params[:J0]*O[:M] + params[:h] ) ) 
    return integrate(X,IM)
end

function exprQ(Model::SK,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    IQ(z) = tanh( params[:β] * ( params[:J] * sqrt(O[:Q]) * z + params[:J0]*O[:M] + params[:h] ) )^2
    return integrate(X,IQ)
end

function lhs!(Model::SK,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew[:M] = exprM(Model,X)
    Onew[:Q] = exprQ(Model,X)
end


########################################################## FREE ENERGY #########################################################

function free_energy(Model::SK,X::TI) where TI <: IntegrationMethod
    @extract Model : params O

    energy = 0.5 * params[:J0] * O[:M]^2 - 0.25 * params[:β] * params[:J]^2 * (1 - O[:Q])^2

    IS(z) =  log( 2 * cosh( params[:β] * ( params[:J] * sqrt(O[:Q]) * z + params[:J0]*O[:M] + params[:h] ) ) )
    return  energy - integrate(X,IS)/params[:β]
end
