struct RFIM <: FPModel
    params::NamedTuple{(:β,:h0,:σ,:J),Tuple{FT,FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RFIM() = RFIM(.1,0.0,.1,1.0,ones(1)) # dummy constructor for testing
RFIM(β,h0,σ,J) = RFIM(β,h0,σ,J, zeros(1)) 
RFIM(β,h0,σ,J,O) = RFIM((β=β,h0=h0,σ=σ,J=J),NamedArray(O,[:M]), Dict{String, FT}())

create_model(Model::RFIM, params) = RFIM(params, Model.O,Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
χI(Model::RFIM,z) = Model.params[:β]*(Model.params[:J]*Model.O[:M] + Model.params[:h0] + sqrt(Model.params[:σ])*z)

# order params SP
function exprM(Model::RFIM,X::TI) where TI <: IntegrationMethod
    IM(z) = tanh( χI(Model,z) )
    return integrate(X,IM)
end

function lhs!(Model::RFIM,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew[:M] = exprM(Model,X)
end

########################################################## FREE ENERGY #########################################################
function free_energy(Model::RFIM,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    energy = 0.5 * params[:J]*O[:M]^2

    IS(z) = log2cosh( χI(Model,z) )
    return  energy - integrate(X,IS)/params[:β]
end


# compute overlap too
function exprQ(Model::RFIM,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    IQ(z) = tanh( χI(Model,z))^2
    return integrate(X,IQ)
end
function compute_aux!(Model::RFIM,X::TI)  where TI <: IntegrationMethod
    Model.aux["Q"] = exprQ(Model, X)
end
