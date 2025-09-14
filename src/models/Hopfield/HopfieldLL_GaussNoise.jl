struct HopLLGaussNoise <: FPModel
    params::NamedTuple{(:β,:σ),Tuple{FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopLLGaussNoise() = HopLLGaussNoise(.1,1.0,ones(2)) # dummy constructor for testing
HopLLGaussNoise(β,σ) = HopLLGaussNoise(β,σ,zeros(2))
HopLLGaussNoise(β,σ,O) = HopLLGaussNoise((β=β,σ=σ),NamedArray(O,[:M;:Q]), Dict("f"=>0.0))

create_model(Model::HopLLGaussNoise, params) = HopLLGaussNoise(params, Model.O,Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

# order params SP
function exprM(Model::HopLLGaussNoise,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    O[:M]==0.0 && return 0.0
    IM(z) = tanh( params[:β] * ( params[:σ] * sqrt(O[:Q]) * z + O[:M]) ) 
    return integrate(X,IM)
end

function exprQ(Model::HopLLGaussNoise,X::TI) where TI <: IntegrationMethod
    @extract Model : params O
    IQ(z) = tanh( params[:β] * ( params[:σ] * sqrt(O[:Q]) * z + O[:M]) )^ 2
    return integrate(X,IQ)
end

function lhs!(Model::HopLLGaussNoise,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew[:M] = exprM(Model,X)
    Onew[:Q] = exprQ(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::HopLLGaussNoise,X::TI) where TI <: IntegrationMethod
    @extract Model : params O

    energy = 0.5 * O[:M]^2 - 0.25 * params[:β] * params[:σ]^2 * (1 - O[:Q])^2

    IS(z) =  log( 2 * cosh( params[:β] * ( params[:σ] * sqrt(O[:Q]) * z + O[:M]) ) )
    return  energy - integrate(X,IS)/params[:β]
end
