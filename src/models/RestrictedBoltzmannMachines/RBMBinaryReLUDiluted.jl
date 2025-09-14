struct RBMBinaryReLUDiluted <: FPModel
    params::NamedTuple{(:β,:α,:l,:η,:g,:λ),Tuple{FT,FT,FT,FT,FT,FT}}
    Signal::BesselVar
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RBMBinaryReLUDiluted() = RBMBinaryReLUDiluted(1.0,0.1,1.0,0.0,0.0,1.0,ones(2))
RBMBinaryReLUDiluted(β,α,l,η,g,λ) = RBMBinaryReLUDiluted(β,α,l,η,g,λ,zeros(2))
RBMBinaryReLUDiluted(β,α,l,η,g,λ,O) = RBMBinaryReLUDiluted((β=β,α=α,l=l,η=η,g=g,λ=λ),BesselVar(l),NamedArray(O,[:M;:Q]),NamedArray(zeros(2),[:M;:Q]), Dict{String,FT}())
create_model(Model::RBMBinaryReLUDiluted, params) = RBMBinaryReLUDiluted(params, Model.Signal,Model.O, Model.Oconj, Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::RBMBinaryReLUDiluted) =  Model.params[:λ] - Model.params[:β] * (1 - Model.O[:Q])
Deltas!(Model::RBMBinaryReLUDiluted) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::RBMBinaryReLUDiluted) = Delta(Model)>0

sanity_checks!(Model::RBMBinaryReLUDiluted) = @assert Model.params[:l] == Model.Signal.l

#### Conjugate order params SP
function aux_variables!(Model::RBMBinaryReLUDiluted)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    aux["a"] = sqrt(params[:β] *  O[:Q] / aux["Δ"])
    aux["b"] = - sqrt(params[:β] / aux["Δ"]) * params[:η]
    aux["Δplus"] =  Model.params[:λ] - Model.params[:β] * (1 - 2 * Model.O[:Q])
    aux["Δminus"] =  Model.params[:λ] - Model.params[:β]
end

function exprMhat!(Model::RBMBinaryReLUDiluted) 
    Model.Oconj[:M] = max(0.0, Model.O[:M] / Model.params[:λ])
end

function exprQhat(Model::RBMBinaryReLUDiluted,X::TI) where TI <: IntegrationMethod
    @extract Model : O Oconj params aux

    first_piece = ( O[:Q] + params[:η]^2 ) / aux["Δ"]^2

    prefactor_M1 = sqrt(2/ ( π*params[:β]*aux["Δ"]) ) * (params[:η]) / (aux["Δplus"]) 
    M1_int(z) = exp(- 0.5 * ( aux["a"]*z + aux["b"])^2 ) / Herf(aux["a"]*z + aux["b"])

    prefactor_M2 = 0.5 * (aux["Δminus"] / aux["Δplus"]) /  ( π * params[:β] * aux["Δ"])
    M2_int(z) = exp(- (aux["a"] *z + aux["b"])^2 ) / Herf(aux["a"]*z + aux["b"])^2
    

    out =  first_piece + prefactor_M1 * integrate(X,M1_int) + prefactor_M2 * integrate(X,M2_int)
    
    return clamp(out, 1e-20, Inf)
end

function exprOconj!(Model::RBMBinaryReLUDiluted,X::TI) where TI <: IntegrationMethod
    aux_variables!(Model)
    exprMhat!(Model) # that depends on the signal term
    Model.Oconj[:Q] = exprQhat(Model,X)
end

# order params SP

χI(Model::RBMBinaryReLUDiluted,z,j::Int) = Model.params[:β]* (Model.params[:g] +  Model.O[:M] * Model.Signal.values[j]  + sqrt(Model.params[:α]*Model.Oconj[:Q]) * z) 

function exprM(Model::RBMBinaryReLUDiluted,X::TI,j::Int) where TI <: IntegrationMethod
    IM(z) = tanh( χI(Model,z,j) )
    integrate(X,IM)
end

function exprQ(Model::RBMBinaryReLUDiluted,X::TI,j::Int) where TI <: IntegrationMethod
    IQ(z) = tanh( χI(Model,z,j) )^2
    integrate(X,IQ)
end

function lhs!(Model::RBMBinaryReLUDiluted,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model,X)
    @assert iszero(Onew)

    for j=1:Signal.Nstates
        Onew[:M] += Signal.probs[j] * (Signal.values[j] / Model.params[:l]) * exprM(Model,X,j)
        Onew[:Q] += Signal.probs[j] * exprQ(Model,X,j)
    end
    #@assert Onew[:M] in Interval(-1.0,1.0)
end

########################################################## FREE ENERGY #########################################################
function exprEntropy(Model::RBMBinaryReLUDiluted,X::TI,j::Int) where TI <: IntegrationMethod
    IS(z) = log2cosh( χI(Model,z,j) )
    integrate(X,IS)
end
 
function free_energy(Model::RBMBinaryReLUDiluted,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux Signal
    exprOconj!(Model,X)

    energy_signal = params[:l] * O[:M] * Oconj[:M] - 0.5*params[:l] *(O[:M]^2 / params[:λ]) * (O[:M]>0.0) 

    energy_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:Q]

    logdetK = (params[:α]/(2*params[:β])) * log(aux["Δ"])  
    other =  0.5*params[:α]* ( O[:Q] + params[:η]^2 ) / aux["Δ"] 

    IE(z) = logHerf(aux["a"]*z + aux["b"])

    entropy = 0.0
    for j=1:Signal.Nstates
        entropy += Signal.probs[j] * exprEntropy(Model,X,j)
    end

    return  energy_signal + energy_overlaps + logdetK - other - (params[:α]/params[:β])*integrate(X,IE) - entropy/params[:β]
end
