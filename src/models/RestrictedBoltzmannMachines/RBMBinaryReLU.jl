struct RBMBinaryReLU{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:β,:α,:p,:η,:g,:λ,:Qhlog),Tuple{FT,FT,FT,FT,FT,FT,Bool}}
    Signal::RT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RBMBinaryReLU(Signal::RT) where RT <: DiscreteExpectation = RBMBinaryReLU(.1, .1,1.0,0.0,0.0,1.0,true,Signal,ones(Signal.L+1)) # dummy constructor for testing
RBMBinaryReLU(β,α,p,η,g,λ,Qhlog,Signal) = RBMBinaryReLU(β,α,p,η,g,λ,Qhlog,Signal,zeros(1+Signal.L))
RBMBinaryReLU(β,α,p,η,g,λ,Qhlog,Signal,O) = RBMBinaryReLU((β=β,α=α,p=p,η=η,g=g,λ=λ,Qhlog=Qhlog),Signal,vcat(setSignal_OP(O,Signal),NamedArray(O[end:end],[:Q])),vcat(setSignal_OP(zeros(Signal.L+1),Signal),NamedArray(zeros(1),[:Q])), Dict{String,FT}())

create_model(Model::RBMBinaryReLU, params) = RBMBinaryReLU(params, Model.Signal,Model.O, Model.Oconj, Dict{String,FT}())

sanity_check_nop!(Model::RBMBinaryReLU{SingleVar}) = @assert n_order_params(Model) == 2
sanity_check_nop!(Model::RBMBinaryReLU{SumVar}) = @assert n_order_params(Model) == 2
sanity_check_nop!(Model::RBMBinaryReLU{VectorVar}) = @assert n_order_params(Model) == Model.Signal.L + 1

function sanity_checks!(Model::RBMBinaryReLU)
    @assert is_free_energy_definite(Model)
    if Model.params[:p] != 1.0
        @assert Model.Signal.ξ0 == [-1;0;1]
        @assert  Model.Signal.p0[1] == 0.5*Model.params[:p]
        @assert  Model.Signal.p0[2] == 1-Model.params[:p]
        @assert  Model.Signal.p0[3] == 0.5*Model.params[:p]
    end
    sanity_check_nop!(Model)
end

function enforce_initial_condition!(Model::RBMBinaryReLU{VectorVar})
    for μ=1:Model.Signal.L 
        Model.O[Mμ(μ)] = rand()*(1-Model.Signal.Mmin) + Model.Signal.Mmin
    end
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::RBMBinaryReLU) =  Model.params[:λ] - Model.params[:β] * Model.params[:p] * (1 - Model.O[:Q])
Deltas!(Model::RBMBinaryReLU) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::RBMBinaryReLU) = Delta(Model)>0

#### Conjugate order params SP
function aux_variables!(Model::RBMBinaryReLU)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    aux["a"] = sqrt(params[:β] * params[:p] * O[:Q] / aux["Δ"])
    aux["b"] = - sqrt(params[:β] / aux["Δ"]) * params[:η]
    aux["Δplus"] =  Model.params[:λ] - Model.params[:β] * Model.params[:p] * (1 - 2 * Model.O[:Q])
    aux["Δminus"] =  Model.params[:λ] - Model.params[:β] * Model.params[:p]
end

exprMhat!(Model::RBMBinaryReLU{SingleVar}) = Model.Oconj[:M] = max(0.0, Model.O[:M] / Model.params[:λ])
exprMhat!(Model::RBMBinaryReLU{SumVar}) = Model.Oconj[:M] = max(0.0, Model.O[:M] / Model.params[:λ]) 

function exprMhat!(Model::RBMBinaryReLU{VectorVar})
    for μ = 1:Model.Signal.L
        Model.Oconj[Mμ(μ)] = max(0.0, Model.O[Mμ(μ)] / Model.params[:λ])
    end
end

function prefactors_Qhat(Model::RBMBinaryReLU)
    @extract Model : O Oconj params aux
    first_piece = ( params[:p] / aux["Δ"]^2 ) * ( params[:p]*O[:Q] + params[:η]^2 )

    prefactor_M1 = sqrt(2/ ( π*params[:β]*aux["Δ"]) ) * (params[:p]*params[:η]) / (aux["Δplus"]) 
    prefactor_M2 = 0.5 * params[:p] / ( π * params[:β] * aux["Δ"]) * (aux["Δminus"] / aux["Δplus"])

    aux["first_piece"] = first_piece
    aux["prefactor_M1"] = prefactor_M1
    aux["prefactor_M2"] = prefactor_M2
end


Herf_argument(Model::RBMBinaryReLU, z) = Model.aux["a"]*z + Model.aux["b"]

function exprQhat(Model::RBMBinaryReLU,X::TI) where TI <: IntegrationMethod
    @extract Model : O Oconj params aux
    aux_variables!(Model)
    prefactors_Qhat(Model)
    M1_int(z) = Herf_argument(Model,z) < 38.6 ? exp(-0.5 * Herf_argument(Model,z)^2 ) / Herf(Herf_argument(Model,z)) : Herf_argument(Model,z) * sqrt(2π)
    M2_int(z) = M1_int(z)^2
    out =  aux["first_piece"]  + aux["prefactor_M2"] * integrate(X,M2_int)
    if params[:η] !=0
        out += aux["prefactor_M1"] * integrate(X,M1_int)
    end
    # if isnan(out)
    #     try a stupid expansion when \eta goes to - infty
    #     out = first_piece - (2 * params[:p] * params[:η]^2) / (aux["Δplus"] * aux["Δ"]) + ( params[:p] * aux["Δminus"] * params[:η]^2 ) / (aux["Δplus"] * aux["Δ"]^2 )
    
    return out
end


function exprQhat_logint(Model::RBMBinaryReLU,X::TI) where TI <: IntegrationMethod
    @extract Model : O Oconj params aux
    aux_variables!(Model)
    prefactors_Qhat(Model)
    
    logM1_int(z) = - 0.5 * Herf_argument(Model, z)^2 - logHerf(Herf_argument(Model,z))
    logM2_int(z) = 2 * logM1_int(z)
    out =  aux["first_piece"] + aux["prefactor_M2"] * integrate_log(X,logM2_int)
    if params[:η] !=0
        out += aux["prefactor_M1"] * integrate_log(X,logM1_int)
    end
    return out
end


function exprOconj!(Model::RBMBinaryReLU,X::TI) where TI <: IntegrationMethod
    exprMhat!(Model) # that depends on the signal term
    lb = 0.0
    if Model.params[:Qhlog]
        Model.Oconj[:Q] = clamp(exprQhat_logint(Model,X),lb,Inf) 
    else
        Model.Oconj[:Q] = clamp(exprQhat(Model,X),lb,Inf)
    end
end

# order params SP
# in these functions j is the index of the component over all the possible configurations of the pattern to be averaged over
Signaldot(Model::RBMBinaryReLU{SingleVar},j::Int) = Model.Oconj[:M] * Model.Signal.values[j]
Signaldot(Model::RBMBinaryReLU{SumVar},j::Int) = Model.Oconj[:M] * Model.Signal.values[j]
Signaldot(Model::RBMBinaryReLU{VectorVar},j::Int) = sum(Model.Oconj[Mμ(μ)] * Model.Signal.values[j,μ] for μ=1:Model.Signal.L)

χI(Model::RBMBinaryReLU,z,j::Int) = Model.params[:β]* (Model.params[:g] +  Signaldot(Model,j)  + sqrt(Model.params[:α]*Model.Oconj[:Q]) * z) 

function exprM(Model::RBMBinaryReLU,X::TI,j::Int) where TI <: IntegrationMethod
    IM(z) = tanh( χI(Model,z,j) )
    integrate(X,IM)
end

function exprQ(Model::RBMBinaryReLU,X::TI,j::Int) where TI <: IntegrationMethod
    IQ(z) = tanh( χI(Model,z,j) )^2
    integrate(X,IQ)
end

function lhs!(Model::RBMBinaryReLU{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model,X)
    @assert iszero(Onew)

    for j=1:Signal.Nstates
        Onew[:M] += Signal.probs[j] * Signal.values[j] * exprM(Model,X,j)
        Onew[:Q] += Signal.probs[j] * exprQ(Model,X,j)
    end
    #@assert Onew[:M] in Interval(-1.0,1.0)
end

function lhs!(Model::RBMBinaryReLU{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model,X)
    @assert iszero(Onew)
    
    for j=1:Signal.Nstates
        Onew[:M] += Signal.probs[j] * (Signal.values[j] / Signal.L) * exprM(Model,X,j) 
        Onew[:Q] += Signal.probs[j] * exprQ(Model,X,j)
    end
end

function lhs!(Model::RBMBinaryReLU{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model,X)
    @assert iszero(Onew)
    
    for j=1:Signal.Nstates
        for μ = 1:Signal.L
            Onew[Mμ(μ)] += Signal.probs[j] * Signal.values[j,μ] * exprM(Model,X,j)
        end
        Onew[:Q] += Signal.probs[j] * exprQ(Model,X,j)
    end
end

########################################################## FREE ENERGY #########################################################

function energetic_term_IE(Model::RBMBinaryReLU,X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O

    IE(z) = logHerf(aux["a"]*z+aux["b"])

    return 0.5 * log(2π) - 0.5 * log(params[:β]) - 0.5 * log(aux["Δ"])  + 0.5*(aux["a"]^2+aux["b"]^2) + integrate(X,IE)
end

function exprEntropy(Model::RBMBinaryReLU,X::TI,j::Int) where TI <: IntegrationMethod
    IQ(z) = log2cosh( χI(Model,z,j) )
    integrate(X,IQ)
end

energySignal(Model::RBMBinaryReLU{SingleVar}) = Model.O[:M] * Model.Oconj[:M] - 0.5*(Model.O[:M]^2 / Model.params[:λ]) * (Model.O[:M]>0.0) 
energySignal(Model::RBMBinaryReLU{SumVar}) = Model.Signal.L * Model.O[:M] * Model.Oconj[:M] - 0.5*Model.Signal.L *(Model.O[:M]^2 / Model.params[:λ]) * (Model.O[:M]>0.0) 
energySignal(Model::RBMBinaryReLU{VectorVar}) = sum(Model.O[Mμ(μ)] * Model.Oconj[Mμ(μ)] - 0.5*(Model.O[Mμ(μ)]^2 / Model.params[:λ]) * (Model.O[Mμ(μ)]>0.0)  for μ =1:Model.Signal.L )


function free_energy(Model::RBMBinaryReLU,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux Signal
    exprOconj!(Model,X)

    energy_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:Q]

    logdetK = (params[:α]/(2*params[:β])) * log(aux["Δ"])  
    other =  0.5*params[:α]* (params[:p] * O[:Q] + params[:η]^2 ) / aux["Δ"] 

    IE(z) = logHerf(aux["a"]*z + aux["b"])

    entropy = 0.0
    for j=1:Signal.Nstates
        entropy += Signal.probs[j] * exprEntropy(Model,X,j)
    end

    return  energySignal(Model) + energy_overlaps + logdetK - other - (params[:α]/params[:β])*integrate(X,IE) - entropy/params[:β]
end



_isRetrieval(Model::RBMBinaryReLU{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isNOTRetrieval(Model::RBMBinaryReLU{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] < 1e-4
_isRetrieval(Model::RBMBinaryReLU{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-3


function _isRetrievalStateDominant(Model::RBMBinaryReLU{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod
    ModelSG = create_model(deepcopy(Model), Model.params)
    ModelSG.O[:Q] = 1.0 
    ModelSG.O[:M] = 0.0
    solve_SPequations!(ModelSG, X; tol = tol, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 10000,force_dict = Dict(:M=>0.0))
    return Model.aux["f"] < ModelSG.aux["f"]
end


function _isRetrievalStateDominant2(Model::RBMBinaryReLU{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod
    ModelSG = create_model(deepcopy(Model), Model.params)
    ModelSG.O[:Q] = 1.0 
    ModelSG.O[:M] = 1e-4
    solve_SPequations!(ModelSG, X; tol = tol, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 10000)
    if isapprox(Model.O[:M],ModelSG.O[:M],atol=1e-4)
        return true
    else
        return Model.aux["f"] < ModelSG.aux["f"]
    end
end

