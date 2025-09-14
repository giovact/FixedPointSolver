struct HopfieldLL{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:β,),Tuple{FT}}
    avg::RT
    Signal::RetrievalHopfield
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end


# constructors
HopfieldLL(avg::RT) where RT <:DiscreteExpectation = HopfieldLL(1.2,avg,ones(n_magnetizations(RetrievalHopfield(),avg))) # dummy constructor for testing
HopfieldLL(β,avg,O) = HopfieldLL((;β=β),avg,RetrievalHopfield(),setSignal_OP(O,avg), Dict{String,FT}())

create_model(Model::HopfieldLL, params) = HopfieldLL(params, Model.avg,RetrievalHopfield(),Model.O, Dict{String,FT}())

function sanity_checks!(Model::HopfieldLL)
    @assert n_order_params(Model) == n_magnetizations(Model.Signal,Model.avg)
end

function enforce_initial_condition!(Model::HopfieldLL{VectorVar})
    for μ=1:Model.avg.L 
        Model.O[Mμ(μ)] = rand()*(1-Model.avg.Mmin) + Model.avg.Mmin
    end
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
# order params SP

χI(Model::HopfieldLL,j::Int) = Model.params[:β] * Model.avg.store[j] 
exprM(Model::HopfieldLL,j::Int) = tanh(χI(Model,j))

function lhs!(Model::HopfieldLL{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : avg
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)
    @assert length(unique(avg.probs))==1
    @assert avg.values[1] == - avg.values[2]
    
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,j)
    end

end

function lhs!(Model::HopfieldLL{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O avg
    @assert iszero(Onew)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,j)
    end
end

function lhs!(Model::HopfieldLL{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O avg
    @assert iszero(Onew)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    for j=1:avg.Nstates
        for μ = 1:avg.L
            Onew[Mμ(μ)] += avg.probs[j] * avg.values[j,μ] * exprM(Model,j)
        end
    end
end

########################################################## FREE ENERGY #########################################################

exprEntropy(Model::HopfieldLL,j::Int) = log2cosh(χI(Model,j))

function free_energy(Model::HopfieldLL,X::TI) where TI <: IntegrationMethod
    @extract Model : params avg

    entropy = 0.0
    for j=1:avg.Nstates
        entropy += avg.probs[j] * exprEntropy(Model,j)
    end
    return  energySignal(Model.Signal,Model.avg,Model.O)  - entropy / params[:β]
end

_isMixedState(Model::HopfieldLL{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-4
_isMixedState(Model::HopfieldLL{VectorVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = all(Model.O[1:end-1] .> 1e-4)


function hessian_eigenvalue_final(Model::HopfieldLL{SumVar})
    @extract Model : params O
    
    if Model.avg.L == 1
        return 1 - params[:β] * ( 1 - tanh(params[:β]*O[:M])^2 )
    elseif iseven(Model.avg.L)
        return 1 - params[:β]
    else
        avg_red = MixedStatesRetrievalBinarySYM(Model.avg.L - 2)
        av_tanh2 = 0.0
        for j=1:avg_red.Nstates
            av_tanh2 += avg_red.probs[j] * tanh(params[:β] * O[:M] * avg_red.values[j]) ^2
        end
        return 1 - params[:β] * (1 - av_tanh2 )
    end
end


function exprQ(Model::HopfieldLL{SumVar})
    @extract Model : avg
    Q = 0.0
    for (j,S) in enumerate(avg.values)
        Q += avg.probs[j] * tanh(χI(Model,j))^2
    end
    return Q
end



function exprR(Model::HopfieldLL{SumVar})
    @extract Model : avg
    avg.L==1 && return 0.0
    R = 0.0
    for (j,S) in enumerate(avg.values)
        R += avg.probs[j] * (S^2-avg.L)*tanh(χI(Model,j))^2
    end
    return R / ( avg.L * (avg.L-1) ) 
end


_isMixtureStable_th(Model::HopfieldLL{SumVar},X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = hessian_eigenvalue_final(Model) > 0

function compute_aux!(Model::HopfieldLL{SumVar},X::TI) where TI <: IntegrationMethod
    @extract Model : avg params aux
    aux["is_SYM_MixtureStable"] = _isMixtureStable_th(Model, X, deepcopy(Model.O), 0.0)
    aux["hessian_final"] = hessian_eigenvalue_final(Model)
    aux["Q"] = exprQ(Model)
    aux["R"] = exprR(Model)
    aux["λ1"] = 1-params[:β]*(1-aux["Q"] + aux["R"])
    aux["λ2"] = 1-params[:β]*(1-aux["Q"]) + params[:β]*(avg.L-1)*aux["R"]
end
