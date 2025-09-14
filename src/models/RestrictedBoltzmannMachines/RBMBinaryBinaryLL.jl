struct RBMBinaryBinaryLL{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:β,:αh),Tuple{FT,FT}}
    avg::RT
    Signal::RetrievalRBMBinaryBinary
    O::NamedVec #vector of orderparameters -> M, Q
    aux::Dict{String,FT} # other stuff, among them free energy
end


# constructors
RBMBinaryBinaryLL(avg::RT) where RT <:DiscreteExpectation = RBMBinaryBinaryLL(1.2,1.0,avg,ones(n_magnetizations(RetrievalRBMBinaryBinary(),avg))) # dummy constructor for testing
RBMBinaryBinaryLL(β,αh,avg,O) = RBMBinaryBinaryLL((β=β,αh=αh),avg,RetrievalRBMBinaryBinary(),setSignal_OP(O,avg), Dict{String,FT}())

create_model(Model::RBMBinaryBinaryLL, params) = RBMBinaryBinaryLL(params, Model.avg,RetrievalRBMBinaryBinary(),Model.O, Dict{String,FT}())

function sanity_checks!(Model::RBMBinaryBinaryLL)
    @assert n_order_params(Model) == n_magnetizations(Model.Signal,Model.avg)
end

function enforce_initial_condition!(Model::RBMBinaryBinaryLL{VectorVar})
    for μ=1:Model.avg.L 
        Model.O[Mμ(μ)] = rand()*(1-Model.avg.Mmin) + Model.avg.Mmin
    end
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
# order params SP

χI(Model::RBMBinaryBinaryLL,j::Int) = Model.params[:β] * Model.params[:αh] * Model.avg.store[j] 
exprM(Model::RBMBinaryBinaryLL,j::Int) = tanh(χI(Model,j))

function lhs!(Model::RBMBinaryBinaryLL{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : avg
    fill_avg!(Model.Signal,Model.avg,Model.O, Model.params[:β])
    @assert iszero(Onew)
    @assert length(unique(avg.probs))==1
    @assert avg.values[1] == - avg.values[2]
    
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,j)
    end

end

function lhs!(Model::RBMBinaryBinaryLL{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O avg
    @assert iszero(Onew)
    fill_avg!(Model.Signal,Model.avg,Model.O, Model.params[:β])
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,j)
    end
end

function lhs!(Model::RBMBinaryBinaryLL{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O avg
    @assert iszero(Onew)
    fill_avg!(Model.Signal,Model.avg,Model.O, Model.params[:β])
    for j=1:avg.Nstates
        for μ = 1:avg.L
            Onew[Mμ(μ)] += avg.probs[j] * avg.values[j,μ] * exprM(Model,j)
        end
    end
end

########################################################## FREE ENERGY #########################################################

exprEntropy(Model::RBMBinaryBinaryLL,j::Int) = log2cosh(χI(Model,j))

function free_energy(Model::RBMBinaryBinaryLL,X::TI) where TI <: IntegrationMethod
    @extract Model : params avg

    entropy = 0.0
    for j=1:avg.Nstates
        entropy += avg.probs[j] * exprEntropy(Model,j)
    end
    return  energySignal(Model.Signal,Model.avg,Model.O)  - entropy / params[:β]
end

_isMixedState(Model::RBMBinaryBinaryLL{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-4
_isMixedState(Model::RBMBinaryBinaryLL{VectorVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = all(Model.O[1:end-1] .> 1e-4)


function hessian_eigenvalue_final(Model::RBMBinaryBinaryLL{SumVar})
    1.0
end


function exprQ(Model::RBMBinaryBinaryLL{SumVar})
    1.0
end



function exprR(Model::RBMBinaryBinaryLL{SumVar})
   1.0
end


_isMixtureStable_th(Model::RBMBinaryBinaryLL{SumVar},X::TI,Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = hessian_eigenvalue_final(Model) > 0

function compute_aux!(Model::RBMBinaryBinaryLL{SumVar},X::TI) where TI <: IntegrationMethod
    @extract Model : avg params aux
    aux["is_SYM_MixtureStable"] = _isMixtureStable_th(Model, X, deepcopy(Model.O), 0.0)
    aux["hessian_final"] = hessian_eigenvalue_final(Model)
    aux["Q"] = exprQ(Model)
    aux["R"] = exprR(Model)
    aux["λ1"] = 1-params[:β]*(1-aux["Q"] + aux["R"])
    aux["λ2"] = 1-params[:β]*(1-aux["Q"]) + params[:β]*(avg.L-1)*aux["R"]
end
