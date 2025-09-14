struct BAMHLT0Mixtures{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:α,:γ),Tuple{FT,FT}}
    avg::RT
    Signal::RetrievalBAM
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors  
BAMHLT0Mixtures(avg::RT) where {RT<:DiscreteExpectation} = BAMHLT0Mixtures(0.1, 1.0, avg, 
                                                            vcat([ones(n_magnetizations_singlelayer(RetrievalBAM(avg),avg));0.1],[ones(n_magnetizations_singlelayer(RetrievalBAM(avg),avg));0.1] ) )
BAMHLT0Mixtures(α, γ, avg) = BAMHLT0Mixtures(α, γ, avg, zeros(2 + n_magnetizations_bothlayers(RetrievalBAM(avg),avg)))

function BAMHLT0Mixtures(α, γ, avg, O)
    Signal = RetrievalBAM(avg)
    @assert length(O) == 2 + n_magnetizations_bothlayers(Signal,avg)
    no_M = n_magnetizations_singlelayer(Signal,avg) 
    return BAMHLT0Mixtures((α=α, γ=γ), 
                    avg,
                    RetrievalBAM(avg),
                    vcat([setBAMsignal_L1(O[1:no_M], avg);NamedArray([O[no_M+1]],[:χ])],[setBAMsignal_L2(O[no_M+2:1+2*no_M], avg);NamedArray([O[2*no_M+2]],[:χb])]), 
                    NamedArray(zeros(2),[:P;:Pb]),
                    Dict{String,FT}())
end

create_model(Model::BAMHLT0Mixtures, params) = BAMHLT0Mixtures(params, Model.avg, Model.Signal,Model.O, Model.Oconj,Dict{String,FT}())

function sanity_checks!(Model::BAMHLT0Mixtures)
    @assert n_order_params(Model) ==  2 + n_magnetizations_bothlayers(Model.Signal,Model.avg)
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

Delta(Model::BAMHLT0Mixtures) = 1 - Model.O[:χ]*Model.O[:χb]
is_free_energy_definite(Model::BAMHLT0Mixtures) = (Delta(Model)>0)

function Deltas!(Model::BAMHLT0Mixtures)
    mineps = 1e-30;
    Model.aux["Δ"] = clamp(Delta(Model),mineps,Inf)
end

#### Conjugate order params SP
function exprOconj!(Model::BAMHLT0Mixtures)
    @extract Model : O Oconj aux params
    Deltas!(Model)
    Oconj[:P] = (1 + O[:χb]^2) / aux["Δ"]^2
    Oconj[:Pb] = (1 + O[:χ]^2) / aux["Δ"]^2  
end

function auxiliary!(Model::BAMHLT0Mixtures)
    Model.aux["A"] = sqrt(Model.params[:α] * Model.Oconj[:P]  / Model.params[:γ])
    Model.aux["Ab"] = sqrt(Model.params[:α] * Model.Oconj[:Pb] * Model.params[:γ])
end

Ξ(Model::BAMHLT0Mixtures, j::Int) =  Model.params[:γ] * Model.Signal.store1[j] * sqrt(0.5) / Model.aux["Ab"]
Ξb(Model::BAMHLT0Mixtures, j::Int) = (1 / Model.params[:γ]) * Model.Signal.store2[j] * sqrt(0.5) / Model.aux["A"]

exprM(Model::BAMHLT0Mixtures, j::Int) = erf(Ξb(Model,j))
exprChi(Model::BAMHLT0Mixtures,j::Int) = sqrt(2/π) * (1/Model.aux["A"]) * exp( - Ξb(Model,j)^2)


exprMb(Model::BAMHLT0Mixtures, j::Int) = erf(Ξ(Model,j))
exprChib(Model::BAMHLT0Mixtures,j::Int) = sqrt(2/π) * (1/Model.aux["Ab"]) * exp( - Ξ(Model,j)^2)


function lhs!(Model::BAMHLT0Mixtures{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)
    @assert length(unique(avg.probs))==1
    @assert avg.values[1] == - avg.values[2]
    exprOconj!(Model)
    auxiliary!(Model)
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)

    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,j)
        Onew[:χ] += avg.probs[j] * exprChi(Model,j)

        Onew[:Mb] += avg.probs[j] * avg.values[j] * exprMb(Model,j)
        Onew[:χb] += avg.probs[j] * exprChib(Model,j)
    end
end

function lhs!(Model::BAMHLT0Mixtures{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)
    exprOconj!(Model)
    auxiliary!(Model)
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)

    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,j)
        Onew[:χ] += avg.probs[j] * exprChi(Model,j)

        Onew[:Mb] += avg.probs[j] * (avg.values[j] / avg.L) * exprMb(Model,j)
        Onew[:χb] += avg.probs[j] * exprChib(Model,j)
    end
end


function lhs!(Model::BAMHLT0Mixtures{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)
    @assert length(unique(avg.probs))==1
    @assert avg.values[1] == - avg.values[2]
    exprOconj!(Model)
    auxiliary!(Model)
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)

    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates

        for μ = 1:avg.L
            Onew[Symbol(:Mb,μ)] += avg.probs[j] * avg.values[j,μ] * exprMb(Model,j)
            Onew[Symbol(:M,μ)] += avg.probs[j] * avg.values[j,μ] * exprM(Model,j)
            Onew[:χ] += avg.probs[j] * exprChi(Model,j)
            Onew[:χb] += avg.probs[j] * exprChib(Model,j)

        end

    end
end
########################################################## FREE ENERGY #########################################################


function exprEntropy(Model::BAMHLT0Mixtures, j::Int) 
    @extract Model : params O Oconj aux Signal


    IS = aux["A"]*sqrt(2/π) *exp(-Ξb(Model,j)^2) + Signal.store2[j]*erf( Ξb(Model,j) )
    IS_bar = aux["Ab"]*sqrt(2/π) *exp(-Ξ(Model,j)^2) + Signal.store1[j]*erf( Ξ(Model,j) )


    return params[:γ]*IS + IS_bar / params[:γ]

end

function free_energy(Model::BAMHLT0Mixtures,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj avg Signal
    exprOconj!(Model)
    auxiliary!(Model)

    energy_chi = 0.5 * params[:α] * ( Oconj[:P]*O[:χ]  +  Oconj[:Pb]*O[:χb] )
    other =  0.5 * params[:α] * ( O[:χ] + O[:χb] )  / Model.aux["Δ"]

    entropy = sum(avg.probs[j] * exprEntropy(Model, j) for j ∈ 1:avg.Nstates)
    return energySignal(Signal,avg,O) + energy_chi -other - entropy

    # there is still an incompatibility between this free energy and the one of BAMHLT0RSstandard or γ ≠ 1
end

_isMixedState(Model::BAMHLT0Mixtures{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-4 && Model.O[:Mb] > 1e-4
_isMixedState(Model::BAMHLT0Mixtures{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-4 && Model.O[:Mb] > 1e-4
_isMixedState(Model::BAMHLT0Mixtures{VectorVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = all(vcat(Model.O[1:Model.avg.L],Model.O[Model.avg.L+2:1+2*Model.avg.L]) .> 1e-4 )
