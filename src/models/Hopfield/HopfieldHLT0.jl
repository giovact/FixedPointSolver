struct HopfieldHLT0{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:α,),Tuple{FT}}
    avg::RT
    Signal::RetrievalHopfield
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopfieldHLT0(avg::RT) where RT <:DiscreteExpectation = HopfieldHLT0(0.1,avg,vcat(ones(avg.L),[1e-2])) # dummy constructor for testing
HopfieldHLT0(α,avg,O) = HopfieldHLT0((;α=α),avg,RetrievalHopfield(),vcat(setSignal_OP(O,avg),NamedArray(O[end:end],[:χ])),NamedArray(zeros(1),[:Q]), Dict{String,FT}())

create_model(Model::HopfieldHLT0, params) = HopfieldHLT0(params, Model.avg,Model.Signal,Model.O, Model.Oconj, Dict{String,FT}())

function sanity_checks!(Model::HopfieldHLT0)
    @assert n_order_params(Model) == 1 + n_magnetizations(Model.Signal,Model.avg)
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::HopfieldHLT0) =  1 - Model.O[:χ]
Deltas!(Model::HopfieldHLT0) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::HopfieldHLT0) = Delta(Model)>0


#### Conjugate order params SP
exprQhat(Model::HopfieldHLT0) = (1 / Model.aux["Δ"]^2 )

function exprOconj!(Model::HopfieldHLT0)
    Deltas!(Model)
    Model.Oconj[:Q] = exprQhat(Model)
end

# order params SP

# in these functions j is the index of the component over all the possible configurations of the pattern to be averaged over
exprM(Model::HopfieldHLT0,j::Int) = erf(Model.avg.store[j] / sqrt(2 * Model.params[:α] * Model.Oconj[:Q]) )
exprchi(Model::HopfieldHLT0,j::Int) = sqrt(2/ ( π * Model.params[:α] * Model.Oconj[:Q]) ) * exp(- 0.5 * (Model.avg.store[j])^2 / ( Model.params[:α] * Model.Oconj[:Q] ))

function lhs!(Model::HopfieldHLT0{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)

    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,j)
        Onew[:χ] += avg.probs[j] * exprchi(Model,j)
    end

end

function lhs!(Model::HopfieldHLT0{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)
    
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,j) 
        Onew[:χ] += avg.probs[j] * exprchi(Model,j)
    end
end

function lhs!(Model::HopfieldHLT0{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)
    
    for j=1:avg.Nstates
        for μ = 1:avg.L
            Onew[Mμ(μ)] += avg.probs[j] * avg.values[j,μ] * exprM(Model,j)
        end
        Onew[:χ] += avg.probs[j] * exprchi(Model,j)
    end
end

########################################################## FREE ENERGY #########################################################

function exprEntropy(Model::HopfieldHLT0,j::Int) 
    @extract Model : params O Oconj aux avg
    dott = avg.store[j]
    sqrt(params[:α] * Oconj[:Q]) * sqrt(2/π) * exp(- 0.5 * (dott)^2 / ( params[:α] * Oconj[:Q] ))  + (dott)*erf( dott / sqrt(2 * params[:α] * Oconj[:Q]) )
end

function free_energy(Model::HopfieldHLT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux avg
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)

    energy_overlaps =  0.5*params[:α]*Oconj[:Q] * O[:χ]

    other = 0.5* (params[:α] / aux["Δ"])

    entropy = 0.0
    for j=1:avg.Nstates
        entropy += avg.probs[j] * exprEntropy(Model,j)
    end

    return  energySignal(Model.Signal,Model.avg,Model.O)  + energy_overlaps  - other - entropy
end

############################## Retrieval CHECK ########################################
isRetrieval(Model::HopfieldHLT0{SingleVar}) = Model.O[:M] > 1e-2
isRetrieval(Model::HopfieldHLT0{SumVar}) = Model.O[:M] > 1e-2
function isRetrieval(Model::HopfieldHLT0{VectorVar})
    s = 0
    for μ = 1:Model.avg.L 
        s+= (Model.O[Mμ(μ)] > 1e-2)
    end
    return s == Model.avg.L
end

_isMixedState(Model::HopfieldHLT0{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-4
_isMixedState(Model::HopfieldHLT0{VectorVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = all(Model.O[1:end-1] .> 1e-4 )