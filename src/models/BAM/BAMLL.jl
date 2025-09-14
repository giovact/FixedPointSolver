struct BAMLL{RT<:DiscreteExpectation} <: FPModel # for the moment patterns are decorrelated and binaryict
    params::NamedTuple{(:β, :γ),Tuple{FT,FT}}
    avg::RT
    Signal::RetrievalBAM
    O::NamedVec #vector of orderparameters -> M1, …, Ml, Mb1, …, Mbl
    aux::Dict{String,FT} # other stuff, among them free energy
end
# %TODO use Named Matrix instead of NamedArray


BAMLL(avg::RT) where {RT<:DiscreteExpectation} = BAMLL(2.0, 1.0, avg, ones(n_magnetizations_bothlayers(RetrievalBAM(avg),avg))) # dummy constructor for testing
BAMLL(β, γ, avg) = BAMLL(β, γ, avg, zeros(n_magnetizations_bothlayers(RetrievalBAM(avg),avg)))

function BAMLL(β, γ, avg, O)
    Signal = RetrievalBAM(avg)
    @assert length(O) == n_magnetizations_bothlayers(Signal,avg)
    no_M = n_magnetizations_singlelayer(Signal,avg) 
    return BAMLL((β=β, γ=γ), 
                    avg,
                    RetrievalBAM(avg),
                    vcat(setBAMsignal_L1(O[1:no_M], avg),setBAMsignal_L2(O[no_M+1:2*no_M], avg)), 
                    Dict{String,FT}())
end

create_model(Model::BAMLL, params) = BAMLL(params, Model.avg, Model.Signal,Model.O, Dict{String,FT}())

function sanity_checks!(Model::BAMLL)
    @assert n_order_params(Model) ==  n_magnetizations_bothlayers(Model.Signal,Model.avg)
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
χ(Model::BAMLL, j::Int) =  (Model.params[:β] * Model.params[:γ]) * Model.Signal.store1[j]
χb(Model::BAMLL, j::Int) = (Model.params[:β] / Model.params[:γ]) * Model.Signal.store2[j]

exprM(Model::BAMLL, j::Int) = tanh(χb(Model,j))
exprMb(Model::BAMLL, j::Int) = tanh(χ(Model,j))


function lhs!(Model::BAMLL{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)
    @assert length(unique(avg.probs))==1
    @assert avg.values[1] == - avg.values[2]
    
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)

    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates
        Onew[:Mb] += avg.probs[j] * avg.values[j] * exprMb(Model,j)
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,j)
    end
end

function lhs!(Model::BAMLL{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)    
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)
    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates
        Onew[:Mb] += avg.probs[j] * (avg.values[j] / avg.L) * exprMb(Model,j)
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,j)
    end
end


function lhs!(Model::BAMLL{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : Signal avg O 
    @assert iszero(Onew)    
    clean!(Signal)
    fill_Signal_BAM_L1!(Signal,avg,O)
    fill_Signal_BAM_L2!(Signal,avg,O)

    # for the moment parallel update, bipartite structure of the model not exploited
    for j=1:avg.Nstates
        for μ = 1:avg.L
            Onew[Symbol(:Mb,μ)] += avg.probs[j] * avg.values[j,μ] * exprMb(Model,j)
            Onew[Symbol(:M,μ)] += avg.probs[j] * avg.values[j,μ] * exprM(Model,j)
        end
    end

end


########################################################## FREE ENERGY #########################################################


exprEntropy(Model::BAMLL, j::Int) = Model.params[:γ] * log2cosh(χb(Model, j)) + (1 / Model.params[:γ]) * log2cosh(χ(Model, j))

function free_energy(Model::BAMLL, X::TI) where {TI<:IntegrationMethod}
    @extract Model: params Signal avg O 

    entropy = sum(avg.probs[j] * exprEntropy(Model, j) for j ∈ 1:avg.Nstates)
    return energySignal(Signal,avg,O) - entropy / params[:β]
end

