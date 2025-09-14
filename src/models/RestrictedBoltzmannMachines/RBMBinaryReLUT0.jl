struct RBMBinaryReLUT0{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:α,:p,:η,:g,:λ),Tuple{FT,FT,FT,FT,FT}}
    Signal::RT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
RBMBinaryReLUT0(Signal::RT) where RT <:DiscreteExpectation = RBMBinaryReLUT0(.1,1.0,0.0,0.0,1.0,Signal,ones(Signal.L+1)) # dummy constructor for testing
RBMBinaryReLUT0(α,p,η,g,λ,Signal) = RBMBinaryReLUT0(α,p,η,g,λ,Signal, ones(Signal.L+1))
RBMBinaryReLUT0(α,p,η,g,λ,Signal,O) = RBMBinaryReLUT0((α=α,p=p,η=η,g=g,λ=λ),Signal,vcat(setSignal_OP(O,Signal),NamedArray(O[end:end],[:χ])),vcat(setSignal_OP(zeros(Signal.L+1),Signal),NamedArray(zeros(1),[:Q])), Dict{String,FT}())

create_model(Model::RBMBinaryReLUT0, params) = RBMBinaryReLUT0(params, Model.Signal,Model.O, Model.Oconj, Dict{String,FT}())

sanity_check_nop!(Model::RBMBinaryReLUT0{SingleVar}) = @assert n_order_params(Model) == 2
sanity_check_nop!(Model::RBMBinaryReLUT0{SumVar}) = @assert n_order_params(Model) == 2
sanity_check_nop!(Model::RBMBinaryReLUT0{VectorVar}) = @assert n_order_params(Model) == Model.Signal.L + 1

function sanity_checks!(Model::RBMBinaryReLUT0)
    #@assert Model.Signal.ξ0 == [-1;0;1]
    #@assert  Model.Signal.p0[1] == 0.5*Model.params[:p]
    #@assert  Model.Signal.p0[2] == 1-Model.params[:p]
    #@assert  Model.Signal.p0[3] == 0.5*Model.params[:p]
    sanity_check_nop!(Model)
end


#enforce_initial_condition!(Model::RBMBinaryReLUT0{SumVar}) = Model.O[:M] = 1 / Model.Signal.L
#function enforce_initial_condition!(Model::RBMBinaryReLUT0{VectorVar})
#    for μ=1:Model.Signal.L 
#        Model.O[Mμ(μ)] = rand()*(1-Model.Signal.Mmin) + Model.Signal.Mmin
#    end
#end


########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::RBMBinaryReLUT0) =  Model.params[:λ] - Model.params[:p] * Model.O[:χ]
Deltas!(Model::RBMBinaryReLUT0) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::RBMBinaryReLUT0) = Delta(Model)>0

function aux_variables!(Model::RBMBinaryReLUT0)
    @extract Model : O Oconj params aux
    Deltas!(Model)
    aux["prefactor_Qhat"] =  (  (params[:p] + params[:η]^2)*Herf(-params[:η]/sqrt(params[:p])) + params[:η] *sqrt(params[:p] / (2*π) )* exp(-0.5*params[:η]^2 / params[:p]) )  
end

#### Conjugate order params SP
exprMhat!(Model::RBMBinaryReLUT0{SingleVar}) = Model.Oconj[:M] = max(0.0, Model.O[:M] / Model.params[:λ])
exprMhat!(Model::RBMBinaryReLUT0{SumVar}) = Model.Oconj[:M] = max(0.0, Model.O[:M] / Model.params[:λ]) 

function exprMhat!(Model::RBMBinaryReLUT0{VectorVar})
    for μ = 1:Model.Signal.L
        Model.Oconj[Mμ(μ)] = max(0.0, Model.O[Mμ(μ)] / Model.params[:λ])
    end
end
exprQhat(Model::RBMBinaryReLUT0) = (Model.params[:p] / Model.aux["Δ"]^2 ) * Model.aux["prefactor_Qhat"]


function exprOconj!(Model::RBMBinaryReLUT0)
    aux_variables!(Model)
    exprMhat!(Model) # that depends on the signal term
    Model.Oconj[:Q] = exprQhat(Model)
   # typeof(Model.Signal)!= PureStates && println(Model.Oconj)
end

# order params SP

# in these functions j is the index of the component over all the possible configurations of the pattern to be averaged over
Signaldot(Model::RBMBinaryReLUT0{SingleVar},j::Int) = Model.O[:M] * Model.Signal.values[j]
Signaldot(Model::RBMBinaryReLUT0{SumVar},j::Int) = Model.O[:M] * Model.Signal.values[j]
Signaldot(Model::RBMBinaryReLUT0{VectorVar},j::Int) = sum(Model.O[Mμ(μ)] * Model.Signal.values[j,μ] for μ=1:Model.Signal.L)

exprM(Model::RBMBinaryReLUT0,j::Int) = erf( (Model.params[:g] +  Signaldot(Model,j)) / sqrt(2 * Model.params[:α] * Model.Oconj[:Q]) )
exprchi(Model::RBMBinaryReLUT0,j::Int) = sqrt(2/ ( π * Model.params[:α] * Model.Oconj[:Q]) ) * exp(- 0.5 * (Model.params[:g] +  Signaldot(Model,j))^2 / ( Model.params[:α] * Model.Oconj[:Q] ))

function lhs!(Model::RBMBinaryReLUT0{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model)
    @assert iszero(Onew)

    for j=1:Signal.Nstates
        Onew[:M] += Signal.probs[j] * Signal.values[j] * exprM(Model,j)
        Onew[:χ] += Signal.probs[j] * exprchi(Model,j)
    end
    @assert Onew[:M] in Interval(-1.0,1.0)
end

function lhs!(Model::RBMBinaryReLUT0{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model)
    @assert iszero(Onew)
    
    for j=1:Signal.Nstates
        Onew[:M] += Signal.probs[j] * (Signal.values[j] / Signal.L) * exprM(Model,j) 
        Onew[:χ] += Signal.probs[j] * exprchi(Model,j)
    end
    #println(Onew[:M], "  ", Onew[:χ])
end

function lhs!(Model::RBMBinaryReLUT0{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj Signal
    exprOconj!(Model)
    @assert iszero(Onew)
    
    for j=1:Signal.Nstates
        for μ = 1:Signal.L
            Onew[Mμ(μ)] += Signal.probs[j] * Signal.values[j,μ] * exprM(Model,j)
        end
        Onew[:χ] += Signal.probs[j] * exprchi(Model,j)
    end
end

########################################################## FREE ENERGY #########################################################
energySignal(Model::RBMBinaryReLUT0{SingleVar}) = Model.O[:M] * Model.Oconj[:M] - 0.5*(Model.O[:M]^2 / Model.params[:λ]) * (Model.O[:M]>0.0) 
energySignal(Model::RBMBinaryReLUT0{SumVar}) = Model.Signal.L * Model.O[:M] * Model.Oconj[:M] - 0.5*Model.Signal.L *(Model.O[:M]^2 / Model.params[:λ]) * (Model.O[:M]>0.0) 
energySignal(Model::RBMBinaryReLUT0{VectorVar}) = sum(Model.O[Mμ(μ)] * Model.Oconj[Mμ(μ)] - 0.5*(Model.O[Mμ(μ)]^2 / Model.params[:λ]) * (Model.O[Mμ(μ)]>0.0)  for μ =1:Model.Signal.L )


function exprEntropy(Model::RBMBinaryReLUT0,j::Int) 
    @extract Model : params O Oconj aux Signal
    dott = Signaldot(Model,j)
    sqrt(params[:α] * Oconj[:Q]) * sqrt(2/π) * exp(- 0.5 * (params[:g] + dott  )^2 / ( params[:α] * Oconj[:Q] ))  + (params[:g] +  dott)*erf( (params[:g] +  dott) / sqrt(2 * params[:α] * Oconj[:Q]) )
end

function free_energy(Model::RBMBinaryReLUT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux Signal
    exprOconj!(Model)

    energy_overlaps =  0.5*params[:α]*Oconj[:Q] * O[:χ]

    other1 = 0.5* (params[:α] / aux["Δ"]) * (params[:p] + params[:η]^2) * Herf( - params[:η] / sqrt(params[:p]))
    other2 = 0.5* (params[:α] / aux["Δ"]) * sqrt(params[:p] / (2*π) ) * params[:η] * exp(-0.5 * params[:η]^2 / params[:p] )

    entropy = 0.0
    for j=1:Signal.Nstates
        entropy += Signal.probs[j] * exprEntropy(Model,j)
    end

    return  energySignal(Model) + energy_overlaps - other1 - other2 - entropy
end

############################## Retrieval CHECK ########################################
isRetrieval(Model::RBMBinaryReLUT0{SingleVar},X::TI,Oold::NamedVec) where {TI <: IntegrationMethod} = Model.O[:M] /Model.params[:p] > 1e-2
isRetrieval(Model::RBMBinaryReLUT0{SumVar},X::TI,Oold::NamedVec) where {TI <: IntegrationMethod} = Model.O[:M] /Model.params[:p] > 1e-2
function isRetrieval(Model::RBMBinaryReLUT0{VectorVar},X::TI,Oold::NamedVec) where {TI <: IntegrationMethod}
    s = 0
    for μ = 1:Model.Signal.L 
        s+= (Model.O[Mμ(μ)] /Model.params[:p] > 1e-2)
    end
    return s == Model.Signal.L
end

     

function free_energy_SG_analytical(Model::RBMBinaryReLUT0{SingleVar})
    ModelSG = deepcopy(Model)

    @extract Model : params
    aux_variables!(ModelSG)
    
    
    ModelSG.O[:χ] =  params[:λ] / (params[:p] + sqrt(0.5 * π * params[:α] * ModelSG.aux["prefactor_Qhat"]))
    ModelSG.O[:M] = 0.0
    
    ModelSG.Oconj[:M] = 0.0
    
    Deltas!(ModelSG) # to recompute Δ with the newvalue of χ
    ModelSG.Oconj[:Q] = params[:p] * ModelSG.aux["prefactor_Qhat"] / ModelSG.aux["Δ"]^2
    
    return free_energy(ModelSG,NoInt())
end