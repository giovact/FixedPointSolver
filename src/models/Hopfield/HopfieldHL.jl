struct HopfieldHL{RT<:DiscreteExpectation} <: FPModel
    params::NamedTuple{(:β,:α),Tuple{FT,FT}}
    avg::RT
    Signal::RetrievalHopfield
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopfieldHL(avg::RT) where RT <:DiscreteExpectation = HopfieldHL(.1, .1,avg,ones(avg.L+1)) # dummy constructor for testing
HopfieldHL(β,α,avg) = HopfieldHL(β,α,avg,zeros(avg.L+1))
HopfieldHL(β,α,avg,O) = HopfieldHL((β=β,α=α),avg,RetrievalHopfield(),vcat(setSignal_OP(O,avg),NamedArray(O[end:end],[:Q])),NamedArray(zeros(1),[:Q]), Dict{String,FT}())

create_model(Model::HopfieldHL, params) = HopfieldHL(params, Model.avg,Model.Signal,Model.O, Model.Oconj, Dict{String,FT}())

function sanity_checks!(Model::HopfieldHL)
    @assert n_order_params(Model) == 1 + n_magnetizations(Model.Signal,Model.avg)
    Model.aux["mineps"] = 1e-14
end

function enforce_initial_condition!(Model::HopfieldHL{VectorVar})
    for μ=1:Model.avg.L 
        #Model.O[Mμ(μ)] = rand()*(1-Model.avg.Mmin) + Model.avg.Mmin
    end
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::HopfieldHL) =  1 - Model.params[:β] * (1 - Model.O[:Q])
Deltas!(Model::HopfieldHL) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::HopfieldHL) = Delta(Model)>0

#### Conjugate order params SP
function exprOconj!(Model::HopfieldHL)
    @extract Model : O Oconj
    Deltas!(Model)
    #Qhat
    Oconj[:Q] = O[:Q] / Model.aux["Δ"]^2
end

# order params SP

χI(Model::HopfieldHL,z,j::Int) = Model.params[:β]* (Model.avg.store[j]  + sqrt(Model.params[:α]*Model.Oconj[:Q]) * z) 

function exprM(Model::HopfieldHL,X::TI,j::Int) where TI <: IntegrationMethod
    IM(z) = tanh( χI(Model,z,j) )
    integrate(X,IM)
end

function exprQ(Model::HopfieldHL,X::TI,j::Int) where TI <: IntegrationMethod
    IQ(z) = tanh( χI(Model,z,j) )^2
    integrate(X,IQ)
end

function expr_tanh_pow(Model::HopfieldHL,X::TI,pow::Int,j::Int) where TI <: IntegrationMethod
    I_tanh_pow(z) = tanh( χI(Model,z,j) )^pow
    integrate(X,I_tanh_pow)
end


function lhs!(Model::HopfieldHL{SingleVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg aux
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)

    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * avg.values[j] * exprM(Model,X,j)
        Onew[:Q] += avg.probs[j] * exprQ(Model,X,j)
    end
    #@assert Onew[:M] in Interval(-1.0,1.0)

    Onew[:M] = threshold_to_0(Onew[:M],aux["mineps"])
end

function lhs!(Model::HopfieldHL{SumVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg aux
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)
    
    for j=1:avg.Nstates
        Onew[:M] += avg.probs[j] * (avg.values[j] / avg.L) * exprM(Model,X,j) 
        Onew[:Q] += avg.probs[j] * exprQ(Model,X,j)
    end
    
    Onew[:M] = threshold_to_0(Onew[:M],aux["mineps"])
end

function lhs!(Model::HopfieldHL{VectorVar},X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj avg aux
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    @assert iszero(Onew)
    
    for j=1:avg.Nstates
        for μ = 1:avg.L
            Onew[Mμ(μ)] += avg.probs[j] * avg.values[j,µ] * exprM(Model,X,j)
        end
        Onew[:Q] += avg.probs[j] * exprQ(Model,X,j)
    end

    for μ = 1:avg.L
        Onew[Mμ(μ)] =  threshold_to_0(Onew[Mμ(μ)],aux["mineps"])
    end
end

########################################################## FREE ENERGY #########################################################
function exprEntropy(Model::HopfieldHL,X::TI,j::Int) where TI <: IntegrationMethod
    IQ(z) = log2cosh( χI(Model,z,j) )
    integrate(X,IQ)
end

function free_energy(Model::HopfieldHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux avg
    exprOconj!(Model)
    fill_avg!(Model.Signal,Model.avg,Model.O)
    energy_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:Q]

    logdetK = (params[:α]/(2*params[:β])) * log(aux["Δ"])  
    other =  0.5*params[:α]* O[:Q] / aux["Δ"] 

    entropy = 0.0
    for j=1:avg.Nstates
        entropy += avg.probs[j] * exprEntropy(Model,X,j)
    end

    return  energySignal(Model.Signal,Model.avg,Model.O) + energy_overlaps + logdetK - other -  entropy/params[:β]

end


_isRetrieval(Model::HopfieldHL{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isMixedState(Model::HopfieldHL{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2
_isMixedState(Model::HopfieldHL{VectorVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = all(Model.O[1:end-1] .> min(tol * 100,1e-2) )

function _isRetrievalStateDominant(Model::HopfieldHL{SingleVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod
    ModelSG = create_model(deepcopy(Model), Model.params)
    ModelSG.O[:Q] = 1.0 
    ModelSG.O[:M] = 0.0
    solve_SPequations!(ModelSG, X; tol = tol, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 10000,force_dict = Dict(:M=>0.0))
    return Model.aux["f"] < ModelSG.aux["f"]
end

function expr_dAT(Model::HopfieldHL,X::TI,j::Int) where TI <: IntegrationMethod
    I_dAT(z) = 1 / cosh( χI(Model,z,j) )^4
    integrate(X,I_dAT)
end


function eig_DAT(Model::HopfieldHL, X::TI) where TI <: IntegrationMethod
    @extract Model: avg O Oconj params aux
    λ_qblock= 0.0
    for j=1:avg.Nstates
        λ_qblock += avg.probs[j] * expr_dAT(Model,X,j)
    end
    aux["Δ"]^2 - params[:α]*params[:β]^2*λ_qblock
end
_is_RS_stable(Model::HopfieldHL, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = eig_DAT(Model,X)>0


############# stability of Symmetric Mixtures, analitical #################################

function compute_eigenvalues_Hessian_SYM(Model::HopfieldHL{SumVar},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux avg
    
    aux["KK"] = (0.5 / aux["Δ"]^2) - (params[:β] *  Oconj[:Q] / aux["Δ"])

    ro = 0.0
    am = 0.0
    aq = 0.0

    av_S2_tanh2 = 0.0
    av_S_tanh2 = 0.0
    av_tanh4 = 0.0
    for j=1:avg.Nstates
        av_S_tanh2 += avg.probs[j] * avg.values[j]  * expr_tanh_pow(Model,X,2,j)
        av_S2_tanh2 += avg.probs[j] * avg.values[j]^2  * expr_tanh_pow(Model,X,2,j)
        av_tanh4 += avg.probs[j] * expr_tanh_pow(Model,X,4,j)
    end


    aux["rd"] = deepcopy(aux["Δ"])
    aux["ro"] = ( params[:β] / (avg.L * (avg.L-1)) ) * ( av_S2_tanh2 - avg.L * O[:Q] )
    aux["am"] = params[:α] * params[:β]^2 * aux["KK"]  * av_S_tanh2 / avg.L
    aux["aq"] = params[:α] * params[:β] * aux["KK"] * (2 * params[:α] * params[:β]^2 * aux["KK"] * (1 - 4 * O[:Q] + 3 * av_tanh4) - 1)
end

function compute_aux!(Model::HopfieldHL{SumVar},X::TI) where TI <: IntegrationMethod
    @extract Model : aux avg

    compute_eigenvalues_Hessian_SYM(Model,X)
    
    aux["λ1"] = aux["rd"] - aux["ro"] 
    aux["sqrt_arg"] = 4 * avg.L * aux["am"]^2 + (aux["rd"] + (avg.L-1) * aux["ro"] - aux["aq"] )^2
    aux["λm"] = 0.5 * (aux["aq"] + aux["rd"] + (avg.L-1) * aux["ro"]  - sqrt(aux["sqrt_arg"]) ) 
    aux["λp"] = 0.5 * (aux["aq"] + aux["rd"] + (avg.L-1) * aux["ro"]  + sqrt(aux["sqrt_arg"]) ) 
    #aux["λ_noq"] = ( aux["rd"] + (avg.L-1) * aux["ro"] ) 
    aux["eig_DAT"] = eig_DAT(Model,X)
    
end



_isSymMixturestable(Model::HopfieldHL{SumVar}, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2 &&  Model.aux["λ1"] >0
