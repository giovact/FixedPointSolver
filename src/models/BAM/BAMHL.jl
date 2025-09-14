struct BAMHL <: FPModel
    params::NamedTuple{(:β,:α,:γ),Tuple{FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BAMHL() = BAMHL(0.1, 0.1, 1.0, ones(4)) # dummy constructor for testing
BAMHL(β,α,γ) = BAMHL(β,α,γ, zeros(4))
BAMHL(β,α,γ,O) = BAMHL((β=β,α=α,γ=γ),NamedArray(O,[:M;:Mb;:Q;:Qb]),NamedArray(zeros(2),[:P;:Pb]),Dict{String,FT}())

create_model(Model::BAMHL, params) = BAMHL(params, Model.O, Model.Oconj,Dict{String,FT}())

function set_initial_condition!(Model::BAMHL,soltype::Symbol; M0::FT=0.9, Q0::FT=0.9,tol::FT = 1e-4)
    vareps = tol
    T= 1/Model.params[:β]
    @assert soltype ∈ [:R; :SG]
    if T <= 1.0
        Model.O[:Q] = max(Q0, 1.0 - T + vareps)
        Model.O[:Qb] = clamp(max(Q0, 1.0 - T + vareps + rand()*vareps*10), 1 - 1e-4)
        Model.O[:M] = soltype==:SG ? 0.0 : max(M0, 1.0 - T + vareps)
        Model.O[:Mb] = soltype==:SG ? 0.0 : max(M0 + rand()*vareps*10, 1.0 - T + vareps)
    else
        Model.O[:Q] = Q0
        Model.O[:Q] = Q0 + rand()*vareps*10
        Model.O[:M] = soltype==:SG ? 0.0 : max(M0, 1.0 - T + vareps)
        Model.O[:Mb] = soltype==:SG ? 0.0 : max(M0 + rand()*vareps*10, 1.0 - T + vareps)
    end 
    return 
end


########################################################## SELF-CONSISTENT EQUATIONS #########################################################

Delta(Model::BAMHL) = 1 - Model.params[:β]^2 * (1 - Model.O[:Q]) * (1 - Model.O[:Qb])
Deltas!(Model::BAMHL) = Model.aux["Δ"] = clamp(Delta(Model),1e-30,Inf)
is_free_energy_definite(Model::BAMHL) = Delta(Model) >0
#### Conjugate order params SP
function exprOconj!(Model::BAMHL)
    @extract Model : O Oconj
    Deltas!(Model)
    #Qhat
    Oconj[:P] = ( O[:Qb] + Model.params[:β]^2 * O[:Q] * (1-O[:Qb])^2 ) / Model.aux["Δ"]^2
    Oconj[:Pb] = ( O[:Q] + Model.params[:β]^2 * O[:Qb] * (1-O[:Q])^2 ) / Model.aux["Δ"]^2
end

# order params SP
function exprM(Model::BAMHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:Mb]==0.0 && return 0.0
    IM(z) = tanh( params[:β] * ( sqrt( params[:α]*Oconj[:P] / params[:γ] )*z + O[:Mb]/params[:γ] ) )
    return integrate(X,IM)
end

function exprMb(Model::BAMHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:M]==0.0 && return 0.0
    IMb(z) = tanh( params[:β] * ( sqrt( params[:α]*Oconj[:Pb]*params[:γ] )*z + O[:M]*params[:γ] ) )
    return integrate(X,IMb)
end


function exprQ(Model::BAMHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    IQ(z) = tanh( params[:β] * ( sqrt( params[:α]*Oconj[:P] / params[:γ] )*z + O[:Mb]/params[:γ] ) ) ^ 2
    return integrate(X,IQ)
end

function exprQb(Model::BAMHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    IQb(z) = tanh( params[:β] * ( sqrt( params[:α]*Oconj[:Pb]*params[:γ] )*z + O[:M]*params[:γ] ) ) ^ 2
    return integrate(X,IQb)
end


function lhs!(Model::BAMHL,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O Oconj
    exprOconj!(Model)
    Onew[:M] = exprM(Model,X)
    Onew[:Mb] = exprMb(Model,X)
    Onew[:Q] = exprQ(Model,X)
    Onew[:Qb] = exprQb(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::BAMHL,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    exprOconj!(Model)
    Δ = Model.aux["Δ"]
    energy = O[:M]* O[:Mb]
    energy_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:P] + 0.5*params[:α]*params[:β]*(1-O[:Qb])*Oconj[:Pb]
    logdetK = (params[:α]/(2*params[:β])) * log(Δ)  - ( (params[:α]*params[:β])/(2*Δ)) * (O[:Q]*(1-O[:Qb]) + O[:Qb]*(1-O[:Q]))
    IS(z)  = log( 2 * cosh( params[:β] * ( sqrt( params[:α]*Oconj[:P] / params[:γ] )*z + O[:Mb]/params[:γ] ) ) )
    ISb(z) = log( 2 * cosh( params[:β] * ( sqrt( params[:α]*Oconj[:Pb] * params[:γ] )*z + O[:M]*params[:γ] ) ) )
    return  energy + energy_overlaps + logdetK - (params[:γ]/params[:β])*integrate(X,IS) - (1/(params[:γ]*params[:β]))*integrate(X,ISb) 
end


#_isRetrieval(Model::BAMHL) = Model.params[:γ] > 1  ? Model.O[:Mb] > 0.8 : Model.O[:M] > 0.8
