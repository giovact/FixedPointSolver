struct BAMHLT0RSstandard <: FPModel
    params::NamedTuple{(:α,:γ),Tuple{FT,FT}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors  
BAMHLT0RSstandard() = BAMHLT0RSstandard(0.1, 1.0, ones(4)) # dummy constructor for testing
BAMHLT0RSstandard(α,γ,O) = BAMHLT0RSstandard((α=α,γ=γ),NamedArray(O,[:M;:χ;:Mb;:χb]), NamedArray(zeros(2),[:P;:Pb]), Dict{String,FT}())
BAMHLT0RSstandard(α,γ) =BAMHLT0RSstandard(α,γ, zeros(4))

create_model(Model::BAMHLT0RSstandard, params) = BAMHLT0RSstandard(params, Model.O, Model.Oconj, copy(Model.aux)) # maybe 

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

Delta(Model::BAMHLT0RSstandard) = 1 - Model.O[:χ]*Model.O[:χb]
is_free_energy_definite(Model::BAMHLT0RSstandard) = (Delta(Model)>0)

function Deltas!(Model::BAMHLT0RSstandard)
    mineps = 1e-30;
    Model.aux["Δ"] = clamp(Delta(Model),mineps,Inf)
end

#### Conjugate order params SP
function exprOconj!(Model::BAMHLT0RSstandard)
    @extract Model : O Oconj aux params
    Deltas!(Model)

    Oconj[:P] = (1 + O[:χb]^2) / aux["Δ"]^2
    Oconj[:Pb] = (1 + O[:χ]^2) / aux["Δ"]^2  

end

# order params SP
function exprM(Model::BAMHLT0RSstandard,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:Mb] == 0.0 && return 0.0
    γb = 1/params[:γ]
    A = sqrt(params[:α]* γb * Oconj[:P])
    B = O[:Mb]*γb

    return erf( B / (A*sqrt(2)) )
end

function exprChi(Model::BAMHLT0RSstandard,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γb = 1/params[:γ]
    A = sqrt(params[:α]* γb * Oconj[:P])
    B = O[:Mb]*γb
    
    return sqrt(2/π) * (1/A) * exp(-B^2 / (2*A^2))
end

################################################################# LAYER BAR ##########################################################################

function exprMbar(Model::BAMHLT0RSstandard,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    
    O[:M] == 0.0 && return 0.0

    γ = params[:γ]
    A_bar = sqrt(params[:α]* γ * Oconj[:Pb])
    B_bar = O[:M]*γ

    return erf( B_bar / (A_bar*sqrt(2)) )
end


function exprChibar(Model::BAMHLT0RSstandard,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γ = params[:γ]
    A_bar = sqrt(params[:α]* γ * Oconj[:Pb])
    B_bar = O[:M]*γ
    
    return sqrt(2/π) * (1 / A_bar) * exp(-B_bar^2 / (2*A_bar^2))
end



function lhs!(Model::BAMHLT0RSstandard,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconj!(Model)

    Onew[:M] = exprM(Model,X)
    Onew[:χ] = exprChi(Model,X)

    Onew[:Mb] = exprMbar(Model,X)
    Onew[:χb] = exprChibar(Model,X)
end

########################################################## FREE ENERGY #########################################################


function free_energy(Model::BAMHLT0RSstandard,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    exprOconj!(Model)

    energy_M = O[:M]*O[:Mb]

    energy_chi = 0.5 * params[:α] * ( Oconj[:P]*O[:χ]  +  Oconj[:Pb]*O[:χb] )
    other =  0.5 * params[:α] * ( O[:χ] + O[:χb] )  / Model.aux["Δ"]

    γb = 1/params[:γ]
    A = sqrt(params[:α]* γb * Oconj[:P])
    B = O[:Mb]*γb

    IS = A*sqrt(2/π) *exp(-B^2 / (2*A^2)) + B*erf( B / (A*sqrt(2)) )


    γ = params[:γ]
    A_bar = sqrt(params[:α]* γ * Oconj[:Pb])
    B_bar = O[:M]*γ

    IS_bar = A_bar*sqrt(2/π) *exp(-B_bar^2 / (2*A_bar^2)) + B_bar*erf( B_bar / (A_bar*sqrt(2)) )

    return energy_M + energy_chi  - other - γ*IS - γb*IS_bar
    # there is still an incompatibility between this free energy and the one of BAMHLT0Mixtures or γ ≠ 1
end