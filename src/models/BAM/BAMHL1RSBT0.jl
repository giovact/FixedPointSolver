struct BAMHL1RSBT0 <: FPModel
    params::NamedTuple{(:α,:γ,:Θ),Tuple{FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BAMHL1RSBT0() = BAMHL1RSBT0(0.1, 1.0, 1.0, [1.0; 1e-3; (1e-3) / 3 ; 1.0; 1e-3; (1e-3) / 4] ) # dummy constructor for testing

BAMHL1RSBT0(α,γ,Θ,O) = BAMHL1RSBT0((α=α,γ=γ,Θ=Θ),NamedArray(O,[:M;:χ;:ΔQ;:Mb;:χb;:ΔQb]), NamedArray(zeros(4),[:P1;:P2;:P1b;:P2b]), Dict{String,FT}())
BAMHL1RSBT0(α,γ,Θ) =BAMHL1RSBT0(α,γ,Θ, zeros(6))

create_model(Model::BAMHL1RSBT0, params) = BAMHL1RSBT0(params, Model.O, Model.Oconj, copy(Model.aux)) # maybe 

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
kappa(Model::BAMHL1RSBT0) = Model.O[:χ] + Model.params[:Θ]*Model.O[:ΔQ]
kappa_bar(Model::BAMHL1RSBT0) = Model.O[:χb] + Model.params[:Θ]*Model.O[:ΔQb]

Delta(Model::BAMHL1RSBT0) = 1 - Model.O[:χ]*Model.O[:χb]
Deltatheta(Model::BAMHL1RSBT0) = 1 - kappa(Model)*kappa_bar(Model)
is_free_energy_definite(Model::BAMHL1RSBT0) = (Delta(Model)>0) && (Deltatheta(Model) >0)

function Deltas!(Model::BAMHL1RSBT0)
    mineps = 1e-30;
    Model.aux["K"] = kappa(Model)
    Model.aux["Kb"] = kappa_bar(Model)
    Model.aux["Δ"] = clamp(Delta(Model),mineps,Inf)
    Model.aux["Δtheta"] = clamp(Deltatheta(Model),mineps,Inf)
end

#### Conjugate order params SP
function exprOconj!(Model::BAMHL1RSBT0)
    @extract Model : O Oconj aux params
    Deltas!(Model)

    Oconj[:P1] = ( (1 - O[:ΔQb]) + (1 - O[:ΔQ])*aux["Kb"]^2  ) / aux["Δtheta"]^2  
    Oconj[:P2] =  Oconj[:P1] +  ( O[:ΔQb] + aux["Kb"]*O[:χb]*O[:ΔQ] )/ (aux["Δ"]*aux["Δtheta"]) 

    Oconj[:P1b] = ( (1 - O[:ΔQ]) + (1 - O[:ΔQb])*aux["K"]^2  ) / aux["Δtheta"]^2  
    Oconj[:P2b] =  Oconj[:P1b] +  ( O[:ΔQ] + aux["K"]*O[:χ]*O[:ΔQb] )/ (aux["Δ"]*aux["Δtheta"]) 

end

# order params SP
function exprM(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γb = 1/params[:γ]

    A = sqrt(params[:α]* γb * (Oconj[:P2] - Oconj[:P1]) )
    B(τ) = O[:Mb]*γb + sqrt( params[:α] *γb* Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*A) + B(τ)/(sqrt(2)*A)
    km(τ) = (sqrt(0.5)*params[:Θ]*A) - B(τ)/(sqrt(2)*A)
    F(τ) = ( 1 + erf(kp(τ)) )  / ( 1 + erf(km(τ)) )

    IM(τ) = 1 / ( 1 + F(τ) * exp(2*params[:Θ]*B(τ)) )
    return 1 - 2 * integrate(X,IM)
end

function exprChi(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γb = 1/params[:γ]

    A = sqrt(params[:α]* γb * (Oconj[:P2] - Oconj[:P1]) )
    B(τ) = O[:Mb]*γb + sqrt( params[:α] *γb* Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*A) + B(τ)/(sqrt(2)*A)
    km(τ) = (sqrt(0.5)*params[:Θ]*A) - B(τ)/(sqrt(2)*A)
    
    
    Ichi(τ) =  exp(-B(τ)^2 / (2*A^2)) / ( exp(params[:Θ]* B(τ) ) * ( 1 + erf(kp(τ)) ) + exp(-params[:Θ]* B(τ))* ( 1 + erf(km(τ)) )   )
    
    return exp(-0.5*params[:Θ]^2*A^2)* sqrt(2/π) * (2 / A) * integrate(X,Ichi) # here the only duda is that there is a gamma bar somewhere in the prefactor
end

function exprdeltaQ(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    
    γb = 1/params[:γ]

    A = sqrt(params[:α]* γb * (Oconj[:P2] - Oconj[:P1]) )
    B(τ) = O[:Mb]*γb + sqrt( params[:α] *γb* Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*A) + B(τ)/(sqrt(2)*A)
    km(τ) = (sqrt(0.5)*params[:Θ]*A) - B(τ)/(sqrt(2)*A)

    IΔQ(τ) = ( ( 1 + erf(kp(τ)) ) * ( 1 + erf(km(τ)) ) )  / ( exp(params[:Θ]* B(τ) ) * ( 1 + erf(kp(τ)) ) + exp(-params[:Θ]*B(τ))*( 1 + erf(km(τ)) )   ) ^2
    return 4 * integrate(X,IΔQ)
end

################################################################# LAYER BAR ##########################################################################

function exprMbar(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γ = params[:γ]

    A_bar = sqrt(params[:α]* γ * (Oconj[:P2b] - Oconj[:P1b]) )
    B_bar(τ) = O[:M]*γ + sqrt( params[:α] *γ* Oconj[:P1b] ) * τ

    kp_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) + B_bar(τ)/(sqrt(2)*A_bar)
    km_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) - B_bar(τ)/(sqrt(2)*A_bar)
    F_bar(τ) = ( 1 + erf(kp_bar(τ)) )  / ( 1 + erf(km_bar(τ)) )

    IMbar(τ) = 1 / ( 1 + F_bar(τ) * exp(2*params[:Θ]*B_bar(τ)) )
    return 1 - 2 * integrate(X,IMbar)
end


function exprChibar(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    γ = params[:γ]

    A_bar = sqrt(params[:α]* γ * (Oconj[:P2b] - Oconj[:P1b]) )
    B_bar(τ) = O[:M]*γ + sqrt( params[:α] *γ* Oconj[:P1b] ) * τ

    kp_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) + B_bar(τ)/(sqrt(2)*A_bar)
    km_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) - B_bar(τ)/(sqrt(2)*A_bar)
      
    Ichibar(τ) =  exp(-B_bar(τ)^2 / (2*A_bar^2)) / ( exp(params[:Θ]* B_bar(τ) ) * ( 1 + erf(kp_bar(τ)) ) + exp(-params[:Θ]* B_bar(τ))* ( 1 + erf(km_bar(τ)) )   )
    
    return exp(-0.5*params[:Θ]^2*A_bar^2)* sqrt(2/π) * (2 / A_bar) * integrate(X,Ichibar) # here the only duda is that there is a gamma somewhere in the prefactor
end

function exprdeltaQbar(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    
    γ = params[:γ]

    A_bar = sqrt(params[:α]* γ * (Oconj[:P2b] - Oconj[:P1b]) )
    B_bar(τ) = O[:M]*γ + sqrt( params[:α] *γ* Oconj[:P1b] ) * τ

    kp_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) + B_bar(τ)/(sqrt(2)*A_bar)
    km_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) - B_bar(τ)/(sqrt(2)*A_bar)

    IΔQbar(τ) = ( ( 1 + erf(kp_bar(τ)) ) * ( 1 + erf(km_bar(τ)) ) )  / ( exp(params[:Θ]* B_bar(τ) ) * ( 1 + erf(kp_bar(τ)) ) + exp(-params[:Θ]*B_bar(τ))*( 1 + erf(km_bar(τ)) )   ) ^2
    return 4 * integrate(X,IΔQbar)
end


function lhs!(Model::BAMHL1RSBT0,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconj!(Model)

    Onew[:M] = exprM(Model,X)
    Onew[:χ] = exprChi(Model,X)
    Onew[:ΔQ] = exprdeltaQ(Model,X)

    Onew[:Mb] = exprMbar(Model,X)
    Onew[:χb] = exprChibar(Model,X)
    Onew[:ΔQb] = exprdeltaQbar(Model,X)
end

########################################################## FREE ENERGY #########################################################


function free_energy(Model::BAMHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    exprOconj!(Model)

    energy_M = O[:M]*O[:Mb]

    energy_chi = 0.5 * params[:α] * ( Oconj[:P2]*O[:χ]  +  Oconj[:P2b]*O[:χb] )
    energy_deltaq =  0.5 * params[:α] * params[:Θ] * ( Oconj[:P1]*O[:ΔQ] + Oconj[:P1b]*O[:ΔQb] )
    
    logterm = ((0.5*params[:α])/params[:Θ]) * ( log( Model.aux["Δtheta"]) - log( Model.aux["Δ"]) )
    
    other =  0.5 * params[:α] * ( (1-O[:ΔQb])*Model.aux["K"] +(1-O[:ΔQ])*Model.aux["Kb"] )  / Model.aux["Δtheta"]


    γb = 1/params[:γ]
    A = sqrt(params[:α]* γb * (Oconj[:P2] - Oconj[:P1]) )
    B(τ) = O[:Mb]*γb + sqrt( params[:α] *γb* Oconj[:P1] ) * τ
    kp(τ) = (sqrt(0.5)*params[:Θ]*A) + B(τ)/(sqrt(2)*A)
    km(τ) = (sqrt(0.5)*params[:Θ]*A) - B(τ)/(sqrt(2)*A)

    IS(τ) = log( exp(params[:Θ]* B(τ) ) * ( 1 + erf(kp(τ)) ) + exp( -params[:Θ]*B(τ))* (1+erf(km(τ)) )    ) 



    γ = params[:γ]
    A_bar = sqrt(params[:α]* γ * (Oconj[:P2b] - Oconj[:P1b]) )
    B_bar(τ) = O[:M]*γ + sqrt( params[:α] *γ* Oconj[:P1b] ) * τ
    kp_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) + B_bar(τ)/(sqrt(2)*A_bar)
    km_bar(τ) = (sqrt(0.5)*params[:Θ]*A_bar) - B_bar(τ)/(sqrt(2)*A_bar)

    IS_bar(τ) = log( exp(params[:Θ]* B_bar(τ) ) * ( 1 + erf(kp_bar(τ)) )  +  exp( -params[:Θ]*B_bar(τ))* (1+erf(km_bar(τ)) )    ) 

    return energy_M + energy_chi + energy_deltaq + logterm - other + (γ+γb)*(log(2.0)/params[:Θ]) - ( γ/params[:Θ])*integrate(X,IS) - ( γb/params[:Θ])*integrate(X,IS_bar)
    
end


# verify phase
#_isRetrieval(Model::BAMHL1RSBT0) = Model.params[:γ] > 1  ? Model.O[:Mb] > 0.8 : Model.O[:M] > 0.8
