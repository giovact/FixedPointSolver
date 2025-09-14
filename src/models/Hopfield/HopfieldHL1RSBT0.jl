struct HopfieldHL1RSBT0 <: FPModel
    params::NamedTuple{(:α,:Θ),Tuple{FT, FT}}
    O::NamedVec #vector of orderparameters -> M
    Oconj::NamedVec #vector of orderparameters -> M
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopfieldHL1RSBT0() = HopfieldHL1RSBT0(.1, 1.0,[1.0; 1e-2; 1e-3]) # dummy constructor for testing
HopfieldHL1RSBT0(α,Θ,O) = HopfieldHL1RSBT0((α=α,Θ=Θ),NamedArray(O,[:M;:χ;:ΔQ]),NamedArray(zeros(2),[:P1;:P2]),Dict{String,FT}())
HopfieldHL1RSBT0(α,Θ) = HopfieldHL1RSBT0(α,Θ,zeros(3))

create_model(Model::HopfieldHL1RSBT0, params) = HopfieldHL1RSBT0(params, Model.O, Model.Oconj,Dict{String,FT}())

#function enforce_initial_condition!(Model::HopfieldHL1RSBT0)
#    mineps = 1e-15
#    Model.O[:χ] = clamp(Model.O[:χ], mineps, 1-mineps)
#    Model.O[:ΔQ] = clamp(Model.O[:ΔQ], mineps, 1-mineps)
#    Model.O[:ΔQ] = clamp(Model.O[:ΔQ], mineps, (1-Model.O[:χ])/Model.params[:Θ] - mineps)
#end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::HopfieldHL1RSBT0) = 1 - Model.O[:χ]
Deltatheta(Model::HopfieldHL1RSBT0) = 1 - Model.O[:χ] - Model.params[:Θ]*Model.O[:ΔQ] 
is_free_energy_definite(Model::HopfieldHL1RSBT0) = (Delta(Model)>0) && (Deltatheta(Model) >0)

function Deltas!(Model::HopfieldHL1RSBT0)
    mineps = 1e-30;
    Model.aux["Δ"] = clamp(Delta(Model),mineps,Inf)
    Model.aux["Δtheta"] = clamp(Deltatheta(Model),mineps,Inf)
end

#### Conjugate order params SP
function exprOconj!(Model::HopfieldHL1RSBT0)
    @extract Model : O Oconj aux params
    Deltas!(Model)
    #P1
    Oconj[:P1] = (1 - O[:ΔQ]) / aux["Δtheta"]^2 
    #P2 
    Oconj[:P2] =  Oconj[:P1] +  Model.O[:ΔQ] / (aux["Δ"]*aux["Δtheta"]) 
end

# order params SP
function exprM(Model::HopfieldHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    a = sqrt(params[:α]*(Oconj[:P2] - Oconj[:P1]) )
    b(τ) = O[:M] + sqrt( params[:α] * Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*a) + b(τ)/(sqrt(2)*a)
    km(τ) = (sqrt(0.5)*params[:Θ]*a) - b(τ)/(sqrt(2)*a)
    
    F(τ) = ( 1 + erf(kp(τ)) )  / ( 1 + erf(km(τ)) )
    IM(τ) = 1 / ( 1 + F(τ) * exp(2*params[:Θ]*b(τ)) )
    return 1 - 2 * integrate(X,IM)
end

function exprChi(Model::HopfieldHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj

    a = sqrt(params[:α]*(Oconj[:P2] - Oconj[:P1]) )
    b(τ) = O[:M] + sqrt( params[:α] * Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*a) + b(τ)/(sqrt(2)*a)
    km(τ) = (sqrt(0.5)*params[:Θ]*a) - b(τ)/(sqrt(2)*a)
    
    Ichi(τ) =  exp(-b(τ)^2 / (2*a^2)) / ( exp(params[:Θ]* b(τ) ) * ( 1 + erf(kp(τ)) ) + exp(-params[:Θ]* b(τ))* ( 1 + erf(km(τ)) )   )
    
    return exp(-0.5*params[:Θ]^2*a^2)* sqrt(2/π) * (2 / a) * integrate(X,Ichi)
end

function exprdeltaQ(Model::HopfieldHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    
    a = sqrt(params[:α]*(Oconj[:P2] - Oconj[:P1]) )
    b(τ) = O[:M] + sqrt( params[:α] * Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*a) + b(τ)/(sqrt(2)*a)
    km(τ) = (sqrt(0.5)*params[:Θ]*a) - b(τ)/(sqrt(2)*a)

    IΔQ(τ) = ( ( 1 + erf(kp(τ)) ) * ( 1 + erf(km(τ)) ) )  / ( exp(params[:Θ]* b(τ) ) * ( 1 + erf(kp(τ)) ) + exp(-params[:Θ]*b(τ))*( 1 + erf(km(τ)) )   ) ^2
    return 4 * integrate(X,IΔQ)
end

function lhs!(Model::HopfieldHL1RSBT0,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconj!(Model)
    Onew[:M] = exprM(Model,X)
    Onew[:χ] = exprChi(Model,X)
    Onew[:ΔQ] = exprdeltaQ(Model,X)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::HopfieldHL1RSBT0,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    exprOconj!(Model)

    energy = 0.5*O[:M]^2
    energy_overlaps =  0.5 * params[:α] * ( Oconj[:P2]*O[:χ] + params[:Θ] *Oconj[:P1]*O[:ΔQ])

    logterm = (0.5*params[:α]/params[:Θ]) * ( log( Model.aux["Δtheta"]) - log( Model.aux["Δ"]) )
    
    other =  (0.5*params[:α]*(1-O[:ΔQ]))  / Model.aux["Δtheta"]

    a = sqrt(params[:α]*(Oconj[:P2] - Oconj[:P1]) )
    b(τ) = O[:M] + sqrt( params[:α] * Oconj[:P1] ) * τ

    kp(τ) = (sqrt(0.5)*params[:Θ]*a) + b(τ)/(sqrt(2)*a)
    km(τ) = (sqrt(0.5)*params[:Θ]*a) - b(τ)/(sqrt(2)*a)

    Gp(τ) =  exp(params[:Θ]* b(τ) ) * ( 1 + erf(kp(τ)) )  +  exp( -params[:Θ]*b(τ))* (1+erf(km(τ)) )
    IS(τ) = log( Gp(τ)  ) 

    return energy + energy_overlaps + logterm - other + log(2)/params[:Θ] - integrate(X,IS)/params[:Θ]
end


_isRetrieval(Model::HopfieldHL1RSBT0) = Model.O[:M] > 0.8
