struct HopfieldHL1RSB <: FPModel
    params::NamedTuple{(:β,:α,:θ),Tuple{FT,FT,FT}}
    O::NamedVec #vector of orderparameters -> M, Q, P
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat, Phat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
HopfieldHL1RSB() = HopfieldHL1RSB(.1, .1,.1,ones(3)) # dummy constructor for testing
HopfieldHL1RSB(β,α,θ) = HopfieldHL1RSB(β,α,θ, zeros(3))
HopfieldHL1RSB(β,α,θ,O) = HopfieldHL1RSB((β=β,α=α,θ=θ),NamedArray(O,[:M,:Q1,:Q2]),NamedArray(zeros(2),[:P1,:P2]),Dict{String,FT}())

create_model(Model::HopfieldHL1RSB, params) = HopfieldHL1RSB(params, Model.O,Model.Oconj,Dict{String,FT}())

#function enforce_initial_condition!(Model::HopfieldHL1RSB)
 #   T = 1/Model.params[:β]
  #  tol = 1e-4
   # if Model.O[:Q2] < 1.0 - T
   #     Model.O[:Q2] = min(1.0, 1.0 - T + 5*tol)
   # end
   # delta = Delta(Model)
   # if Model.O[:Q1] < Model.O[:Q2] - (T*delta)/Model.params[:θ]
   #     Model.O[:Q1] = min(1.0, Model.O[:Q2] - (T*delta)/Model.params[:θ] + tol)
   # end
#end

BufferedIntegral(Model::HopfieldHL1RSB,X::TI) where TI <: IntegrationMethod = X

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

Delta(Model::HopfieldHL1RSB) = 1 - Model.params[:β] * (1 - Model.O[:Q2])
Deltatheta(Model::HopfieldHL1RSB) = 1 - Model.params[:β] * (1 - Model.O[:Q2]) - Model.params[:β] * Model.params[:θ] * (Model.O[:Q2] - Model.O[:Q1]) 

is_free_energy_definite(Model::HopfieldHL1RSB) = (Delta(Model)>0) && (Deltatheta(Model) >0) && (Model.O[:Q2] >= Model.O[:Q1])

function Deltas!(Model::HopfieldHL1RSB)
    mineps = 1e-30;
    Model.aux["Δ"] = max(Delta(Model),mineps)
    Model.aux["Δtheta"] = max(Deltatheta(Model),(1-Model.params[:θ])*mineps)
end

#### Conjugate order params SP
function exprOconj!(Model::HopfieldHL1RSB)
    @extract Model : O Oconj aux params
    Deltas!(Model)
    #P1
    ########### IMPORTANT ############# w.r.t. Linda's paper, I am removing a beta factor
    Oconj[:P1] = O[:Q1] / aux["Δtheta"]^2 
    #P2 
    Oconj[:P2] =  Oconj[:P1] +  ( Model.O[:Q2] - Model.O[:Q1] ) / (aux["Δ"]*aux["Δtheta"]) 
end

#### auxiliary functions and inner integrals
function Ξ(Model::HopfieldHL1RSB,z,t)
    @extract Model : params O Oconj
    return params[:β] * ( O[:M] + sqrt(params[:α]*(Oconj[:P2] - Oconj[:P1])) * z  + sqrt(params[:α]*Oconj[:P1])*t )
end

function K0(Model::HopfieldHL1RSB,X::TI,t) where TI <: IntegrationMethod
    I0(z) = cosh(Ξ(Model,z,t))^Model.params[:θ]
    return integrate(X,I0)
end

function K1(Model::HopfieldHL1RSB,X::TI,t) where TI <: IntegrationMethod
    I1(z) = tanh(Ξ(Model,z,t)) * cosh(Ξ(Model,z,t))^Model.params[:θ]
    return integrate(X,I1)
end

function K2(Model::HopfieldHL1RSB,X::TI,t) where TI <: IntegrationMethod
    I2(z) = tanh(Ξ(Model,z,t))^2 * cosh(Ξ(Model,z,t))^Model.params[:θ]
    return integrate(X,I2)
end

#### SP equation for M
function exprM(Model::HopfieldHL1RSB,X::TI) where TI <: IntegrationMethod
    Model.O[:M]==0 && return 0.0 # is that correct? probably yes
    IM(t) = K1(Model,X,t) / K0(Model,X,t)
    return integrate(X,IM)
end

### SP equation for Q1
function exprQ1(Model::HopfieldHL1RSB,X::TI) where TI <: IntegrationMethod
    # limiting cases y=1 or gamma=0.0 to be yet handled
    IQ1(t) = ( K1(Model,X,t) / K0(Model,X,t) ) ^2
    return integrate(X,IQ1)
end

### SP equation for Q2
function exprQ2(Model::HopfieldHL1RSB,X::TI) where TI <: IntegrationMethod
    IQ2(t) = K2(Model,X,t) / K0(Model,X,t)
    return integrate(X,IQ2)
end

function lhs!(Model::HopfieldHL1RSB,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconj!(Model)
    Onew[:M] = exprM(Model,X)
    mineps = 1e-50
    Onew[:Q1] = clamp(exprQ1(Model,X), mineps, Inf)
    Onew[:Q2] = clamp(exprQ2(Model,X), Onew[:Q1] + mineps, Inf)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::HopfieldHL1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    exprOconj!(Model)
    energy = 0.5*O[:M]^2 
    energy_overlaps1 = - ( 0.5 * params[:α] * O[:Q1]) / aux["Δtheta"]
    energy_overlaps2 = 0.5 * params[:α] * params[:β] * Oconj[:P2] * (1 - O[:Q2] )
    energy_overlaps3 = 0.5 * params[:α] * params[:β] * params[:θ] * ( Oconj[:P2]*O[:Q2] - Oconj[:P1]*O[:Q1] )

    energy_overlaps = energy_overlaps1 + energy_overlaps2 + energy_overlaps3
    
    #logterms =  - ( ( 0.5 * params[:α]) / (params[:β] * params[:θ]) ) * log(1 + params[:β]*params[:θ]*(O[:Q2] - O[:Q1]) / aux["Δtheta"] ) + ( ( 0.5 * params[:α]) / (params[:β] ) ) * log(aux["Δ"]) 
    logterms = ( ( 0.5 * params[:α]) / (params[:β]*params[:θ]) )*log(aux["Δtheta"]) + ( ( 0.5 * params[:α]) / (params[:β] ) ) *(1 - 1/params[:θ])* log(aux["Δ"]) 
    
    IS(t) = log(K0(Model, X,t))
    return  energy + energy_overlaps + logterms - integrate(X,IS)/(params[:β]*params[:θ])
    
end