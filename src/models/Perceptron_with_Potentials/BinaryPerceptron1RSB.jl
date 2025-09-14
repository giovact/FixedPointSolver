struct BinaryPerceptron1RSB{PT<:Potential} <: Perceptron
    params::NamedTuple{(:β,:α,:θ),Tuple{FT,FT,FT}}
    V::PT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BinaryPerceptron1RSB(rule::PT) where PT<:Potential = BinaryPerceptron1RSB{PT}()
BinaryPerceptron1RSB(β,α,θ,rule::PT) where PT<:Potential = BinaryPerceptron1RSB(β,α,θ,rule,zeros(3))
BinaryPerceptron1RSB(β,α,θ,rule::PT,O) where PT<:Potential = BinaryPerceptron1RSB((β=β,α=α,θ=θ),rule,NamedArray(O,[:R;:Q1;:Q0]),NamedArray(zeros(3),[:R;:Q1;:Q0]), Dict{String, FT}())

create_model(Model::BinaryPerceptron1RSB, params) = BinaryPerceptron1RSB(params, Model.V,Model.O, Model.Oconj, Dict{String, FT}())

set_integrationmethod(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod = X

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::BinaryPerceptron1RSB)
    @extract Model : O aux params
    
    mineps = 1e-50

    aux["sqrt_Q1mR2"] = sqrt(clamp(O[:Q1] - O[:R]^2, mineps, Inf)) 
    aux["sqrt_1mQ0"] = sqrt(clamp((1 - O[:Q0]), mineps, Inf))
    
    aux["a"] = O[:R] / aux["sqrt_Q1mR2"]
    aux["a2"] = O[:R]^2 / (O[:Q1]-O[:R]^2)
    
    aux["b2"] = O[:Q1] / (1-O[:Q0])
    aux["b"] = sqrt(clamp(aux["b2"], mineps, Inf))

    aux["c2"] = (O[:Q0] - O[:Q1]) / (1 - O[:Q0])
    aux["c"] = sqrt(clamp(aux["c2"], mineps, Inf))

    aux["sqrtQ0mQ1"] =  sqrt( clamp( (O[:Q0] - O[:Q1])  , mineps, Inf) ) 

    aux["betafactor"] = 1-exp(-params[:β]) # only for Gibbs
end


#### Conjugate order params SP
function K0energ(Model::BinaryPerceptron1RSB,X::TI,t) where TI <: IntegrationMethod
    I0energ(τ) = Ξenerg(Model,τ,t)^Model.params[:θ]
    return integrate(X,I0energ)
end

function exprOconjall!(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj
    aux_variables!(Model)
    Oconj[:R] = exprRhat(Model, X)
    Oconj[:Q1] = exprQ1hat(Model, X)
    Oconj[:Q0] = clamp(exprQ0hat(Model, X),  Oconj[:Q1] + 1e-40, Inf)
end

# order params SP 
function K0ent(Model::BinaryPerceptron1RSB,X::TI,τ) where TI <: IntegrationMethod
    I0ent(z) = (cosh(χI(Model,z,τ)))^Model.params[:θ]
    return integrate(X,I0ent)
end

function K1ent(Model::BinaryPerceptron1RSB,X::TI,τ) where TI <: IntegrationMethod
    I1ent(z) = tanh(χI(Model,z,τ)) * cosh(χI(Model,z,τ))^Model.params[:θ]
    return integrate(X,I1ent)
end

function K2ent(Model::BinaryPerceptron1RSB,X::TI,τ) where TI <: IntegrationMethod
    I2ent(z) = tanh(χI(Model,z,τ))^2 * cosh(χI(Model,z,τ))^Model.params[:θ]
    return integrate(X,I2ent)
end

#### SP equation for R
function exprR(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    IR(τ) = K1ent(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IR)
end

### SP equation for Q1
function exprQ1(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    IQ1(τ) = ( K1ent(Model,X,τ) / K0ent(Model,X,τ) ) ^2
    return integrate(X,IQ1)
end

### SP equation for Q0
function exprQ0(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    IQ0(τ) = K2ent(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IQ0)
end

function lhs!(Model::BinaryPerceptron1RSB,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconjall!(Model,X)
    Onew[:R] = exprR(Model, X)
    Onew[:Q1] = exprQ1(Model, X)
    Onew[:Q0] = exprQ0(Model, X)
end

########################################################## FREE ENERGY #########################################################
function GE_in_free_energy(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    GEinner(t) = Herf(aux["a"]*t)*log(K0energ(Model,X,t)) 
    return (2 / params[:θ] )*integrate(X,GEinner)
end

function GI_in_free_energy(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    GIinner(τ) = log(K0ent(Model, X,τ))
    return integrate(X,GIinner) / Model.params[:θ] + log(2)
end

function free_energy(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    
    exprOconjall!(Model,X)
    constraints = - O[:R]*Oconj[:R] + 0.5 * Oconj[:Q0] * (O[:Q0] - 1) - 0.5*params[:θ]*(O[:Q0]*Oconj[:Q0] - O[:Q1]*Oconj[:Q1])

    return - (constraints + GI_in_free_energy(Model,X) + params[:α] * GE_in_free_energy(Model,X)) # this quantity is beta*f
end

########################################################## COMPLEXITY #########################################################

function Klog_ent(Model::BinaryPerceptron1RSB,X::TI,τ) where TI <: IntegrationMethod
    I0ent_logterm(z) = logcosh(χI(Model,z,τ)) * (cosh(χI(Model,z,τ)))^Model.params[:θ]
    return integrate(X,I0ent_logterm) 
end

function Klog_energ(Model::BinaryPerceptron1RSB,X::TI,t) where TI <: IntegrationMethod 
    @extract Model : params O Oconj aux
    I2energ(τ) = log(Ξenerg(Model,τ,t)) * Ξenerg(Model,τ,t)^Model.params[:θ]
    return integrate(X,I2energ)
end

function complexity(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    exprOconjall!(Model,X)
    
    # first part, constraints
    Σ_constraints = 0.5 * (params[:θ]^2) *(O[:Q0]*Oconj[:Q0] - O[:Q1]*Oconj[:Q1])

    # second part, entropy
    GI_outerlog_int(τ) = log(K0ent(Model, X,τ))
    GI_innerlog_int(τ) = Klog_ent(Model, X,τ) / K0ent(Model, X,τ)

    ΣI = integrate(X,GI_outerlog_int) - params[:θ] * integrate(X,GI_innerlog_int)

    # third part, energy
    GE_outerlog_int(t) = Herf(aux["a"]*t) * log( K0energ(Model, X, t) )
    GE_innerlog_int(t) = Herf(aux["a"]*t) * Klog_energ(Model, X, t) / K0energ(Model, X, t)

    ΣE = 2*params[:α] * ( integrate(X,GE_outerlog_int) - params[:θ] * integrate(X,GE_innerlog_int) ) 

    return Σ_constraints + ΣI + ΣE
end

function compute_aux!(Model::BinaryPerceptron1RSB,X::TI) where TI <: IntegrationMethod
    Model.aux["energy"] = energy_V(Model,X)
    Model.aux["entropy"] = Model.params[:β]*Model.aux["energy"] - Model.aux["f"]
    Model.aux["complexity"] = complexity(Model, X)
end

