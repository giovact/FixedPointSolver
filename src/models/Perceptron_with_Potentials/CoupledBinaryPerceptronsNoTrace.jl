struct CoupledBinaryPerceptronsNoTrace{PT<:Potential} <: Perceptron
    params::NamedTuple{(:β,:α,:γ, :y),Tuple{FT,FT,FT, Int}}
    V::PT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
CoupledBinaryPerceptronsNoTrace(rule::PT) where PT<:Potential = CoupledBinaryPerceptronsNoTrace{PT}()
CoupledBinaryPerceptronsNoTrace(β,α,γ,y,rule::PT) where PT<:Potential = CoupledBinaryPerceptronsNoTrace(β,α,γ,y,rule,zeros(3))
CoupledBinaryPerceptronsNoTrace(β,α,γ,y,rule::PT,O) where PT<:Potential = CoupledBinaryPerceptronsNoTrace((β=β,α=α,γ=γ,y=y),rule,NamedArray(O,[:R;:Q;:P]),NamedArray(zeros(3),[:R;:Q;:P]),Dict{String, FT}())

create_model(Model::CoupledBinaryPerceptronsNoTrace, params) = CoupledBinaryPerceptronsNoTrace(params, Model.V,Model.O, Model.Oconj, Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::CoupledBinaryPerceptronsNoTrace)
    @extract Model : O aux params
    
    mineps = 1e-50

    aux["sqrt_QmR2"] = sqrt(clamp(O[:Q] - O[:R]^2, mineps, Inf)) 
    aux["sqrt_1mP"] = sqrt(clamp((1 - O[:P]), mineps, Inf))
    
    aux["a"] = O[:R] / aux["sqrt_QmR2"]
    aux["a2"] = O[:R]^2 / (O[:Q]-O[:R]^2)
    
    aux["b2"] = O[:Q] / (1-O[:P])
    aux["b"] = sqrt(clamp(aux["b2"], mineps, Inf))

    aux["c2"] = (O[:P] - O[:Q]) / (1 - O[:P])
    aux["c"] = sqrt(clamp(aux["c2"], mineps, Inf))

    aux["sqrtPmQ"] =  sqrt( clamp( (O[:P] - O[:Q])  , mineps, Inf) ) 

    aux["betafactor"] = 1-exp(-params[:β]) # only for Gibbs
end

#### Conjugate order params SP
function K0energ(Model::CoupledBinaryPerceptronsNoTrace,X::TI,t) where TI <: IntegrationMethod
    I0energ(τ) = Ξenerg(Model,τ,t)^Model.params[:y]
    return integrate(X,I0energ)
end

function exprOconjall!(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj
    aux_variables!(Model)
    Oconj[:R] = exprRhat(Model, X)
    Oconj[:Q] = exprQhat(Model, X)
    Oconj[:P] = clamp(exprPhat(Model, X),  Oconj[:Q] + 1e-40, Inf)
end

# ORDER PARAMS SP #########
function K0ent(Model::CoupledBinaryPerceptronsNoTrace,X::TI,τ) where TI <: IntegrationMethod
    I0ent(z) = (cosh(χI(Model,z,τ)))^Model.params[:y]
    return integrate(X,I0ent)
end

function K1ent(Model::CoupledBinaryPerceptronsNoTrace,X::TI,τ) where TI <: IntegrationMethod
    I1ent(z) = tanh(χI(Model,z,τ)) * cosh(χI(Model,z,τ))^Model.params[:y]
    return integrate(X,I1ent)
end

function K2ent(Model::CoupledBinaryPerceptronsNoTrace,X::TI,τ) where TI <: IntegrationMethod
    I2ent(z) = tanh(χI(Model,z,τ))^2 * cosh(χI(Model,z,τ))^Model.params[:y]
    return integrate(X,I2ent)
end

#### SP equation for R
function exprR(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    IR(τ) = K1ent(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IR)
end

### SP equation for Q1
function exprQ(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    IQ1(τ) = ( K1ent(Model,X,τ) / K0ent(Model,X,τ) ) ^2
    return integrate(X,IQ1)
end

### SP equation for Q0
function exprP(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    IQ0(τ) = K2ent(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IQ0)
end

function lhs!(Model::CoupledBinaryPerceptronsNoTrace,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconjall!(Model,X)
    Onew[:R] = exprR(Model, X)
    Onew[:Q] = exprQ(Model, X)
    Onew[:P] = exprP(Model, X)
end

########################################################## FREE ENERGY #########################################################
function GE_in_free_energy(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    GEinner(t) = Herf(aux["a"]*t)*log(K0energ(Model,X,t))
    return (2 / params[:y] )*integrate(X,GEinner)
end

function GI_in_free_energy(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    GIinner(τ) = log(K0ent(Model, X,τ)) 
    return integrate(X,GIinner) / Model.params[:y] + log(2)
end

function free_energy(Model::CoupledBinaryPerceptronsNoTrace,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    exprOconjall!(Model,X)

    term_fromHS = - 0.5 * ((params[:γ] / params[:y]) + (Oconj[:P] - Oconj[:Q] ) ) # this is simply -Jtilde / 2
    #constraints = - O[:R]*Oconj[:R] + 0.5 * Oconj[:Q] * (O[:Q] - 1) - 0.5*(params[:y] - 1)*(O[:P]*Oconj[:P] - O[:Q]*Oconj[:Q])
    constraints = - O[:R]*Oconj[:R] - 0.5*(params[:y] - 1)*O[:P]*Oconj[:P] + 0.5 * Oconj[:Q] * (params[:y] * O[:Q] - 1)
    G = constraints + term_fromHS + GI_in_free_energy(Model,X) + params[:α] * GE_in_free_energy(Model,X)
    return - G # this quantity is beta*f
end

function compute_aux!(Model::CoupledBinaryPerceptronsNoTrace,X::TI)  where TI <: IntegrationMethod
    Model.aux["U_V"] = energy_V(Model,X)
    Model.aux["U_gamma"] = - 0.5 * Model.params[:γ] * ( (Model.params[:y]-1) / Model.params[:y] ) * Model.O[:P] # this is NOT divided by beta
    Model.aux["energy"] = Model.aux["U_V"] - Model.aux["U_gamma"]/Model.params[:β]
    Model.aux["entropy"] = Model.params[:β]*Model.aux["U_V"] + Model.aux["U_gamma"] - Model.aux["f"]
end
