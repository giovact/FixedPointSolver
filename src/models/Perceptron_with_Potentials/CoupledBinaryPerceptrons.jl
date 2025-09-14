struct CoupledBinaryPerceptrons{PT<:Potential} <: Perceptron
    params::NamedTuple{(:β,:α,:γ, :y),Tuple{FT,FT,FT, Int}}
    V::PT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
CoupledBinaryPerceptrons(rule::PT) where PT<:Potential = CoupledBinaryPerceptrons{PT}()
CoupledBinaryPerceptrons(β,α,γ,y,rule::PT) where PT<:Potential = CoupledBinaryPerceptrons(β,α,γ,y,rule,zeros(3))
CoupledBinaryPerceptrons(β,α,γ,y,rule::PT,O) where PT<:Potential = CoupledBinaryPerceptrons((β=β,α=α,γ=γ,y=y),rule,NamedArray(O,[:R;:Q;:P]),NamedArray(zeros(3),[:R;:Q;:P]), HomogeneousSpin1Half(y),Dict{String, FT}())

create_model(Model::CoupledBinaryPerceptrons, params) =  CoupledBinaryPerceptrons(params, Model.V,Model.O, Model.Oconj,HomogeneousSpin1Half(params[:y]), Dict{String, FT}())

set_integrationmethod(Model::CoupledBinaryPerceptrons,X::TI) where TI <: IntegrationMethod = TracedSpin1HalfIntegral(Model,X)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::CoupledBinaryPerceptrons)
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
function K0energ(Model::CoupledBinaryPerceptrons,X::TI,t) where TI <: IntegrationMethod
    I0energ(τ) = Ξenerg(Model,τ,t)^Model.params[:y]
    return integrate(X,I0energ)
end

function exprOconjall!(Model::CoupledBinaryPerceptrons,X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj
    aux_variables!(Model)
    Oconj[:R] = Model.V.storage == true ? 0.0 : exprRhat(Model, X)
    Oconj[:Q] = exprQhat(Model, X)
    Oconj[:P] = clamp(exprPhat(Model, X),  Oconj[:Q] + 1e-40, Inf)
end

# ORDER PARAMS SP #########

function lhs!(Model::CoupledBinaryPerceptrons,X::TracedSpin1HalfIntegral{TI},Onew::NamedVec) where TI <: IntegrationMethod

    exprOconjall!(Model,X.X)
    inner_trace(Model, X)
    
    Onew[:R] = integrate(X.X,X.m)
    Onew[:P] = integrate(X.X,X.chi)
    Onew[:Q] = integrate(X.X,X.m .^ 2)

    if Model.V.storage==true
        Onew[:R]=0.0
    end
end

########################################################## FREE ENERGY #########################################################
function GE_in_free_energy(Model::CoupledBinaryPerceptrons,X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    GEinner(t) = Herf(aux["a"]*t)*log(K0energ(Model,X.X,t))
    return (2 / params[:y] )*integrate(X.X,GEinner)
end

function GI_in_free_energy(Model::CoupledBinaryPerceptrons,X::TI) where TI <: IntegrationMethod
    inner_trace(Model, X)
    return integrate(X.X, log.(X.Z) )/ Model.params[:y]

end

function free_energy(Model::CoupledBinaryPerceptrons,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    
    exprOconjall!(Model,X.X)
    
    constraints = - O[:R]*Oconj[:R] - 0.5*(params[:y] - 1)*O[:P]*Oconj[:P] + 0.5 * Oconj[:Q] * (params[:y] * O[:Q] - 1)

    return - (constraints + GI_in_free_energy(Model,X) + params[:α] * GE_in_free_energy(Model,X)) # this quantity is beta*f
end

function compute_aux!(Model::CoupledBinaryPerceptrons,X::TI)  where TI <: IntegrationMethod
    Model.aux["U_V"] = energy_V(Model,X.X)
    Model.aux["U_gamma"] = - 0.5 * Model.params[:γ] * ( (Model.params[:y]-1) / Model.params[:y] ) * Model.O[:P] # this is NOT divided by beta
    Model.aux["energy"] = Model.aux["U_V"] - Model.aux["U_gamma"]/Model.params[:β]
    Model.aux["entropy"] = Model.params[:β]*Model.aux["U_V"] + Model.aux["U_gamma"] - Model.aux["f"]
end
