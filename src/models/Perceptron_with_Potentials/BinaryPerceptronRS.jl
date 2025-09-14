struct BinaryPerceptronRS{PT<:Potential} <: Perceptron
    params::NamedTuple{(:β,:α),Tuple{FT,FT}}
    V::PT
    O::NamedVec 
    Oconj::NamedVec 
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BinaryPerceptronRS(rule::PT) where PT<:Potential = BinaryPerceptronRS{PT}()
BinaryPerceptronRS(β,α,rule::PT) where PT<:Potential = BinaryPerceptronRS(β,α,rule,zeros(2))
BinaryPerceptronRS(β,α,rule::PT,O) where PT<:Potential = BinaryPerceptronRS((β=β,α=α),rule,NamedArray(O,[:R;:Q]),NamedArray(zeros(2),[:R;:Q]), Dict{String, FT}())

create_model(Model::BinaryPerceptronRS, params) = BinaryPerceptronRS(params, Model.V,Model.O, Model.Oconj, Dict{String, FT}())


########################################################## SELF-CONSISTENT EQUATIONS #########################################################
function aux_variables!(Model::BinaryPerceptronRS)
    @extract Model : O aux params
    
    mineps = 1e-50

    aux["sqrt_QmR2"] = sqrt(clamp(O[:Q] - O[:R]^2, mineps, Inf)) 
    aux["sqrt_1mQ"] = sqrt(clamp((1 - O[:Q]), mineps, Inf))
    aux["a"] = O[:R] / aux["sqrt_QmR2"] 

    aux["b2"] = clamp(O[:Q] / (1-O[:Q]), mineps/100, Inf)
    aux["b"] = sqrt(clamp(aux["b2"], mineps, Inf))

    aux["v"]  = sqrt(clamp((O[:Q] - O[:R]^2) / (1 - O[:Q]), mineps, Inf))
    aux_potential!(Model)
end

function exprOconjall!(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    aux_variables!(Model)
    Oconj[:R] = Model.V.storage == true ? 0.0 : exprRhat(Model, X)
    Oconj[:Q] = exprQhat(Model, X)
end

# order params SP
function exprR(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod
    IR(z) = tanh(χI(Model,z))
    return integrate(X,IR)
end

function exprQ(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj
    IQ(z) = tanh( χI(Model,z) )^2
    return integrate(X,IQ)
end

function lhs!(Model::BinaryPerceptronRS,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconjall!(Model,X)
    Onew[:R] = exprR(Model,X)
    Onew[:Q] = exprQ(Model,X)

    if Model.V.storage==true
        Onew[:R]=0.0
    end
end

########################################################## FREE ENERGY #########################################################

function GE_in_free_energy(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    IE(t) = Herf(aux["a"]*t)*log( KV(Model,t) ) 
    return 2*integrate(X,IE)
end

function free_energy(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod 
    @extract Model : params O Oconj
    exprOconjall!(Model,X)
    constraints = 0.5 * Oconj[:Q] * (O[:Q] - 1) - O[:R]*Oconj[:R]

    IS(z) = log2cosh( χI(Model,z) )
    GI = integrate(X,IS)
    
    G = constraints + GI + params[:α]*GE_in_free_energy(Model, X)
    
    return - G # this quantity is beta*f
end



function compute_aux!(Model::BinaryPerceptronRS,X::TI)  where TI <: IntegrationMethod
    Model.aux["energy"] = energy(Model,X)
    Model.aux["entropy"] = entropy(Model,X)

    Model.aux["λE"] = λE(Model,X)
    Model.aux["λI"] = λI(Model,X)
    Model.aux["dAT"] = Model.aux["λE"] * Model.aux["λI"] - 1
    
    Model.aux["Fannealed"] = annealed_free_energy(Model, X)
end

# function for DAT line
function λI(Model::BinaryPerceptronRS,X::TI) where TI <: IntegrationMethod
    @extract Model : params aux O Oconj
    λI_inner(z) = cosh(χI(Model,z))^(-4)
    return integrate(X,λI_inner)
end
