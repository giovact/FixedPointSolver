struct BinaryPerceptron1RSBDyntheta1withWick{PT<:Potential} <: Perceptron
    params::NamedTuple{(:β,:α),Tuple{FT,FT}}
    V::PT
    O::NamedVec #vector of orderparameters -> M, Q
    Oconj::NamedVec #vector of conjugate orderparameters -> Qhat
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
BinaryPerceptron1RSBDyntheta1withWick(rule::PT) where PT<:Potential = BinaryPerceptron1RSBDyntheta1withWick{PT}()
BinaryPerceptron1RSBDyntheta1withWick(β,α,rule::PT) where PT<:Potential = BinaryPerceptron1RSBDyntheta1withWick(β,α,rule,zeros(1))
BinaryPerceptron1RSBDyntheta1withWick(β,α,rule::PT,O) where PT<:Potential = BinaryPerceptron1RSBDyntheta1withWick((β=β,α=α),rule,NamedArray(O,[:Q0]),NamedArray(zeros(1),[:Q0]), Dict{String, FT}())

create_model(Model::BinaryPerceptron1RSBDyntheta1withWick, params) = BinaryPerceptron1RSBDyntheta1withWick(params, Model.V,Model.O, Model.Oconj, Dict{String, FT}())

set_integrationmethod(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI) where TI <: IntegrationMethod = X

function enforce_initial_condition!(Model::BinaryPerceptron1RSBDyntheta1withWick)
    XRS = LegendreQuadrature(10000;bound=8.0)
    ModelRS = BinaryPerceptronRS(Model.params[:β],Model.params[:α], Model.V, [1e-1; 1e-1])
    conv, iter, vareps = solve_SPequations!(ModelRS, XRS; tol = 1e-10, K = FixedPoint(0.9), niter = 5000,printprogress=false,showinfo =false)
    Model.aux["R"] = ModelRS.O[:R]
    Model.aux["Rhat"] = ModelRS.Oconj[:R]

    Model.aux["Q1"] = ModelRS.O[:Q]
    Model.aux["Q1hat"] = ModelRS.Oconj[:Q]

#    @info "Preprocessing RS -> (R,Q)=($(Model.aux["R"]), $(Model.aux["Q1"])) -> conv=$conv in $iter iterations "

    #Model.O[:Q0] = min(1.0-1e-5, Model.aux["Q1"] + 0.5* (1 - Model.aux["Q1"]))
    Model.O[:Q0] = 1.0-1e-3
    #Model.O[:Q0] = min(1.0-1e-5, Model.aux["Q1"] + 1e-2)
end


########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function aux_variables!(Model::BinaryPerceptron1RSBDyntheta1withWick)
    @extract Model : O aux params
    
    mineps = 1e-50

    aux["sqrt_Q1mR2"] = sqrt(clamp(aux["Q1"] - aux["R"]^2, mineps, Inf)) 
    aux["sqrt_1mQ0"] = sqrt(clamp((1 - O[:Q0]), mineps, Inf))
    
    aux["a"] = aux["R"] / aux["sqrt_Q1mR2"]
    aux["a2"] = aux["R"]^2 / (aux["Q1"]-aux["R"]^2)
    
    aux["b2"] = aux["Q1"] / (1-O[:Q0])
    aux["b"] = sqrt(clamp(aux["b2"], mineps, Inf))

    aux["c2"] = (O[:Q0] - aux["Q1"]) / (1 - O[:Q0])
    aux["c"] = sqrt(clamp(aux["c2"], mineps, Inf))

    aux["sqrtQ0mQ1"] =  sqrt( clamp( (O[:Q0] - aux["Q1"])  , mineps, Inf) ) 

    aux["betafactor"] = 1-exp(-params[:β]) # only for Gibbs

    aux["ce"] = 0.5 * params[:β] / (sqrt(2π) * aux["sqrt_1mQ0"] )
end


#### Conjugate order params SP

function exprOconjall!(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI) where TI <: IntegrationMethod
    @extract Model : Oconj aux
    aux_variables!(Model)
    Oconj[:Q0] = clamp(exprQ0hat(Model, X),  aux["Q1"] + 1e-40, Inf)
end

# order params SP 
function K0ent(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI,τ) where TI <: IntegrationMethod
    I0ent(z) = cosh(χI(Model,z,τ))
    return integrate(X,I0ent)
end

function KTanhSinh(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI,τ) where TI <: IntegrationMethod
    I1ent(z) = tanh(χI(Model,z,τ)) * sinh(χI(Model,z,τ))
    return integrate(X,I1ent)
end

### SP equation for Q0
function exprQ0(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI) where TI <: IntegrationMethod
    IQ0(τ) = KTanhSinh(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IQ0)
end

function lhs!(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconjall!(Model,X)
    Onew[:Q0] = exprQ0(Model, X)
end

########################################################## FREE ENERGY #########################################################
function free_energy(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI) where TI <: IntegrationMethod
    return 1.0
end



#######################################################################################################################################################################
############################################################################# 1 step RSB functions ####################################################################
#######################################################################################################################################################################


# dummy constructor for testing
BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron}() = BinaryPerceptron1RSBDyntheta1withWick(1.0,1.0,TSPerceptron())

# Ξenerg
function Rexp(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return 2 * aux["ce"] * exp(-0.5 * (aux["b"]*t + aux["c"]*τ)^2 ) - params[:β]^2 * Texp(Model,τ,t) * Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ)
end

function Texp(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return exp(0.5*(1-O[:Q0])*params[:β]^2 - params[:β]*(sqrt(aux["Q1"])*t + aux["sqrtQ0mQ1"]*τ))
end

function Ξenerg(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return Herf(aux["b"]*t + aux["c"]*τ) + Texp(Model,τ,t)*Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ)
end



function K0energ(Model::BinaryPerceptron1RSBDyntheta1withWick,X::TI,t) where TI <: IntegrationMethod
    I0energ(τ) = Ξenerg(Model,τ,t)
    return integrate(X,I0energ)
end

function Kenerg_Rexp(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    IKenerg_Rexp(τ) = Rexp(Model,τ,t)
    return integrate(X,IKenerg_Rexp)
end

function Kenerg_Rexp_logΞ(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    IKenerg_Rexp_logΞ(τ) = Rexp(Model,τ,t) * log(Ξenerg(Model,τ,t))
    return integrate(X,IKenerg_Rexp_logΞ)
end

function Kenerg_ΞlogΞ(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    IKenerg_ΞlogΞ(τ) = Ξenerg(Model,τ,t) * log(Ξenerg(Model,τ,t))
    return integrate(X,IKenerg_ΞlogΞ)
end


function Kenerg_ΞinvTphi2(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux 
    Ienerg_ΞinvTphi2(τ) = ( Texp(Model, τ,t) * Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ))^2 / Ξenerg(Model,τ,t)
    return integrate(X,Ienerg_ΞinvTphi2)
end


function exprQ0hat(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O
    term_1(t) = Kenerg_Rexp_logΞ(Model, X, t) / K0energ(Model, X, t)
    term_2(t) = Kenerg_ΞlogΞ(Model, X, t) / K0energ(Model, X, t)
    term_3(t) = Kenerg_Rexp(Model, X, t) / K0energ(Model, X, t)
    term_4(t) = Kenerg_ΞinvTphi2(Model, X, t) / K0energ(Model, X, t)
    
    IQ0hat(t) = Herf(aux["a"]*t) * ( term_1(t) - term_2(t)*term_3(t) - 0.5*params[:β]^2 * term_4(t)) 
    return 4*params[:α]*integrate(X,IQ0hat)
end

χI(Model::BinaryPerceptron1RSBDyntheta1withWick{TSPerceptron},z,τ) = Model.aux["Rhat"] + sqrt(clamp(Model.Oconj[:Q0] - Model.aux["Q1hat"], 1e-20,Inf))*z  + sqrt(Model.aux["Q1hat"])*τ