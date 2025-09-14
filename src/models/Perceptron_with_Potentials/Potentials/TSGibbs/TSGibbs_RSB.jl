#######################################################################################################################################################################
############################################################################# 1 step RSB functions ####################################################################
#######################################################################################################################################################################

# dummy constructor for testing
BinaryPerceptron1RSB{TSGibbs}() = BinaryPerceptron1RSB(1.0,1.0,1.0,TSGibbs(),[1e-1; 1e-1;2e-1] )


# Ξenerg
function Ξenerg(Model::BinaryPerceptron1RSB{TSGibbs},τ,t)
    @extract Model : params  aux 
    # consider to write it in the following way
    ## exp(-params[:β]) + (1-exp(-params[:β]))*Herf(aux["b"]*t + aux["c"]*τ)
    return Herf(aux["b"]*t + aux["c"]*τ) + exp(-params[:β]) * Herf( - aux["b"]*t - aux["c"]*τ)
end

function K1energ(Model::BinaryPerceptron1RSB{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    I1energ(τ) = exp(- 0.5 *(Model.aux["b"]*t + Model.aux["c"]*τ)^2 ) * Ξenerg(Model,τ,t)^(Model.params[:θ]-1)
    return integrate(X,I1energ)
end

function K2energ(Model::BinaryPerceptron1RSB{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    I2energ(τ) = exp(- (Model.aux["b"]*t + Model.aux["c"]*τ)^2 ) * Ξenerg(Model,τ,t)^(Model.params[:θ]-2)
    return integrate(X,I2energ)
end

function exprRhat(Model::BinaryPerceptron1RSB{TSGibbs},X::TI) where TI <: IntegrationMethod
    Model.V.storage == true && return 0.0
    @extract Model : params aux 
    IRhat(t) = exp(-0.5*aux["a"]^2*t^2) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) ) 
    return (params[:α] / π )  * (aux["b"] / aux["sqrt_Q1mR2"]) * aux["betafactor"] *integrate(X,IRhat)
end

function exprQ1hat(Model::BinaryPerceptron1RSB{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O 
    IQ1hat(t) = Herf(aux["a"]*t) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) )^2
    return (params[:α] / π ) * (  aux["betafactor"]^2 / (1 - O[:Q0]) ) * integrate(X,IQ1hat)
end

function exprQ0hat(Model::BinaryPerceptron1RSB{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O
    IQ0hat(t) = Herf(aux["a"]*t) * ( K2energ(Model, X,t) /  K0energ(Model, X, t) )
    return (params[:α] / π ) * (  aux["betafactor"]^2 / (1 - O[:Q0]) ) *integrate(X,IQ0hat)
end


χI(Model::BinaryPerceptron1RSB{TSGibbs},z,τ) = Model.Oconj[:R]*(1-Model.V.storage) + Model.aux["sqrtQ0mQ1"]*z  + sqrt(Model.Oconj[:Q1])*τ

function exprR(Model::BinaryPerceptron1RSB{TSGibbs},X::TI) where TI <: IntegrationMethod
    Model.V.storage == true && return 0.0
    IR(τ) = K1ent(Model,X,τ) / K0ent(Model,X,τ)
    return integrate(X,IR)
end

#energy
function Kenergy_V(Model::BinaryPerceptron1RSB{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux
    IE_inner(τ) =  Herf( - aux["b"]*t - aux["c"]*τ) * Ξenerg(Model,τ,t)^(Model.params[:θ]-1)
    return integrate(X,IE_inner)
end

function energy_V(Model::BinaryPerceptron1RSB{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    params[:β]==Inf && return 0.0
    Ienergy(t) = Herf(aux["a"]*t) * ( Kenergy_V(Model,X,t) /  K0energ(Model,X,t))
    return 2*params[:α]* exp(-params[:β])* integrate(X,Ienergy)
end
