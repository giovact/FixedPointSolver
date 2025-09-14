
#######################################################################################################################################################################
############################################################################# 1 step RSB functions ####################################################################
#######################################################################################################################################################################


# dummy constructor for testing
BinaryPerceptron1RSB{TSPerceptron}() = BinaryPerceptron1RSB(1.0,1.0,1.0,TSPerceptron(),[1e-1; 1e-1;2e-1] )

# Ξenerg
function Texp(Model::BinaryPerceptron1RSB{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return exp(0.5*(1-O[:Q0])*params[:β]^2 - params[:β]*(sqrt(O[:Q1])*t + aux["sqrtQ0mQ1"]*τ))

end
function Ξenerg(Model::BinaryPerceptron1RSB{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return Herf(aux["b"]*t + aux["c"]*τ) + Texp(Model,τ,t)*Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ)
end

function K1energ(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I1energ(τ) = Texp(Model,τ,t)*Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ)*Ξenerg(Model,τ,t)^(Model.params[:θ]-1)
    return integrate(X,I1energ)
end

function K2energ(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I2energ(τ) = ( Texp(Model,τ,t)*Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ) ) ^2 *Ξenerg(Model,τ,t)^(Model.params[:θ]-2)
    return integrate(X,I2energ)
end

function exprRhat(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params O aux 
    IRhat(t) = exp(-0.5*aux["a2"]*t^2) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) ) 
    return sqrt(2/π)*params[:α]*params[:β]*(sqrt(O[:Q1]) / aux["sqrt_Q1mR2"] ) *integrate(X,IRhat)
end

function exprQ1hat(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O 
    IQ1hat(t) = Herf(aux["a"]*t) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) )^2
    return 2*params[:α]*params[:β]^2 * integrate(X,IQ1hat)
end

function exprQ0hat(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O
    IQ0hat(t) = Herf(aux["a"]*t) * ( K2energ(Model, X,t) /  K0energ(Model, X, t) )
    return 2*params[:α]*params[:β]^2*integrate(X,IQ0hat)
end


χI(Model::BinaryPerceptron1RSB{TSPerceptron},z,τ) = Model.Oconj[:R] + sqrt(Model.Oconj[:Q0] - Model.Oconj[:Q1])*z  + sqrt(Model.Oconj[:Q1])*τ



# energy
function Kenergy_V(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux
    IE_inner(τ) =  ( (exp(-0.5*(aux["b"]*t + aux["c"]*τ)^2) /sqrt(2*π)) + (aux["b"]*t + aux["c"]*τ - params[:β]*aux["sqrt_1mQ0"])*Texp(Model,τ,t)*Herf(sqrt(1-O[:Q0])*params[:β] - aux["b"]*t - aux["c"]*τ) ) * Ξenerg(Model,τ,t)^(Model.params[:θ]-1)
    return integrate(X,IE_inner)
end

function energy_V(Model::BinaryPerceptron1RSB{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    Ienergy(t) = Herf(aux["a"]*t) * ( Kenergy_V(Model,X,t) /  K0energ(Model,X,t))
    return 2*params[:α]*aux["sqrt_1mQ0"]*integrate(X,Ienergy)
end