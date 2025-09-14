#######################################################################################################################################################################
############################################################# Replicated Case with ferromagnetic coupling #############################################################
#######################################################################################################################################################################


# dummy constructor for testing
CoupledBinaryPerceptrons{TSPerceptron}() = CoupledBinaryPerceptrons(1.0,1.0,1.0,2,TSPerceptron(),[1e-1; 1e-1;2e-1] )

# Ξenerg
function Texp(Model::CoupledBinaryPerceptrons{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return exp(0.5*(1-O[:P])*params[:β]^2 - params[:β]*(sqrt(O[:Q])*t + aux["sqrtPmQ"]*τ))
end
function Ξenerg(Model::CoupledBinaryPerceptrons{TSPerceptron},τ,t)
    @extract Model : params O Oconj aux 
    return Herf(aux["b"]*t + aux["c"]*τ) + Texp(Model,τ,t)*Herf(sqrt(1-O[:P])*params[:β] - aux["b"]*t - aux["c"]*τ)
end

function K1energ(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I1energ(τ) = Texp(Model,τ,t)*Herf(sqrt(1-O[:P])*params[:β] - aux["b"]*t - aux["c"]*τ)*Ξenerg(Model,τ,t)^(Model.params[:y]-1)
    return integrate(X,I1energ)
end

function K2energ(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I2energ(τ) = ( Texp(Model,τ,t)*Herf(sqrt(1-O[:P])*params[:β] - aux["b"]*t - aux["c"]*τ) ) ^2 * Ξenerg(Model,τ,t)^(Model.params[:y]-2)
    return integrate(X,I2energ)
end

function exprRhat(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params O aux 
    IRhat(t) = exp(-0.5*aux["a"]^2*t^2) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) ) 
    return sqrt(2/π)*params[:α]*params[:β]*(sqrt(O[:Q]) / aux["sqrt_QmR2"] ) *integrate(X,IRhat)
end

function exprQhat(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O 
    IQhat(t) = Herf(aux["a"]*t) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) )^2
    return 2*params[:α]*params[:β]^2 * integrate(X,IQhat)
end

function exprPhat(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O
    IPhat(t) = Herf(aux["a"]*t) * ( K2energ(Model, X,t) /  K0energ(Model, X, t) )
    return 2*params[:α]*params[:β]^2*integrate(X,IPhat)
end

function inner_trace(Model::CoupledBinaryPerceptrons{TSPerceptron},z::FT)
    @extract Model : params O Oconj
    heff = Oconj[:R] + z*sqrt(Oconj[:Q])
    Jeff = (params[:γ] / params[:y]) + (Oconj[:P] - Oconj[:Q] ) 
    return trace_CurieWeiss(Model.InnerCW,heff,Jeff)
end

function Kenergy_V(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux
    IE_inner(τ) =  ( (exp(-0.5*(aux["b"]*t + aux["c"]*τ)^2) * aux["sqrt_1mP"] / sqrt(2*π)) + (sqrt(O[:Q])*t + aux["sqrtPmQ"]*τ - params[:β]*(O[:P]-O[:Q]))*Texp(Model,τ,t)*Herf(params[:β]* aux["sqrt_1mP"] - aux["b"]*t - aux["c"]*τ) ) * Ξenerg(Model,τ,t)^(Model.params[:y]-1)
    return integrate(X,IE_inner)
end

function energy_V(Model::CoupledBinaryPerceptrons{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    Ienergy(t) = Herf(aux["a"]*t) * ( Kenergy_V(Model,X,t) /  K0energ(Model,X,t))
    return 2*params[:α]*integrate(X,Ienergy)
end