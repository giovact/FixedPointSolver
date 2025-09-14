#######################################################################################################################################################################
############################################################# Replicated Case with ferromagnetic coupling #############################################################
#######################################################################################################################################################################


# dummy constructor for testing
CoupledBinaryPerceptronsNoTrace{TSGibbs}() = CoupledBinaryPerceptronsNoTrace(1.0,1.0,1.0,2,TSGibbs(),[1e-1; 1e-1;2e-1] )


function Ξenerg(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},τ,t)
    @extract Model : params O Oconj aux 
    return Herf(aux["b"]*t + aux["c"]*τ) + exp(-params[:β]) * Herf( - aux["b"]*t - aux["c"]*τ)
end

function K1energ(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I1energ(τ) = exp(- 0.5 *(Model.aux["b"]*t + Model.aux["c"]*τ)^2 ) * Ξenerg(Model,τ,t)^(Model.params[:y]-1)
    return integrate(X,I1energ)
end

function K2energ(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux 
    I2energ(τ) = exp(- (Model.aux["b"]*t + Model.aux["c"]*τ)^2 ) * Ξenerg(Model,τ,t)^(Model.params[:y]-2)
    return integrate(X,I2energ)
end

function exprRhat(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI) where TI <: IntegrationMethod
    Model.V.storage == true && return 0.0
    @extract Model : params O aux 
    IRhat(t) = exp(-0.5*aux["a"]^2*t^2) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) ) 
    return (params[:α] / π )  * (aux["b"] / aux["sqrt_QmR2"]) * aux["betafactor"] *integrate(X,IRhat)
end

function exprQhat(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O 
    IQhat(t) = Herf(aux["a"]*t) * ( K1energ(Model, X,t) /  K0energ(Model, X, t) )^2
    return (params[:α] / π ) * (  aux["betafactor"]^2 / (1 - O[:P]) ) * integrate(X,IQhat)
end

function exprPhat(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : aux params O
    IPhat(t) = Herf(aux["a"]*t) * ( K2energ(Model, X,t) /  K0energ(Model, X, t) )
    return (params[:α] / π ) * (  aux["betafactor"]^2 / (1 - O[:P]) ) *integrate(X,IPhat)
end


χI(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},z,τ) = Model.Oconj[:R] + sqrt((Model.params[:γ] /Model.params[:y]) +  Model.Oconj[:P] - Model.Oconj[:Q])*z + sqrt(Model.Oconj[:Q])*τ


function Kenergy_V(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI,t) where TI <: IntegrationMethod
    @extract Model : params O aux
    IE_inner(τ) =  Herf( - aux["b"]*t - aux["c"]*τ) * Ξenerg(Model,τ,t)^(Model.params[:y]-1)
    return integrate(X,IE_inner)
end

function energy_V(Model::CoupledBinaryPerceptronsNoTrace{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    params[:β]==Inf && return 0.0
    Ienergy(t) = Herf(aux["a"]*t) * ( Kenergy_V(Model,X,t) /  K0energ(Model,X,t))
    return 2*params[:α]* exp(-params[:β])* integrate(X,Ienergy)
end
