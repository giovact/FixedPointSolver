# dummy constructor for testing
BinaryPerceptronRS{TSPerceptron}() = BinaryPerceptronRS(1.0,1.0,TSPerceptron(),ones(2)*1e-2 )

function aux_potential!(Model::BinaryPerceptronRS{TSPerceptron})
    Model.aux["jj"] = Model.params[:β]*sqrt(1 - Model.O[:Q])
end

# KV this is what you called K in the notes, at least for the RS
function KV(Model::BinaryPerceptronRS{TSPerceptron},t)
    @extract Model : params aux O
    Herf(aux["b"]*t) + exp(0.5*aux["jj"]^2 - params[:β]* sqrt(O[:Q]) *t )*Herf(aux["jj"] - aux["b"]*t)
end

# Equation for Rhat
function exprRhat(Model::BinaryPerceptronRS{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    #IconjR(z) = exp(0.5*jj^2 - params[:β]* aux["sqrt_QmR2"] *z )*Herf(jj - v*z) /  ( Herf(v*z) + exp(0.5*jj^2 - params[:β]* aux["sqrt_QmR2"] *z )*Herf(jj - v*z) )
    IconjR(z) = 1 /  ( 1 +  ( Herf(aux["v"]*z) / Herf(aux["jj"] - aux["v"]*z)) * exp(params[:β]* aux["sqrt_QmR2"] *z - 0.5*aux["jj"]^2  ) )
    return ( ( 2 * params[:α] * params[:β]) / sqrt(2*π) ) *integrate(X,IconjR)
end

# Equation for Qhat
function exprQhat(Model::BinaryPerceptronRS{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params aux  O 
    #IconjQ(t) = Herf(aux["a"]*t) * ( exp(0.5*aux["jj"]^2 - params[:β]*sqrt(O[:Q])*t )*Herf(aux["jj"] - aux["b"]*t) /  ( Herf(aux["b"]*t) + exp(0.5*aux["jj"]^2 - params[:β]* sqrt(O[:Q] ) *t )*Herf(aux["jj"] - aux["b"]*t) ) ) ^2
    IconjQ(t) = Herf(aux["a"]*t)  /  ( 1 +  ( Herf(aux["b"]*t) / Herf(aux["jj"] - aux["b"]*t)) * exp(params[:β]*sqrt(O[:Q]) *t - 0.5*aux["jj"]^2  ) )^2
    return 2 * params[:α] * params[:β]^2 * integrate(X,IconjQ)
end

# inner entropic formula
χI(Model::BinaryPerceptronRS{TSPerceptron}, z) = sqrt( Model.Oconj[:Q] )*z  + Model.Oconj[:R]


# energy
function energy(Model::BinaryPerceptronRS{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    Ienergy(t) = Herf(aux["a"]*t) * (aux["sqrt_1mQ"] * exp(-0.5*aux["b2"]*t^2) / sqrt(2π)  - (params[:β]*(1-O[:Q]) - sqrt(O[:Q])*t) *exp( 0.5*aux["jj"]^2 - params[:β]* sqrt(O[:Q] ) *t) *  Herf(aux["jj"]-aux["b"]*t)) / KV(Model,t)

    return 2*params[:α] * integrate(X,Ienergy)
    #I0 = (O[:Q]*Oconj[:Q] - O[:R] * Oconj[:R]) / params[:β]
    #I1_int(t) = Herf(aux["a"]*t)  /  ( 1 +  ( Herf(aux["b"]*t) / Herf(aux["jj"] - aux["b"]*t)) * exp(params[:β]*sqrt(O[:Q]) *t - 0.5*aux["jj"]^2  ) )
    #I1 = - 2 * params[:α]*params[:β] * integrate(X,I1_int)

    #I2_int(t) = Herf(aux["a"]*t) * exp(-0.5 * aux["b2"] * t^2) / KV(Model, t)
    #I2 = (2*params[:α]/sqrt(2*π)) * (aux["sqrt_1mQ"] + 1/aux["sqrt_1mQ"]) * integrate(X,I2_int)
    #return I0 + I1 + I2
end

function entropy(Model::BinaryPerceptronRS{TSPerceptron}, X::TI) where TI <:IntegrationMethod
    return Model.params[:β]*Model.aux["energy"] - Model.aux["f"]
end

annealed_free_energy(Model::BinaryPerceptronRS{TSPerceptron},X::TI) where TI <: IntegrationMethod = -1.23456




################################## DAT line #############################################
function λE(Model::BinaryPerceptronRS{TSPerceptron},X::TI) where TI <: IntegrationMethod
    @extract Model : params aux O Oconj
    λE_inner(t) = Herf(aux["a"]*t) * ( (params[:β] / (sqrt(2*π) * aux["sqrt_1mQ"])) * ( exp(-0.5*aux["b2"]*t^2) / KV(Model,t)) - params[:β]^2 *  Herf(aux["b"]*t) * exp(0.5*aux["jj"]^2 - params[:β]* sqrt(O[:Q]) *t )*Herf(aux["jj"] - aux["b"]*t) /KV(Model,t)^2  ) ^2
    return  2*params[:α] * integrate(X,λE_inner)
end
