# dummy constructor for testing
BinaryPerceptronRS{TSGibbs}() = BinaryPerceptronRS(1.0,1.0,TSGibbs(),ones(2)*1e-2 )

function aux_potential!(Model::BinaryPerceptronRS{TSGibbs})
    if Model.V.storage==true
        Model.aux["a"] = 0.0
    end
    Model.aux["1minus_exp_minusbeta"] = 1-exp(-Model.params[:β])
end

# KV this is what you called K in the notes, at least for the RS
function KV(Model::BinaryPerceptronRS{TSGibbs},t)
    @extract Model : params aux O
    return exp(-params[:β]) + (1-exp(-params[:β]))*Herf(aux["b"]*t)
end

# Equation for Rhat
function exprRhat(Model::BinaryPerceptronRS{TSGibbs},X::TI)  where TI <: IntegrationMethod 
    @extract Model : params aux  O
    Model.V.storage == true && return 0.0

    IconjR(z) = exp(- 0.5 * aux["v"]^2 * z^2 )  / (exp(-params[:β]) + (1-exp(-params[:β]))*Herf(aux["v"]*z))
    return (params[:α] / π ) * (aux["1minus_exp_minusbeta"] / aux["sqrt_1mQ"] ) * integrate(X,IconjR)
end

# Equation for Qhat
function exprQhat(Model::BinaryPerceptronRS{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : params aux O
    IconjQ(t) = Herf( aux["a"]*t) * exp(- aux["b2"] * t^2 )  / ( KV(Model,t) )^2
    return (params[:α] / π ) * (aux["1minus_exp_minusbeta"] / (1-O[:Q]) ) * integrate(X,IconjQ)
end

# inner entropic formula
χI(Model::BinaryPerceptronRS{TSGibbs}, z) = sqrt( Model.Oconj[:Q] )*z  + Model.Oconj[:R]*(1 - Model.V.storage)

# equation for R 
function exprR(Model::BinaryPerceptronRS{TSGibbs},X::TI) where TI <: IntegrationMethod
    Model.V.storage == true && return 0.0
    IR(z) = tanh(χI(Model,z))
    return integrate(X,IR)
end

# energy
function energy(Model::BinaryPerceptronRS{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : params aux 
    IU(t) = Herf(aux["a"]*t) *Herf(-aux["b"]*t) / KV(Model,t)
    return  2*params[:α]*exp(-params[:β]) * integrate(X,IU)
end

function entropy(Model::BinaryPerceptronRS{TSGibbs}, X::TI) where TI <:IntegrationMethod
    Model.params[:β] == Inf && return - Model.aux["f"]
    return Model.params[:β]*Model.aux["energy"] - Model.aux["f"]
end

# annealed free energy
function annealed_free_energy(Model::BinaryPerceptronRS{TSGibbs},X::TI) where TI <: IntegrationMethod
    @extract Model : params
    return - (log(2)+ params[:α]* log(KV(Model,0.0))) / params[:β]
end

################################## DAT line #############################################
function λE(Model::BinaryPerceptronRS{TSGibbs},X::TI) where TI <: IntegrationMethod
   return 0.0
end

