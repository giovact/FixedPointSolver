struct HopfieldExample <: FPModel
    params  ::NamedTuple{(:β,:α),Tuple{FT,FT}}  # *REQUIRED*  # NamedTuple containing the control parameters of the model (here β, α) 
    O       ::NamedVec                                    # *REQUIRED*  # NamedArray containing the order parameters (here M,Q) 
    Oconj   ::NamedVec                                                  # NamedArray containing the order parameters (here M,Q)
    aux     ::Dict{String,FT}                        # *REQUIRED*  # dictionary used to store additional observables (among them, the free energy) 
end

# constructors
HopfieldExample() = HopfieldExample(.1, .1,ones(2)) # dummy constructor for testing
HopfieldExample(β,α) = HopfieldExample(β,α,zeros(2))
HopfieldExample(β,α,O) = HopfieldExample((β=β,α=α),NamedArray(O,[:M;:Q]),NamedArray(zeros(1),[:Q]), Dict{String,FT}())

create_model(Model::HopfieldExample, params) = HopfieldExample(params,Model.O, Model.Oconj, Dict{String,FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
Delta(Model::HopfieldExample) =  1 - Model.params[:β] * (1 - Model.O[:Q])
Deltas!(Model::HopfieldExample) = Model.aux["Δ"] = clamp(Delta(Model), 1e-30, Inf)
is_free_energy_definite(Model::HopfieldExample) = Delta(Model)>0

#### Conjugate order params SP
function exprOconj!(Model::HopfieldExample)
    @extract Model : O Oconj
    Deltas!(Model)
    #Qhat
    Oconj[:Q] = O[:Q] / Model.aux["Δ"]^2
end

# order params SP
function exprM(Model::HopfieldExample,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    O[:M] == 0.0 && return 0.0
    IM(z) = tanh(  params[:β]*(O[:M]  + sqrt(params[:α]*Oconj[:Q]) * z) )
    integrate(X,IM)
end

function exprQ(Model::HopfieldExample,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj
    IQ(z) = tanh(  params[:β]*(O[:M]  + sqrt(params[:α]*Oconj[:Q]) * z) )^2
    integrate(X,IQ)
end


function lhs!(Model::HopfieldExample,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprOconj!(Model)

    Onew[:M] = exprM(Model,X)
    Onew[:Q] = exprQ(Model,X)
end


function free_energy(Model::HopfieldExample,X::TI) where TI <: IntegrationMethod
    @extract Model : params O Oconj aux
    exprOconj!(Model)
    
    energy_signal= 0.5 * O[:M]^2
    fourier_overlaps =  0.5*params[:α]*params[:β]*(1-O[:Q])*Oconj[:Q]

    logdetK = (params[:α]/(2*params[:β])) * log(aux["Δ"])  
    other =  0.5*params[:α]* O[:Q] / aux["Δ"] 

    IS(z) = log2cosh(  params[:β]*(O[:M]  + sqrt(params[:α]*Oconj[:Q]) * z) )
    entropy = integrate(X,IS)

    return  energy_signal + fourier_overlaps + logdetK - other - entropy/params[:β]
end

function compute_aux!(Model::HopfieldExample,X::TI) where TI <: IntegrationMethod
    Model.aux["additional_observables"] = 42.0        
end


_isRetrieval(Model::HopfieldExample, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod = Model.O[:M] > 1e-2

function _isRetrievalStateDominant(Model::HopfieldExample, X::TI, Oold::NamedVec,tol::FT) where TI <: IntegrationMethod
    ModelSG = create_model(deepcopy(Model), Model.params)
    ModelSG.O[:Q] = 1.0 
    ModelSG.O[:M] = 0.0
    solve_SPequations!(ModelSG, X; tol = tol, printprogress=false, showinfo=false, K = FixedPoint(0.95), niter = 10000,force_dict = Dict(:M=>0.0))
    return Model.aux["f"] < ModelSG.aux["f"]
end
