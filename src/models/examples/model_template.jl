
"""
Basic Template Model 

struct TemplateModel <: FPModel
    # required fields 
    params::NamedTuple{(:a,:b,:c,:d),Tuple{FT,FT,Int,Bool}} # control parameters (e.g. temperature, global external fields) -> 
    O::NamedVec #vector of orderparameters 
    aux::Dict{String,FT} # other stuff, among them free energy
end

# constructors
TemplateModel() = TemplateModel(42.0,0.0,42,true,ones(10)) # dummy constructor for testing

function TemplateModel(a,b,c,d,O)  # constructor with parameters a,b,c,d and input initial condition O (which must be an array)
    TemplateModel((a=a,b=b,c=c,d=d),
                    NamedArray(O,[:O1;O2;...;:On]),  # declare NamedArray from the initial condition O (given as an input) and the names of the order parameters
                    NamedArray(O,[:O1;O2;...;:On]),  # declare NamedArray from the initial condition O (given as an input) and the names of the order parameters
                    Dict{String, FT}()) # empty dictionary -> aux
end

create_model(Model::TemplateModel, params) = TemplateModel(params, Model.O,Model.Oconj,Dict{String, FT}())
# this function is needed to scan the model over a parameter (among params)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function expr_Op_1(Model::TemplateModel,X::TI) where TI <: IntegrationMethod
    ....
end

......

function expr_Op_n(Model::TemplateModel,X::TI) where TI <: IntegrationMethod
    ....
end


function lhs!(Model::TemplateModel,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : O
    Onew will be of the same kind of O, i.e. a NamedArray with the same names
    Onew[:O1] = expr_Op_1(Model,X)
    Onew[:O2] = expr_Op_2(Model,X)
    ..........................
    Onew[:On] = expr_Op_n(Model,X)
end

########################################################## FREE ENERGY #########################################################
function free_energy(Model::TemplateModel,X::TI) where TI <: IntegrationMethod
    # must be a scalar of type FT
    
    return 12345.0
end

########################################################## Auxiliary Observables #########################################################

function compute_aux!(Model::TemplateModel,X::TI)  where TI <: IntegrationMethod
    Model.aux["additional_observable"] = ...
end
"""
