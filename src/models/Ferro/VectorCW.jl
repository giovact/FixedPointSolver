struct VectorCW <: FPModel
    params::NamedTuple{(:J,),Tuple{FT}}
    H::Vec
    O::NamedVec
    minusbetaF::Vec
    Msquared::Vec
    aux::Dict{String,FT}
end

#constructors
VectorCW() = VectorCW(0.5,zeros(2),ones(2)) # dummy constructor for testing
VectorCW(J,H,O) = VectorCW((;J=J), H,NamedArray(O,collect([Symbol(:M,i) for i=1:length(H)])),zeros(length(H)),zeros(length(H)),Dict("J"=>0.0))

create_model(Model::VectorCW, params) = VectorCW(params, Model.H, Model.O, zeros(length(Model.H)),zeros(length(H)),Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function lhs!(Model::VectorCW,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
   Onew .= (Model.params[:J]+Model.aux["J"]).*Model.O .+ Model.H .|> tanh 
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::VectorCW,X::TI) where TI <: IntegrationMethod
    @extract Model : params H minusbetaF O Msquared aux
    Msquared .= O.^2
    minusbetaF .= (params[:J]+aux["J"])*Model.O .+ Model.H .|> log2cosh 
    minusbetaF .-= 0.5 * (params[:J]+aux["J"]) .* Msquared
    return -1.0 
end
