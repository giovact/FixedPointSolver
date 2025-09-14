struct VariationalMF{TMatrix<:AbstractMatrix} <: FPModel # parallel solution of MF/TAP equations
    N::Int
    params::NamedTuple{(:β,),Tuple{FT}}
    J::TMatrix
    H::Vec
    O::NamedVec
    localH::Vec
    aux::Dict{String,FT}
end

#constructors
function VariationalMF() 
    N = 10
    J = triu(rand(N,N))
    J += J'
    J[diagind(J)] .= 0.0
    H = rand(N) .- 0.5
    return VariationalMF(0.5, J, H, sign.(H))
end

function VariationalMF(β,J,H,O) 
    !issymmetric(J) && throw("J not symmetric")
    @assert size(J,1) == length(H)
    @assert length(O) == length(O)
    return VariationalMF(length(H),(;β=β), J,H,NamedArray(O,collect([Symbol(:M,i) for i=1:length(H)])),zeros(length(H)),Dict{String,FT}())
end

create_model(Model::VariationalMF, params) = VariationalMF(N,params, Model.J, Model.H, Model.O,Model.localH,Dict{String, FT}())

########################################################## SELF-CONSISTENT EQUATIONS #########################################################

function lhs!(Model::VariationalMF,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : J H params O localH
    mul!(localH,J,O,params[:β],0.0)
    localH .+= params[:β] * H
    Onew .= tanh.(localH) 
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::VariationalMF,X::TI) where TI <: IntegrationMethod
    return -1.0 
end
