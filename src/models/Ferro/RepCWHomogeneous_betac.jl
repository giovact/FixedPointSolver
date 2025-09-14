struct RepCWHomogeneous_betac <: FPModel
    params::NamedTuple{(:J,:γ,:y,:gy,:βgamma),Tuple{FT,FT,Int,Int,Bool}}
    O::NamedVec
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} 
end

# constructors
RepCWHomogeneous_betac() = RepCWHomogeneous_betac(1.0,1.0,2,2,true,ones(1)) # dummy constructor for testing
RepCWHomogeneous_betac(J,γ,y,gy,βgamma) = RepCWHomogeneous_betac(J,γ,y,gy,βgamma, ones(1))
RepCWHomogeneous_betac(J,γ,y,gy,βgamma,O) = RepCWHomogeneous_betac((J=J,γ=γ, y=y,gy=gy,βgamma=βgamma),NamedArray(O,[:βc]), HomogeneousSpin1Half(y), Dict{String, FT}())

create_model(Model::RepCWHomogeneous_betac, params) = RepCWHomogeneous_betac(params, Model.O, HomogeneousSpin1Half(Model.params[:y]), Model.aux)

########################################################## SELF-CONSISTENT EQUATIONS #########################################################


function effective_parameters!(Model::RepCWHomogeneous_betac)
    @extract Model : params O aux

    aux["heff"] = 0.0
    aux["Jeff"] = params[:βgamma] ? (O[:βc] * params[:γ] / params[:gy]) : params[:γ] / params[:gy]
end

function lhs!(Model::RepCWHomogeneous_betac,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    @extract Model : params O aux
    effective_parameters!(Model)
    Onew[:βc] = 1 / ( 1 + (params[:y]-1)*trace_CurieWeiss(Model.InnerCW,Model.aux["heff"],Model.aux["Jeff"])[3] ) 
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::RepCWHomogeneous_betac,X::TI) where TI <: IntegrationMethod
    return -1.
end
