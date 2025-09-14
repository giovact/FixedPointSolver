struct RepCWHomogeneous <: FPModel
    params::NamedTuple{(:β,:J,:h,:γ,:y,:gy,:βgamma),Tuple{FT,FT,FT,FT,Int,Int,Bool}}
    O::NamedVec
    InnerCW::HomogeneousSpin1Half
    aux::Dict{String,FT} 
end

# constructors
RepCWHomogeneous() = RepCWHomogeneous(1.0,0.1,0.1,0.1,2,2,false,ones(1)) # dummy constructor for testing
RepCWHomogeneous(β,J,h,γ,y,gy,βgamma) = RepCWHomogeneous(β,J,h,γ,y,gy,βgamma,zeros(1))
RepCWHomogeneous(β,J,h,γ,y,gy,βgamma,O) = RepCWHomogeneous((β=β,J=J,h=h,γ=γ,y=y,gy=gy,βgamma=βgamma),NamedArray(O,[:M]), HomogeneousSpin1Half(y), Dict{String, FT}())

create_model(Model::RepCWHomogeneous, params) = RepCWHomogeneous(params, Model.O, HomogeneousSpin1Half(Model.params[:y]), Model.aux)
set_integrationmethod(Model::RepCWHomogeneous,X::TI) where TI <: IntegrationMethod = X

########################################################## SELF-CONSISTENT EQUATIONS #########################################################


function effective_parameters!(Model::RepCWHomogeneous)
    @extract Model : params O aux

    aux["heff"] = params[:β]*params[:J]*O[:M] + params[:β]*params[:h]
    aux["Jeff"] = params[:βgamma] ? (params[:β] * params[:γ] / params[:gy]) : params[:γ] / params[:gy]
end


function lhs!(Model::RepCWHomogeneous,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    effective_parameters!(Model)
    Onew[:M] = trace_CurieWeiss(Model.InnerCW,Model.aux["heff"],Model.aux["Jeff"])[2]
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::RepCWHomogeneous,X::TI) where TI <: IntegrationMethod
    @extract Model : params O

    effective_parameters!(Model)
    return 0.5 * params[:J] * O[:M]^2 - log(trace_CurieWeiss(Model.InnerCW,Model.aux["heff"],Model.aux["Jeff"])[1] ) / (params[:β] * params[:y])

end

_isCWmetastable(Model::RepCWHomogeneous, X::TI, Oold::NamedVec)  where TI <: IntegrationMethod = sign(Model.params[:h] * Model.O[:M]) < 0