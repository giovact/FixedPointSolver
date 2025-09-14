struct RepCW <: FPModel
    params::NamedTuple{(:β,:J,:h,:γ,:y,:gy,:βgamma),Tuple{FT,FT,FT,FT,Int,Int,Bool}}
    O::NamedVec
    InnerCW::Spin1Half
    aux::Dict{String,FT} 
end

# constructors
RepCW() = RepCW(1.0,0.1,0.1,0.1,2,2,false,ones(2)) # dummy constructor for testing
RepCW(β,J,h,γ,y,gy,βgamma) = RepCW(β,J,h,γ,y,gy,βgamma,zeros(y))
RepCW(β,J,h,γ,y,gy,βgamma,O) = RepCW((β=β,J=J,h=h,γ=γ,y=y,gy=gy,βgamma=βgamma),NamedArray(O,vcat([Symbol(:M,u) for u=1:y]) ), Spin1Half(y), Dict{String, FT}())


create_model(Model::RepCW, params) = RepCW(params, Model.O, Spin1Half(Model.params[:y]), Model.aux)
set_integrationmethod(Model::RepCW,X::TI) where TI <: IntegrationMethod = X

########################################################## SELF-CONSISTENT EQUATIONS #########################################################


function effective_parameters!(Model::RepCW)
    @extract Model : params O aux

    fill!(Model.InnerCW.H,0.0)
    fill!(Model.InnerCW.J,0.0)
    for u=1:params[:y]
        Model.InnerCW.H[u] = params[:β]*params[:J]*O[Symbol(:M,u)] + params[:β]*params[:h]
        for v=1:params[:y]
            if u!=v
                Model.InnerCW.J[u,v] =  params[:βgamma] ? (params[:β] * params[:γ] / params[:gy]) : params[:γ] / params[:gy]
            end
        end
    end
end


function lhs!(Model::RepCW,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    effective_parameters!(Model)
    traceSpin1Half(Model.InnerCW)
    for u=1:Model.params[:y]
        Onew[Symbol(:M,u)] = Model.InnerCW.M[u]
    end
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::RepCW,X::TI) where TI <: IntegrationMethod
    @extract Model : params O

    effective_parameters!(Model)
    Z = traceSpin1Half(Model.InnerCW)

    return 0.5 * (params[:J] / params[:y]) * dot(O,O) - log(Z) / (params[:β] * params[:y])

end

_isCWmetastable(Model::RepCW, X::TI, Oold::NamedVec)  where TI <: IntegrationMethod = sign(Model.params[:h] * Model.O[:M]) < 0