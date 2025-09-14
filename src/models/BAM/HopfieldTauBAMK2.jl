struct HopfieldTauBAMK2 <: FPModel # for the moment patterns are decorrelated and binary
    params::NamedTuple{(:β,:γ,:τ),Tuple{FT,FT, FT}}
    O::NamedVec #vector of orderparameters -> M #
    aux::Dict{String,FT} # other stuff, among them free energy
    patterns::Matrix{Int} # assuming binary
end

#constructors
HopfieldTauBAMK2() = HopfieldTauBAMK2(0.1,1.0,0.1,ones(4)) # dummy constructor for testing
HopfieldTauBAMK2(β,γ,τ) = HopfieldTauBAMK2(β,γ,τ,zeros(4))
HopfieldTauBAMK2(β,γ,τ,O) = HopfieldTauBAMK2((β=β,γ=γ,τ=τ), NamedArray(O,[:M1;:M2;:Mb1;:Mb2]),Dict{String,FT}(),hcat(collect([tovec(x,2) for x= 0 : 3])...) )

create_model(Model::HopfieldTauBAMK2, params) = HopfieldTauBAMK2(params, Model.O,Dict{String,FT}(),patterns)

function set_initial_condition!(Model::HopfieldTauBAMK2; Mmin::FT=0.9)
    @extract Model : O
    for (µ,_) in enumerate(O)
        O[µ] = rand()*(1-Mmin) + Mmin
    end 
end

########################################################## SELF-CONSISTENT EQUATIONS #########################################################
function Ξ(Model::HopfieldTauBAMK2, ξ)
    @extract Model : O params
    return Model.params[:β] * dot( (1-params[:τ])*O[[:M1;:M2]] .+ (params[:τ]/params[:γ])*O[[:Mb1;:Mb2]], ξ)
end

function Ξb(Model::HopfieldTauBAMK2, ξb)
    @extract Model : O params
    return Model.params[:β] * dot( (1-params[:τ])*O[[:Mb1;:Mb2]] .+ (params[:τ]*params[:γ])*O[[:M1;:M2]], ξb)
end



function exprMall!(Model::HopfieldTauBAMK2,Onew::NamedVec ) 
    @extract Model : patterns O
    fill!(Onew, 0.0)
    for ξ in eachcol(patterns)
        Onew[1:2] .+= ξ .* tanh(Ξ(Model,ξ))
        Onew[3:4] .+= ξ .* tanh(Ξb(Model,ξ))
    end
    Onew ./= 4 # assuming to be independent
end

function lhs!(Model::HopfieldTauBAMK2,X::TI,Onew::NamedVec) where TI <: IntegrationMethod
    exprMall!(Model,Onew)
end

########################################################## FREE ENERGY #########################################################

function free_energy(Model::HopfieldTauBAMK2,X::TI) where TI <: IntegrationMethod
    @extract Model : params patterns O
    energy1 = (0.5*params[:γ])*(1-params[:τ])*(O[:M1]^2 + O[:M2]^2)
    energy2 = (0.5/params[:γ])*(1-params[:τ])*(O[:Mb1]^2 + O[:Mb2]^2)
    energymix = params[:τ]*(O[:M1]*O[:Mb1] + O[:M2]*O[:Mb2])

    entropy1 = sum( logcosh(Ξ(Model, ξ)) for ξ in eachcol(patterns)) / 4
    entropy2 = sum( logcosh(Ξb(Model, ξ)) for ξ in eachcol(patterns)) / 4
    return energy1 + energy2 + energymix - (params[:γ]/params[:β])*entropy1 - (1/(params[:γ]*params[:β]))*entropy2 
end
