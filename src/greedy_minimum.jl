tuplejoin(t1::Tuple, t2::Tuple, t3...) = tuplejoin((t1..., t2...), t3...)
tuplejoin(t::Tuple) = t

function greedy_search_free_energy(Model0::TM,X::TI; 
                                    thetavec::Vector{FT} = collect(0.1:0.2:1.0),
                                    ranges::Dict = Dict(O=>collect(range(0.0, stop = 1.0, length = 3)) for O in names(Model0.O,1) ),
                                    showinfo::Bool = true,
                                    todiscard::Function = (x-> isnan(x) || isinf(x))) where {TM <: FPModel, TI<:IntegrationMethod}
   
    nop = n_order_params(Model0)
    @assert nop == length(keys(ranges))

    free_energy_values = Matrix{FT}(undef, 0, 1 + n_order_params(Model0) + 1) # first is theta, then nop Order params, then free energy

    violates_bounds = 0
    discardedbyfunction = 0
    
    meshlengths = Tuple(length(thetavec))
    for O in names(Model0.O,1)
        meshlengths = tuplejoin(meshlengths, Tuple(length(ranges[O])))
    end

    @showprogress for CI in CartesianIndices( meshlengths )
        thetav = thetavec[CI.I[1]]
        Model = create_model(Model0, reset_parameter(Model0.params, :θ,thetav ) )
        
        for (iO,Oname) in enumerate(names(Model.O,1))
            Model.O[Oname] = ranges[Oname][CI.I[1+iO]]
        end
        if is_free_energy_definite(Model)
            fv = free_energy(Model, X)
            if !todiscard(fv)
                free_energy_values = cat(free_energy_values, [thetav Vector(Model.O)' fv], dims = 1)
            else
                discardedbyfunction+=1
            end
        else
            violates_bounds +=1
        end
    end

    ntot = prod(meshlengths)

    percentage_remaining =  (1 - (violates_bounds + discardedbyfunction)/ntot)*100
    showinfo && @info " Number of original mesh points is $ntot -> discarded $violates_bounds because do not satisfy bounds"
    showinfo && @info "Among theese, $discardedbyfunction do not satisfy todiscard function, in total the percentage of remaining is $percentage_remaining %"
    return free_energy_values
end



function greedy_oldsearch_free_energy(Model0::TM,X::TI;
                                    fix_op = Dict(),
                                    npointseachO::Int = 10, 
                                    thetavec::Vector{FT} = collect(0.1:0.2:1.0),
                                    ranges::Dict = Dict(O=>collect(range(0.0, stop = 1.0, length = npointseachO)) for O in names(Model0.O,1) ),
                                    todiscard::Function = (x-> isnan(x) || isinf(x))) where {TM <: FPModel, TI<:IntegrationMethod}
   
    list_OP = filter( x -> x ∉ collect(keys(fix_op)), names(Model0.O,1) )
    Oprange = collect(range(0.0, stop = 1.0, length = npointseachO))

    # i can enfore fix_op[:M] =
    free_energy_values = Matrix{FT}(undef, 0, 1+n_order_params(Model0)+1)
    namesOPtofix = collect(keys(fix_op))
    violates_bounds = 0
    discardedbyfunction = 0
    @showprogress for CI in CartesianIndices( tuplejoin(Tuple(length(thetavec)),Tuple(fill(npointseachO,length(list_OP) )) ))
        thetav = thetavec[CI.I[1]]
        Model = create_model(Model0, reset_parameter(Model0.params, :θ,thetav ) )
        itofix = 1
        for Oname in names(Model.O,1)
            if Oname ∈ namesOPtofix
                Model.O[Oname] = fix_op[Oname]
            else
                Model.O[Oname] = Oprange[CI.I[1+itofix]]
                itofix+=1
            end
        end
        if is_free_energy_definite(Model)
            fv = free_energy(Model, X)
            if !todiscard(fv)
                free_energy_values = cat(free_energy_values, [thetav Vector(Model.O)' fv], dims = 1)
            else
                discardedbyfunction+=1
            end
        else
            violates_bounds +=1
        end
    end

    ntot = length(thetavec)*npointseachO^length(list_OP)

    percentage_remaining =  (1 - (violates_bounds + discardedbyfunction)/ntot)*100
    @info " Number of original mesh points is $ntot -> discarded $violates_bounds because do not satisfy bounds"
    @info "Among theese, $discardedbyfunction do not satisfy todiscard function, in total the percentage of remaining is $percentage_remaining %"
    return free_energy_values
end
