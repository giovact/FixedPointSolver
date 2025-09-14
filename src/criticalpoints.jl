struct ScanCritical
    param::Symbol
    p0::FT
    Δp::FT
    Δptol::FT
    p_damp::FT
    bounds::Vector{FT}
    protocol::Symbol
    increase::Bool
    ScanCritical(param, p0, Δp, Δptol, p_damp,bounds,protocol, increase) = protocol ∉ [:follow; :restart] ? throw("select valid protocol") : new(param, p0, Δp, Δptol, p_damp,bounds,protocol, increase)
end

struct ScannedCriticalModel{TM<:FPModel}
    Model::TM
    scan::ScanCritical
    data::Matrix{FT}
    names::Vector{String}
    idxfinal::Int
    idxfinallastconv::Int
    pfinal::FT
    pfinallastconv::FT
    free_energy_final::FT
    free_energy_finallastconv::FT
end


function scanning_info(scan::ScanCritical, Modelstart)
    scan.protocol == :restart && return string(" RESTART - each time from -> ", print_op(Modelstart) )
    scan.protocol == :follow && return string(" FOLLOW fps;")
end


function parameter_change_proposal(p_current,Δp,increase,pbounds,p_damp)
    ptrial = p_current + ( (-1) ^ (1 - increase)) * Δp  # propose a change in parameter p

    if ptrial ∉ pbounds 
        return parameter_change_proposal(p_current,Δp * p_damp,increase,pbounds,p_damp)
    else
        return ptrial, Δp 
    end
end

function findcriticalpoint(Modelstart::TM,scan::ScanCritical,keepgoing::Function,X::TI; # the function keepgoing depends on the Model, and X
                                niter::Int = 10000,
                                K::TU = FixedPoint(0.9),
                                tol::FT = 1e-5,
                                Kfirstrun::TU=FixedPoint(0.95),
                                eps_reset::FT = -Inf,
                                printprogress::Bool = true,
                                fsave::String = "",
                                tosave::Vector{String} = String[],
                                force_dict = Dict{Symbol, FT}(),
                                ) where {TM<:FPModel, TI<:IntegrationMethod, TU<:Solver }
    
    @extract scan : param p0 Δp Δptol p_damp bounds protocol increase
    
    Model = create_model(deepcopy(Modelstart), reset_parameter(Modelstart.params, param,p0) )
    @assert Model.params[param] == p0 #debug

    O0 = deepcopy(Modelstart.O)
    
    infotoprint = string("Scanning over parameter ",param, "  ", print_params(Model; separator = " - ", exclude =[param]), " to find critical point of function ",keepgoing,"  ",scanning_info(scan, Modelstart))
    prog = printprogress ? ProgressUnknown( infotoprint ) : nothing
    
    # do the first run 
    teps = @elapsed conv, iter, vareps = solve_SPequations!(Model,X; K = scan.protocol == :restart ? K : Kfirstrun, tol = tol, niter = niter,showinfo=false,printprogress=false,force_dict = force_dict)
    obs_tosave = ifelse(isempty(tosave),zeros(0),copy([Model.aux[oo] for oo in tosave]))
    #initialize data matrix 
    

    p = p0
    plastconv = p0
    pbounds = Interval(bounds[1], bounds[2])
    if p0 ∉ pbounds 
        throw("Error, initial parameter not in specified bounds. We have p0=$p0 and pbounds=$bounds")
    end
    Otemp = deepcopy(Model.O)
    data = FT.(collect(copy([p0 Vector(Model.O)' Model.aux["f"] conv iter keepgoing(Model,X,Otemp,tol) get_Oconj(Model)' Vector(obs_tosave)' teps])))
    
    names_keys = [string(scan.param); string.(names(Model.O,1)); "f"; "conv";"iter"; "keepgoing";]
    if has_conj_order_params(Model)
        names_keys = vcat(names_keys, get_Oconjnames(Model))
    end
    names_keys = vcat(names_keys,[tosave; "elapsed"])

    idx_iter, idxfinal, idxfinal_lastconv = 1, 1,1
    while Δp > Δptol
        ptrial,Δp = parameter_change_proposal(p,Δp,increase,pbounds,p_damp)
        Otemp .= deepcopy(Model.O)
        if protocol == :restart
            Model.O .= deepcopy(O0)
        elseif protocol == :follow
            reset_orderparams!(Model.O,eps_reset)
        end

        Model = create_model(Model, reset_parameter(Model.params, param,ptrial ) ) # maybe this line can be put before the if, it should change nothing but ok

        teps = @elapsed conv, iter, vareps = solve_SPequations!(Model,X; K = K, tol = tol, niter = niter,showinfo=false,printprogress=false,force_dict = force_dict)
        obs_tosave = ifelse(isempty(tosave),zeros(0),copy([Model.aux[oo] for oo in tosave]))
        @assert Model.params[param] == ptrial
        
        goon = keepgoing(Model,X,Otemp,tol)
        data = cat(data, collect(copy([ptrial Vector(Model.O)' Model.aux["f"] conv iter goon get_Oconj(Model)' Vector(obs_tosave)' teps])) , dims = 1)
        
        # il problema è che ti devi salvare pure i prametri d'ordine di prima, altrimenti lui riparte da quelli se è follow -> deve essere follow ma con storia
        # condition here 
        # saving if fsave >0
        (length(fsave)>0) && writedlm(fsave,data)
        if goon # such that you can go further
            p = ptrial
            idxfinal = idx_iter + 1
            if conv == 1
                idxfinal_lastconv = idx_iter+1
                plastconv = ptrial
            end
        else
            Δp *= p_damp
            if protocol == :follow
                Model.O .= deepcopy(Otemp)
            end
        end
        printprogress && ProgressMeter.next!(prog; showvalues = [print_op_withparam(Model,param); (:f,Model.aux["f"]); (:iter,iter); (:conv, conv); (:Δp,Δp)] )
        idx_iter+=1
    end

    return ScannedCriticalModel(Model, scan, data,names_keys,idxfinal,idxfinal_lastconv,p,plastconv,data[idxfinal,1+length(Model.O)+1],data[idxfinal_lastconv,1+length(Model.O)+1])
end


##### PLOTTING

function showtable(S::ScannedCriticalModel)
    dataC = S.data
    hl_1 = Highlighter(f = (dataC,i,j) -> dataC[i,3+n_order_params(S.Model)]==1, crayon = crayon"green")
    hl_0 = Highlighter(f = (dataC,i,j) -> dataC[i,3+n_order_params(S.Model)]==0, crayon = crayon"yellow")
    hl_m1 = Highlighter(f = (dataC,i,j) -> dataC[i,3+n_order_params(S.Model)]==-1, crayon = crayon"red")

    hl_final = Highlighter(f = (dataC,i,j) -> i==S.idxfinal, crayon = Crayon(background = :light_yellow, bold = true))
    hl_finallastconv = Highlighter(f = (dataC,i,j) -> i==S.idxfinal, crayon = Crayon(background = :light_green, bold = true))

    pretty_table(S.data, header = S.names,  header_crayon = crayon"blue bold", highlighters = (hl_finallastconv,hl_final,hl_m1,hl_0,hl_1))
end

