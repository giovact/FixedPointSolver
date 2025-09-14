
function plotres(X::ScannedModel)
    @extract X : Model scan data
    @extract scan : param pvalues protocol reverseordering
    dims = has_conj_order_params(Model) ? (2,2) : (1,3)
    sizes =  has_conj_order_params(Model) ? (10,10) : (6.6*3,4.8)
    fig, axs = subplots(dims[1], dims[2], figsize = sizes )
    fig.subplots_adjust(wspace = 0.3)
    for (i,x) in enumerate(names(Model.O)[1])
        axs[1].plot(pvalues, data[:,i], label = string(x))
    end
    ax2r = axs[3].twinx()
    axs[3].plot(pvalues, data[:,n_order_params(Model)+3], color = "blue")
    ax2r.plot(pvalues, data[:,n_order_params(Model)+2], linestyle = "dashed", color = "red", lw = 0.8)
    
    axs[2].plot(pvalues, data[:,n_order_params(Model)+1], color = "tab:green", label ="f")
   
    axs[1].legend()
    axs[2].set_ylabel("free energy")
    axs[3].set_ylabel("iter")
    ax2r.set_ylabel("conv")
    ax2r.tick_params(axis="y", colors="tab:red")
    ax2r.spines["right"].set_color("tab:red")

    if has_conj_order_params(Model)
        for (i,x) in enumerate(names(Model.Oconj)[1])
            axs[4].plot(pvalues, data[:,n_order_params(Model)+3+i], label = string(x,"_conj"))
        end
        axs[4].set_title("Conjugate Order Params")
        axs[4].legend()

    end

    fig.suptitle(string(typeof(Model), "  ", params_to_textstring(Model;separator = " - ", exclude = [param], equalsymbol = "=")[1:end-3]) ) 
    for iax=1:length(axs)
        axs[iax].set_xlabel(param)
    end
    axs[1].set_title("Order Params")
    axs[2].set_title("Free energy")
    
    
end




function plotres(data)
    dims = sum(n_conj_order_params(single.Model) for single in data) > 0  ? (2,2) : (1,3)
    sizes =  sum(n_conj_order_params(single.Model) for single in data) > 0  ? (15,8) : (6.6*3,4.8)
    
    
    fig, axs = subplots(dims[1], dims[2], figsize = sizes )
    ax2r = axs[3].twinx()

    fig.subplots_adjust(wspace = 0.3)
    
    str_title = ""
    for (ii, singleout) in enumerate(data)
        for (io,o) in enumerate(names(singleout.Model.O)[1])
            axs[1].plot(singleout.scan.pvalues, singleout.data[:,io], label = string(o, "- ",ii))
        end
        axs[2].plot(singleout.scan.pvalues, singleout.data[:,n_order_params(singleout.Model)+1], label =string(ii))
        axs[3].plot(singleout.scan.pvalues, singleout.data[:,n_order_params(singleout.Model)+3], label = ii)
        ax2r.plot(singleout.scan.pvalues, singleout.data[:,n_order_params(singleout.Model)+2], lw = 0.8, label = ii)

        if has_conj_order_params(singleout.Model)
            for (io,o) in enumerate(names(singleout.Model.Oconj)[1])
                axs[4].plot(singleout.scan.pvalues, singleout.data[:,n_order_params(singleout.Model)+3+io], label = string(o,"_conj", "-",ii))
            end
            axs[4].set_title("Conjugate Order Params")
            axs[4].legend()
        end
        str_title = string(str_title, " \n ",ii,"-->",string(typeof(singleout.Model), "  ", params_to_textstring(singleout.Model;separator = " - ", exclude = [singleout.scan.param], equalsymbol = "=")[1:end-3], " protocol: ", singleout.scan.protocol, " - init -> ", print_op(singleout.Model)) )

    end
    
    axs[2].set_ylabel("free energy")
    axs[3].set_ylabel("iter")
    ax2r.set_ylabel("conv")
    ax2r.tick_params(axis="y", colors="tab:red")
    ax2r.spines["right"].set_color("tab:red")
 
    fig.suptitle(str_title )
    params = collect([singleout.scan.param for singleout in data])
    @assert length(unique(params)) == 1 
    for iax in eachindex(axs)
        axs[iax].legend()
        axs[iax].set_xlabel(params[1])
    end
    axs[1].set_title("Order Params")    
end

