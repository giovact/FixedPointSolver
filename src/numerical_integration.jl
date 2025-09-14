abstract type IntegrationMethod
end

# TODO: add QUADGK see the library in julia, it has a gauss quadrature too

############################################## Dummy ########################################

struct NoInt <: IntegrationMethod
    points::Vector{FT}
end
NoInt() = NoInt(zeros(0))

print_integration_params(X::NoInt) = "_NoInt_"

############################################## HermiteQuadrature (new version) ########################################

struct HermiteQuadrature <: IntegrationMethod
    n::Int
    points::Vector{FT}
    weights::Vector{FT}
    sum_weights::FT
end

function HermiteQuadrature(n::Int)
    x,w = gausshermite(n)
    return HermiteQuadrature(n,sqrt(2) .* x, w,sqrt(π))
end

function integrate(X::HermiteQuadrature,f)
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += f(z)*X.weights[iz]
    end
    return Σ /X.sum_weights
end

function integrate(X::HermiteQuadrature,f,V0::Vector{FT})
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += f(z)*X.weights[iz] / V0[iz]
    end
    return Σ / X.sum_weights
end

integrate(X::HermiteQuadrature, V::Vector{FT}) = dot(V,X.weights) / X.sum_weights

print_integration_params(X::HermiteQuadrature) = string("_Int_",typeof(X), "_npoints_", X.n)



############################################## Legendre Quadrature  ########################################

struct LegendreQuadrature <: IntegrationMethod
    n::Int
    extr::FT
    points::Vector{FT}
    weights::Vector{FT}
    logweights::Vector{FT}
    sumweights::FT
end

function LegendreQuadrature(n::Int;bound::FT = 5.0)
    x,w = gausslegendre(n)
    y = x .* bound
    logweights_including_DZ = log(bound) .+ log.(w)  .- 0.5 * y .^ 2 .- ( 0.5 * log(2π) ) 
    #weights_including_DZ = bound * ( w .* exp.(-0.5*y.^2) ) / sqrt(2*π)
    weights_including_DZ = exp.(logweights_including_DZ)
    return LegendreQuadrature(n,bound, y, weights_including_DZ , logweights_including_DZ, sum(weights_including_DZ))
end

function integrate(X::LegendreQuadrature,f)
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += f(z)*X.weights[iz]
    end
    return Σ / X.sumweights
end

function integrate_log(X::LegendreQuadrature,logf)
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += exp( logf(z) + X.logweights[iz] ) 
    end
    return Σ /X.sumweights
end

function integrate_double(X::LegendreQuadrature,f)
    Σ = 0.0
    for (iz1,z1) in enumerate(X.points)
        for (iz2,z2) in enumerate(X.points)
            Σ += f(z1,z2)*X.weights[iz1]*X.weights[iz2]
        end
    end
    return Σ / X.sumweights
end

function integrate(X::LegendreQuadrature,f,V0::Vector{FT})
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += f(z)*X.weights[iz] / V0[iz]
    end
    return Σ / X.sumweights
end

integrate(X::LegendreQuadrature, V::Vector{FT}) = dot(V,X.weights) / X.sumweights

print_integration_params(X::LegendreQuadrature) = string("_Int_",typeof(X), "_npoints_", X.n, "_extr_", X.extr)

############################################## UniformInterp ########################################
struct UniformInterp <: IntegrationMethod
    n::Int
    k::Int
    extr::FT
    points::Vector{FT}
    weights::Vector{FT}
    sum_weights::FT
    fvalues::Vector{FT}
end

function UniformInterp(n::Int;k::Int = 2)
    extr = sqrt(2*log(k*n))
    x = range(-extr,stop = extr,length = n)
    return UniformInterp(n, k, extr, x, exp.(- x.^2 ./2) , sum(exp.(- x.^2 ./2)), similar(x))
end

function integrate(X::UniformInterp,f)
    Σ = 0.0
    for (iz,z) in enumerate(X.points)
        Σ += f(z)*X.weights[iz]
    end
    return Σ /  X.sum_weights
end
# actually the for loop makes it even quicker

function integrate(X::UniformInterp,f, V0::Vector{FT})
    X.fvalues .= f.(X.points) ./ V0
    return dot(X.weights, X.fvalues ) / X.sum_weights
    #Σ = 0.0
    #for (iz,z) in enumerate(X.points)
    #    Σ += f(z)*X.weights[iz] / V0[iz]
    #end
    #return Σ / X.sum_weights
end

# should be dispatched correctly
integrate(X::UniformInterp, V::Vector{FT}) = dot(V,X.weights) / X.sum_weights

print_integration_params(X::UniformInterp) = string("_Int_",typeof(X), "_npoints_", X.n, "_extr_",X.extr)


######################################################### Trapezoid from Trapz package ##############################################

struct Trapezoid <: IntegrationMethod
    extr::FT
    n::Int
    points::Vector{FT}
    gauss_prefactor::Vector{FT}
    fvalues::Vector{FT}
end

function Trapezoid(extr::FT,n::Int)
    x = range(-extr,stop = extr,length = n)
    return Trapezoid(extr,n, x, exp.(-0.5 .* x.^2) ./ sqrt(2*π), similar(x))
end

function integrate(X::Trapezoid,f)
    X.fvalues .= f.(X.points) .* X.gauss_prefactor
    return trapz(X.points, X.fvalues)
end

function integrate(X::Trapezoid,f,V0::Vector{FT})
    X.fvalues .= f.(X.points) .* X.gauss_prefactor ./ V0
    return trapz(X.points, X.fvalues)
end


###################################################################################################################################


struct BufferedIntegral{TI <: IntegrationMethod} <: IntegrationMethod
    Xinner::TI
    auxV::Dict{String,Vector{FT}}
end

BufferedIntegral(X::TI) where {TI <: IntegrationMethod} = BufferedIntegral(X,Dict{String,Vector{FT}}())

set_integrationmethod(Model::TM, X::TI) where {TM <: FPModel,TI <: IntegrationMethod} = X


######################################## CUBA ######################################

abstract type CubatureWrap <: IntegrationMethod
end

struct HAdapt <: CubatureWrap
    extr::FT

end
HAdapt(b=6.0) = HAdapt(b) 
integrate(X::HAdapt,f) = Cubature.hquadrature(x-> f(x)*exp(-0.5*x^2),-X.extr,X.extr)[1] / sqrt(2π)

#############################

struct PAdapt <: CubatureWrap
    extr::FT
end
PAdapt(b=6.0) = PAdapt(b) 
integrate(X::PAdapt,f) = Cubature.pquadrature(x-> f(x)*exp(-0.5*x^2),-X.extr,X.extr)[1] / sqrt(2π)

#############################


# what you could do is to allocate the arrays in form of a matrix of Xi (z,t) and then integrate. Might be more efficient in terms of allocation
# the only problem is to check wheter this is posisble with a inner structure inside each model. 

############################################## support from NIntegration not working ##############################################
#using NIntegration

#struct NIntfromLibrary <: IntegrationMethod
#    bound::FT
#    reltol::FT
#    abstol::FT
#    maxeval::Int
#end

#NIntfromLibrary(xbound) = NIntfromLibrary(xbound, 1e-6, eps(),Int(1e5))

#function integrate(X::NIntfromLibrary,f)
#    fp(z) = f(z)*exp(- 0.5* z^2 )
#    return nintegrate(fp ,(-X.bound,), (X.bound,); reltol = X.reltol, abstol = X.abstol, maxevals= X.maxeval) / sqrt(2π)
#end


function set_integration_method(m::String, n::Int, b::FT)
    m== "LegendreQuadrature" && return LegendreQuadrature(n; bound = b)
    m== "UniformInterp" && return UniformInterp(n)
    m== "NoInt" && return NoInt()
end





###################################################################################################################################
###################################################################################################################################
struct DebuggerIntegral{TI <: IntegrationMethod} <: IntegrationMethod
    Xinner::TI
    fvalues::Vector{FT}
end
DebuggerIntegral(X::TI) where {TI <: IntegrationMethod} = DebuggerIntegral(X,X.n)

function integrate(X::DebuggerIntegral,f)
    Σ = 0.0
    for (iz,z) in enumerate(X.Xinner.points)
        Σ += f(z)*X.Xinner.weights[iz]
        X.values[iz] = f(z)
    end
    return Σ / X.Xinner.sumweights
end

function integrate_log(X::DebuggerIntegral,logf)
    Σ = 0.0
    for (iz,z) in enumerate(X.Xinner.points)
        Σ += exp( logf(z) + X.Xinner.logweights[iz] ) 
        X.Xinner.values[iz] = exp( logf(z) + X.Xinner.logweights[iz] ) 

    end
    return Σ /X.Xinner.sumweights
end


###################################################################################################################################
###################################################################################################################################
