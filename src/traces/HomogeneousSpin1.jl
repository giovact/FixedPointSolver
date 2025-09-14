

############################################################

struct HomogeneousSpin1
    N::Int
    p0::SingleVar
    MandQ::Vector{Tuple{Int,Int}} 
    P0::Dict{Tuple{Int,Int},FT}
    dos::Dict{Tuple{Int,Int},BigInt}
    logdos::Dict{Tuple{Int,Int},FT} # occhio che qui sta cosa da errori grandi per N grande (tra exp(logdos) e dos). Quindi boh
end

function HomogeneousSpin1(N::Int,p0::SingleVar)
    @assert N>1
    @assert p0.ξ0 == [-1;0;1]
    @assert p0.probs[1] == p0.probs[3]

    MandQ = collect([(Q,M) for Q in collect(range(0,stop = N,step = 1)) for M in collect(range(-Q, stop = Q, step = 2))] )

    @assert length(MandQ) == Int((N+1)*(N+2)/2)
    HomogeneousSpin1(N,p0,MandQ,
                        Dict((Q,M)=> (p0.probs[1])^Q * p0.probs[2]^(N-Q) for (Q,M) in MandQ ),  # P0
                        Dict((Q,M)=> BigCombinatorics.Binomial(N,Q) * BigCombinatorics.Binomial(Q,div(Q+M,2))  for (Q,M) in MandQ ), # DOS
                        Dict((Q,M)=> logabsbinomial(N,Q)[1] + logabsbinomial(Q,div(Q+M,2))[1]  for (Q,M) in MandQ ) # logDOS
    )

end

function trace_spin1(C::HomogeneousSpin1,h1,h2,J)  
    #@assert J >= 0
    Z = 0.0
    m = 0.0
    m2 = 0.0
    χ = 0.0

    @assert sum(C.dos[k] for k in keys(C.dos)) == 3^C.N
    for (Q,M) in C.MandQ
        ### you still need to include (actually subtract) the ground state energy for a more stable evaluation
        w = C.P0[(Q,M)] * exp( h1*M + 0.5 * J * M^2 + (h2 - 0.5*J)*Q + C.logdos[(Q,M)])
        Z += w
        m += M * w
        m2 += Q * w
        χ += w * (M^2-Q)
    end
    return Z, m/(C.N*Z), m2/(C.N*Z), χ/(Z*C.N*(C.N-1))
end


struct TracedSpin1Integral{TI <: IntegrationMethod} <: IntegrationMethod
    X::TI
    Z::Dict{Int,Vector{FT}}
    m::Dict{Int,Vector{FT}}
    m2::Dict{Int,Vector{FT}}
    chi::Dict{Int,Vector{FT}}
end

function clean!(X::TracedSpin1Integral{TI}) where TI <: IntegrationMethod
    for x0 in keys(X.Z)
        fill!(X.Z[x0], 0.0)    
        fill!(X.m[x0], 0.0)
        fill!(X.m2[x0], 0.0)
        fill!(X.chi[x0], 0.0)
    end
end

function inner_trace(Model::TM, X::TracedSpin1Integral{TI}) where {TM<:FPModel,TI <: IntegrationMethod}
    @assert typeof(Model)==Rank1SparseCoupled
    clean!(X)
    for (ix0,x0) in enumerate(Model.π0.values)
        for (iz, z) in enumerate(X.X.points)
            Z, m, m2, χ = inner_trace(Model,z,x0)
            X.Z[x0][iz] = Z
            X.m[x0][iz] = m
            X.m2[x0][iz] = m2
            X.chi[x0][iz] = χ
        end
    end
end


################################################### DEPRECATED #####################################################
function tovecq(i,q,n)
    x=zeros(n)
    for j=1:n
        x[j]=div(i-1,q^(n-j))+1
        i-=(x[j]-1)*(q^(n-j))
    end
    return Int.(x)
end


function trace_spin1(N::Int,π::SingleVar,h1::FT,h2::FT,J1::FT)
    #trace

    m=zeros(N)
    C=zeros(N,N)
	Z=0;
    conf = zeros(Int,N)
    
    howmany = 0
    for CI in CartesianIndices(Tuple(fill(π.Nstates,N)))
        conf .= π.ξ0[collect(Iterators.flatten(CI.I))]
        howmany+=1
        p = prod(π.probs[CI.I[i]] for i=1:N ) * exp( h1*sum(conf) + h2*sum(conf .^ 2) + J1 * sum(conf[i]*conf[j] for i=1:N for j=i+1:N)  )
        m .+= conf*p
        C .+= (conf*conf') *p
        Z += p
    end

    @assert howmany==π.Nstates^N
    m ./= Z
    C ./= Z
    return Z, m[1], C[1,1], C[1,2]
    #return m,C, log(Z)
end


function uniform_trace(N::Int,π::SingleVar,h1::FT,h2::FT,J1::FT)
    #trace

    m=0.0
    m2=0.0
    χ=0.0
	Z=0.0;
    conf = zeros(Int,N)
    
    howmany = 0
    for CI in CartesianIndices(Tuple(fill(π.Nstates,N)))
        conf .= π.ξ0[collect(Iterators.flatten(CI.I))]
        howmany+=1
        p = prod(π.probs[CI.I[i]] for i=1:N ) * exp( h1*sum(conf) + h2*sum(conf .^ 2) + J1 * sum(conf[i]*conf[j] for i=1:N for j=i+1:N)  )
        m += conf[1]*p
        m2 += conf[1]^2*p
        χ += (conf[1]*conf[2])*p
        Z += p
    end

    @assert howmany==π.Nstates^N
    m /= Z
    m2 /= Z
    χ /= Z
    return Z, m,m2, χ
end