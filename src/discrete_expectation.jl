abstract type DiscreteExpectation end

####################################################################################################################################################################
########################################################## SIngle Variable ξ #######################################################################################
####################################################################################################################################################################

# Expectation of E_ξ p(ξ) f(ξ) with ξ in finite alphabet. 
struct SingleVar <: DiscreteExpectation # you can change name in DISCRETE EXPECTATION
    L::Int
    ξ0::Vector{Int}
    p0::Vector{FT}
    values::Vector{Int} # vaues taken by the variable to be averaged over
    probs::Vector{FT} # probability taken by the variable to be averaged over
    Nstates::Int
    store::Vector{FT}
end

# as a generic rule of thumb. Instead of writing the SingleVar as a sum over probabilities, you might as well write it as you always did by analytically removing the average over well
# BETTER FOR consistency with the other two methods in the limit L -> 1

function SingleVar(ξ,p) 
    @assert sum(p)==1
    @assert length(ξ) == length(p)
    return SingleVar(1,ξ,p,ξ,p,length(ξ),zeros(2))
end
setSignal_OP(O,R::SingleVar) = NamedArray(O[1:1],[:M] )



####################################################################################################################################################################
########################################################## Vector of Variables ξ1,…,ξL #######################################################################################
####################################################################################################################################################################
### TODO: implement a rescaling for each pattern component, useful to track different modes with different strenght -> only here in VectorVar
# Expectation of E_(ξ1,…,ξL)  f(ξ1,…,ξL) ∏i p(ξμ) with ξμ in finite alphabet.  assuming ξμ i.i.d
struct VectorVar <: DiscreteExpectation
    L::Int 
    ξ0::Vector{Int}
    α::Vector{FT} # strengths, default to 1  
    p0::Vector{FT}
    values::Matrix{Int}
    probs::Vector{FT}
    Nstates::Int
    Mmin::FT
    store::Vector{FT}
end

function VectorVar(ξ,p, L;Mmin::FT = 0.9,strenghts::Vector{FT}=ones(L)) # Mmin is needed to shuffle the initial condition in [Mmin, 1.0]
    q = length(ξ)
    @assert sum(p)==1
    @assert length(ξ) == length(p)

    values = zeros(Int,q^L,L)
    probs = zeros(q^L)
    idx = Tuple(q*ones(Int,L))
    for (iii,CI) in enumerate(CartesianIndices(idx))
        values[iii,:] .= collect([ξ[CI.I[i]] for i =1:L])
        probs[iii] = prod(p[CI.I[i]] for i=1:L)
    end
    @assert isapprox(sum(probs),1.0)
    return VectorVar(L,ξ,strenghts,p,values, probs,q^L,Mmin,zeros(q^L))
end

setSignal_OP(O,R::VectorVar) = NamedArray(O[1:R.L],vcat([Mμ(i) for i=1:R.L]) )



####################################################################################################################################################################
########################################################## Sum of Variables ξ1,…,ξL #######################################################################################
####################################################################################################################################################################

# Expectation of E_(S)  f(S) p(S) with S=∑_μ ξμ  assuming ξμ i.i.d -> build directly the p(S) using multinomial theorem
struct SumVar <: DiscreteExpectation
    L::Int # number of patterns
    ξ0::Vector{Int} # value of components for each pattern -> size q
    p0::Vector{FT} # value of probabilities for each pattern  -> size q
    table_of_coefficients_list::Vector{Vector{Int}} # the size of this is given by the iterator
    multinomialcoeff_list::Vector{Int}
    sumξ_list::Vector{Int}
    #FINAL
    values::Vector{Int} # final values of sum_μ=1:L (ξ_μ) # unique values!! its size is L(q-1)+1
    probs::Vector{FT} # value of probabilities for each unique value of the sum  -> size L(q-1)+1
    Nstates::Int
    store::Vector{FT}
end


function SumVar(ξ0,p0,L)
    
    q = length(ξ0)
    @assert sum(p0)==1
    @assert length(ξ0) == length(p0)
    sumξ_list = Vector{Int}(undef,0)
    multinomialcoeff_list = Vector{Int}(undef,0)
    probs_list = Vector{FT}(undef,0)
    table_of_coefficients_list = Vector{Vector{Int}}(undef,0)
    howmany = 0
    for qtuple in eachmultinomial(q,L) # qtuple has size q -> each element counts the number of times element (k) of the pattern appears in the L-length sequence
        @assert sum(qtuple) == L
        push!(table_of_coefficients_list,collect(qtuple))
        push!(sumξ_list,sum(ξ0[k]*qtuple[k] for k=1:q ))

        multinomial_coeff = multinomial(qtuple)
        push!(multinomialcoeff_list,multinomial_coeff)
        push!(probs_list, multinomial_coeff * prod(p0 .^ qtuple))
        howmany +=1
    end
    @assert length(probs_list) == length(sumξ_list)
    values = sort(unique(sumξ_list))
    @assert length(values)==L*(q-1)+1
    probs = zeros(length(values))
    # sum probabilities for equal ξfinal, in the end it is the same.
    for (i,ξfinal) in enumerate(values)
        probs[i] = sum(probs_list[findall(x->x==ξfinal,sumξ_list)])
    end
    @assert isapprox(sum(probs),1.0)
    out = SumVar(L,ξ0,p0,table_of_coefficients_list,multinomialcoeff_list,sumξ_list,values,probs,length(probs),zeros(length(values))) # the last value was the incorrect one, vaffanculo
    L<=10 && check_compatibility(out)
    return out
end

setSignal_OP(O,R::SumVar) = NamedArray(O[1:1],[:M] )


function check_compatibility(Rsym)
    RMulti = VectorVar(Rsym.ξ0,Rsym.p0,Rsym.L)

    probs_check = zeros(length(Rsym.probs))

    for i=1:length(RMulti.probs)
        sumξ = sum(RMulti.values[i,:])
        @assert sumξ in Rsym.values
        idxv = findall(x->x==sumξ,Rsym.values)
        @assert length(idxv)==1
        idx = idxv[1]
        probs_check[idx] += RMulti.probs[i]
    end
    @assert isapprox(probs_check,Rsym.probs)
end




# dummy constructors
DummySingleVar()  = SingleVar([-1;1], [0.5;0.5])
DummySumVar() = SumVar([-1;1],[0.5;0.5],5)
DummyVectorVar() = VectorVar([-1;1],[0.5;0.5],5)

# constructors
PureStatesRetrievalBinary() = SingleVar([-1;1],[0.5;0.5])
MixedStatesRetrievalBinarySYM(L) = SumVar([-1;1],[0.5;0.5],L)
MixedStatesRetrievalBinary(L;Mmin = 0.98) = VectorVar([-1;1],[0.5;0.5],L;Mmin = Mmin)


PureStatesRetrievalSparsePatterns(p) = SingleVar([-1;0;1],[0.5*p;1-p;0.5*p])
MixedStatesRetrievalSparsePatternsBinarySYM(p,L) = SumVar([-1;0;1],[0.5*p;1-p;0.5*p],L)
MixedStatesRetrievalSparsePatternsBinary(p,L;Mmin = 0.95) = VectorVar([-1;0;1],[0.5*p;1-p;0.5*p],L,Mmin = Mmin)


"""
struct MixedStatesSym <: DiscreteExpectation
    L::Int 
    ξ::Vector{Int}
    probs::Vector{FT}
end

function MixedStatesSym(p,ξ, L) # modify to have the average w.r.t. the mean
    #println(ξ, "  ", O)
    Lm1 = L-1
    patterns = zeros(Int,Lm1,length(ξ)^Lm1)
    probabilities = zeros(length(ξ)^Lm1)
    idx = Tuple(3*ones(Int,Lm1))
    for (iii,CI) in enumerate(CartesianIndices(idx))
        patterns[:,iii] .= collect([ξ[CI.I[i]] for i =1 :Lm1])
        n0 = sum( patterns[:,iii] .== 0)
        probabilities[iii] = (1-p)^n0 * (0.5*p)^(Lm1 - n0) 
    end
    sumξ = unique(sum(patterns, dims = 1))[:]
    probξ = zeros(length(sumξ))
    for (iii,ww) in enumerate(eachcol(patterns))
        idx = findall(x -> x==sum(ww), sumξ)[1]
        probξ[idx] += probabilities[iii]
    end
    @assert isapprox(sum(probξ),1.0)
    MixedStatesSym(L,sumξ, probξ)
end

Signal_Op(O,R::MixedStatesSym) = NamedArray(O[1:1],[:M] )
Signaldot(R::MixedStatesSym,O,j) = O * R.ξ[j]

"""





struct BesselVar 
    l::FT
    W::Int
    values::Vector{Int} # vaues taken by the variable to be averaged over
    probs::Vector{FT} # probability taken by the variable to be averaged over
    Nstates::Int
end

function BesselVar(l;Wmax = 70)
    values = collect(-Wmax:1:Wmax)
    probs = collect([exp(-l)*SpecialFunctions.besseli(w,l) for w in values])
    Nstates = length(values)
    @assert isapprox(sum(probs),1.0)
    @assert all(!isnan, probs)
    return BesselVar(l,Wmax,values, probs, Nstates)
end


