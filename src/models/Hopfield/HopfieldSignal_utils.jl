

Mμ(μ) = Symbol(string("M",μ))


struct RetrievalHopfield 
end


Signaldot(Signal::RetrievalHopfield,avg::SingleVar,O::NamedVec,j::Int) = O[:M] * avg.values[j]
Signaldot(Signal::RetrievalHopfield,avg::SumVar,O::NamedVec,j::Int) = O[:M] * avg.values[j]
Signaldot(Signal::RetrievalHopfield,avg::VectorVar,O::NamedVec,j::Int) = sum(O[Mμ(μ)] * avg.values[j,μ] for μ=1:avg.L)

function fill_avg!(Signal::RetrievalHopfield,avg::RT,O::NamedVec) where RT <: DiscreteExpectation
    for j = 1 : avg.Nstates
        avg.store[j] = Signaldot(Signal,avg,O,j)
    end
end



energySignal(Signal::RetrievalHopfield,avg::SingleVar,O::NamedVec) = 0.5 * O[:M]^2 
energySignal(Signal::RetrievalHopfield,avg::SumVar,O::NamedVec) =  0.5 * avg.L * O[:M]^2 
energySignal(Signal::RetrievalHopfield,avg::VectorVar,O::NamedVec) = 0.5 * sum(avg.α[μ] * O[Mμ(μ)]^2 for μ =1:avg.L )


n_magnetizations(Signal::RetrievalHopfield,avg::SingleVar) =  1
n_magnetizations(Signal::RetrievalHopfield,avg::SumVar) = 1
n_magnetizations(Signal::RetrievalHopfield,avg::VectorVar) = avg.L

