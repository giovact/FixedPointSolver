

Mμ(μ) = Symbol(string("M",μ))


struct RetrievalRBMBinaryBinary
end


function Signaldot(Signal::RetrievalRBMBinaryBinary,avg::SingleVar,O::NamedVec,β::FT,j::Int) 
    return tanh(β*O[:M]) * avg.values[j]
end
Signaldot(Signal::RetrievalRBMBinaryBinary,avg::SumVar,O::NamedVec,β::FT,j::Int) = tanh(β*O[:M]) * avg.values[j]
Signaldot(Signal::RetrievalRBMBinaryBinary,avg::VectorVar,O::NamedVec,β::FT,j::Int) = sum(tanh(β*O[Mμ(μ)]) * avg.values[j,μ] for μ=1:avg.L)

function fill_avg!(Signal::RetrievalRBMBinaryBinary,avg::RT,O::NamedVec,β::FT) where RT <: DiscreteExpectation
    for j = 1 : avg.Nstates
        avg.store[j] = Signaldot(Signal,avg,O,β,j)
    end
end



energySignal(Signal::RetrievalRBMBinaryBinary,avg::SingleVar,O::NamedVec) = 0.5 * O[:M]^2 
energySignal(Signal::RetrievalRBMBinaryBinary,avg::SumVar,O::NamedVec) =  0.5 * avg.L * O[:M]^2 
energySignal(Signal::RetrievalRBMBinaryBinary,avg::VectorVar,O::NamedVec) = 0.5 * sum(avg.α[μ] * O[Mμ(μ)]^2 for μ =1:avg.L )


n_magnetizations(Signal::RetrievalRBMBinaryBinary,avg::SingleVar) =  1
n_magnetizations(Signal::RetrievalRBMBinaryBinary,avg::SumVar) = 1
n_magnetizations(Signal::RetrievalRBMBinaryBinary,avg::VectorVar) = avg.L


