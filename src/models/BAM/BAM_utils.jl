setBAMsignal_L1(O,R::SingleVar) = NamedArray(O,[:M])
setBAMsignal_L1(O,R::SumVar) = NamedArray(O,[:M] )
setBAMsignal_L1(O,R::VectorVar) = NamedArray(O,vcat([Symbol(:M,i) for i=1:R.L]))


setBAMsignal_L2(O,R::SingleVar) = NamedArray(O,[:Mb])
setBAMsignal_L2(O,R::SumVar) = NamedArray(O,[:Mb] )
setBAMsignal_L2(O,R::VectorVar) = NamedArray(O,vcat([Symbol(:Mb,i) for i=1:R.L]))



struct RetrievalBAM
    store1::Vec
    store2::Vec
end
RetrievalBAM(avg::RT) where RT <: DiscreteExpectation = RetrievalBAM(zeros(avg.Nstates),zeros(avg.Nstates) )  

function clean!(in::RetrievalBAM)
    fill!(in.store1,0.0)
    fill!(in.store2,0.0) 
end


Signaldot_L1(Signal::RetrievalBAM,avg::SingleVar,O::NamedVec,j::Int) = O[:M] * avg.values[j]
Signaldot_L1(Signal::RetrievalBAM,avg::SumVar,O::NamedVec,j::Int) = O[:M] * avg.values[j]
Signaldot_L1(Signal::RetrievalBAM,avg::VectorVar,O::NamedVec,j::Int) = sum(O[Symbol(:M,μ)] * avg.values[j,μ] for μ=1:avg.L)


Signaldot_L2(Signal::RetrievalBAM,avg::SingleVar,O::NamedVec,j::Int) = O[:Mb] * avg.values[j]
Signaldot_L2(Signal::RetrievalBAM,avg::SumVar,O::NamedVec,j::Int) = O[:Mb] * avg.values[j]
Signaldot_L2(Signal::RetrievalBAM,avg::VectorVar,O::NamedVec,j::Int) = sum(O[Symbol(:Mb,μ)] * avg.values[j,μ] for μ=1:avg.L)


function fill_Signal_BAM_L1!(Signal::RetrievalBAM,avg::RT,O::NamedVec) where RT <: DiscreteExpectation    
    for j = 1 : avg.Nstates
        Signal.store1[j] = Signaldot_L1(Signal,avg,O,j)
    end
end

function fill_Signal_BAM_L2!(Signal::RetrievalBAM,avg::RT,O::NamedVec) where RT <: DiscreteExpectation
    for j = 1 : avg.Nstates
        Signal.store2[j] = Signaldot_L2(Signal,avg,O,j)
    end
end



energySignal(Signal::RetrievalBAM,avg::SingleVar,O::NamedVec) = O[:M] * O[:Mb]
energySignal(Signal::RetrievalBAM,avg::SumVar,O::NamedVec) =  avg.L * O[:M] * O[:Mb] 
energySignal(Signal::RetrievalBAM,avg::VectorVar,O::NamedVec) = sum(O[Symbol(:M,μ)] * O[Symbol(:Mb,μ)] for μ ∈ 1:avg.L)


n_magnetizations_singlelayer(Signal::RetrievalBAM,avg::SingleVar) =  1
n_magnetizations_singlelayer(Signal::RetrievalBAM,avg::SumVar) = 1
n_magnetizations_singlelayer(Signal::RetrievalBAM,avg::VectorVar) =  avg.L



n_magnetizations_bothlayers(Signal::RetrievalBAM,avg) = 2 * n_magnetizations_singlelayer(Signal,avg)