
struct Spin1Half
    N::Int
    J::Matrix{FT}
    H::Vector{FT}
    χ::Matrix{FT}
    M::Vector{FT}
end

Spin1Half(N::Int) = Spin1Half(N, zeros(N,N), zeros(N), zeros(N,N), zeros(N))
Spin1Half(Model::TM) where {TM<:FPModel} = Spin1Half(Model.params[:y])
Spin1Half(N::Int, J::Matrix,H::Vector{FT}) = Spin1Half(N,FT.(J),H, zeros(N,N), zeros(N))
function traceSpin1Half(X::Spin1Half)
    #trace
    @extract X : N J H χ M
    fill!(M,0.0)
    fill!(χ,0.0)
    
    v = zeros(Int,N)
    Z = 0.0
    for ix=0:2^N-1
        v .= tovec(ix,N)
        p = exp(dot(v, (H .+ 0.5 * J*v) ) )
        M .+= v*p
        χ .+= (v*v') *p
        Z += p
    end

    M ./= Z
    χ ./= Z
    return Z
end
