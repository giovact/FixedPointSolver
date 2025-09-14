abstract type Solver
end

struct FixedPoint <: Solver
    ρ::FT
end

# each function update should return the norm between Onew and Old after the update
function update!(xold,xnew,K::FixedPoint)
    r=norm(xnew - xold,Inf) ;
        xold .*= K.ρ ;
        xold .+= (1.0 - K.ρ) .* xnew;
    return r
end
