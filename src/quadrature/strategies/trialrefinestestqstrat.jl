struct TrialRefinesTestQStrat{S}
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::TrialRefinesTestQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end

# function quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd,
#     qs::TrialRefinesTestQStrat)

#     hits = _numhits(τ, σ)
#     if hits > 0
#         return TrialRefinesTestQRule(qs.conforming_qstrat)
#     end

#     return quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)    
# end

function momintegrals!(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::TrialRefinesTestQStrat, test_space, trial_space, zlocal)

    hits = _numhits(τ, σ)
    if hits > 0
        return momintegrals!(zlocal, a, test_space, i, τ, trial_space, j, σ, 
            TrialRefinesTestQRule(qs.conforming_qstrat))  
    end

    return momintegrals!(a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat, test_space, 
        trial_space, zlocal)   
end