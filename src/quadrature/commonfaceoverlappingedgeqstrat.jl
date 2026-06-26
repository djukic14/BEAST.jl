struct CommonFaceOverlappingEdgeQStrat{S} <: AbstractQuadStrat
    conforming_qstrat::S
end

function quaddata(a, X, Y, tels, bels, qs::CommonFaceOverlappingEdgeQStrat)
    return quaddata(a, X, Y, tels, bels, qs.conforming_qstrat)
end


function quadrule(a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::CommonFaceOverlappingEdgeQStrat)

    return quadrule(ReturnQuadrule(), a, 𝒳, 𝒴, i, τ, j, σ, qd, qs)
end

function quadrule(f::QuadruleCallback, a, 𝒳, 𝒴, i, τ, j, σ, qd,
    qs::CommonFaceOverlappingEdgeQStrat)

    if CompScienceMeshes.overlap(τ, σ)
        return quadrule(f, a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)
    end

    for (i,λ) in pairs(faces(τ))
        for (j,μ) in pairs(faces(σ))
            if CompScienceMeshes.overlap(λ, μ)
                return f(NonConformingTouchQRule(qs.conforming_qstrat, i, j))
    end end end

    # Either positive distance, common face, or common vertex, which can
    # be handled directly by the parent quadrature strategy
    return quadrule(f, a, 𝒳, 𝒴, i, τ, j, σ, qd, qs.conforming_qstrat)
end
