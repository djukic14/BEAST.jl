quadraturebuffer(quadstrat) = (;)

function quadraturebuffer(qs::Union{TestRefinesTrialQStrat,TrialRefinesTestQStrat,NonConformingIntegralOpQStrat,CommonFaceOverlappingEdgeQStrat})
    return quadraturebuffer(qs.conforming_qstrat)
end

function quadraturebuffer(qs::NonConfTestBaryRefOfTrialQStrat)
    return quadraturebuffer(qs.conforming_qstrat)
end

function _sauterschwab_buffer(n)
    return (;
        I = Vector{Int}(undef, n),
        J = Vector{Int}(undef, n),
        K = Vector{Int}(undef, n),
        L = Vector{Int}(undef, n),
    )
end

function quadraturebuffer(::Union{DoubleNumSauterQstrat,DoubleNumWiltonSauterQStrat,SelfSauterOtherwiseDNumQStrat,CommonFaceVertexSauterCommonEdgeWiltonPostitiveDistanceNumQStrat})
    return (;
        edge = _sauterschwab_buffer(2),
        triangle = _sauterschwab_buffer(3),
        quadrilateral = _sauterschwab_buffer(4),
    )
end



function quadraturebuffer(::Union{SauterSchwabQuadrature.CommonVertex,SauterSchwabQuadrature.CommonEdge,SauterSchwabQuadrature.CommonFace})
    return _sauterschwab_buffer(3)
end

function quadraturebuffer(::Union{SauterSchwabQuadrature.CommonVertexQuad,SauterSchwabQuadrature.CommonEdgeQuad,SauterSchwabQuadrature.CommonFaceQuad})
    return _sauterschwab_buffer(4)
end
