function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumWiltonSauterQStrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space, test_charts, (qs.outer_rule_far, qs.outer_rule_near))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule_far, qs.inner_rule_near))

    leg = (
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_vert, 0, 1)),
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_edge, 0, 1)),
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_face, 0, 1)),)

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

"""
    quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs::DoubleNumWiltonSauterQStrat)

Compute the quadrature rule for the given integral operator and reference spaces.

# Arguments
- `op::IntegralOperator`: The integral operator.
- `g::RefSpace`: The reference space for the test functions.
- `f::RefSpace`: The reference space for the trial functions.
- `i`: Test index parameter.
- `τ`: Test element.
- `j`: Trial index parameter.
- `σ`: Trial element.
- `qd`: Quadrature data.
- `qs::DoubleNumWiltonSauterQStrat`: Quadrature strategy.

# Returns
The computed quadrature rule as a `FiveQuadStratsSumType``.

"""
function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd,
    qs::DoubleNumWiltonSauterQStrat)
    hits, dmin2 = _hitsanddmin2(τ, σ)

    commonfacerule = SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    commonedgerule = SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    commonvertexrule = SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
    wiltonrule = WiltonSERule(
        qd.tpoints[2, i], DoubleQuadRule(qd.tpoints[2, i], qd.bpoints[2, j])
    )
    doublequadrule = DoubleQuadRule(qd.tpoints[1, i], qd.bpoints[1, j])

    # decide which quad rule should be used and convert to appropriate sum type
    quadrule = _sumtypequadrule(
        op,
        hits,
        σ,
        dmin2,
        commonfacerule,
        commonedgerule,
        commonvertexrule,
        wiltonrule,
        doublequadrule,
    )

    return quadrule
end


"""
    _sumtypequadrule(op, hits, σ, dmin2, commonfacerule, commonedgerule, commonvertexrule, wiltonrule, doublequadrule)

This function determines the appropriate sum type quadrature rule based on the number of hits, the element volume, and other parameters.

## Arguments
- `op`: The operator.
- `hits`: The number of hits.
- `σ`:
- `dmin2`:
- `commonfacerule`: The quadrature rule for common face.
- `commonedgerule`: The quadrature rule for common edge.
- `commonvertexrule`: The quadrature rule for common vertex.
- `wiltonrule`: The Wilton quadrature rule.
- `doublequadrule`: The double quadrature rule.

## Returns
- `FiveQuadStratsSumType`: The appropriate quadrature rule.

"""
function _sumtypequadrule(
    op, hits, σ, dmin2, commonfacerule, commonedgerule, commonvertexrule, wiltonrule,
    doublequadrule,
)
    # types of all quad rules are needed for type inference of return value of this
    # function
    sumtype = FiveQuadStratsSumType{
        typeof(commonfacerule),
        typeof(commonedgerule),
        typeof(commonvertexrule),
        typeof(wiltonrule),
        typeof(doublequadrule),
    }
    hits == 3 && return convert(
        sumtype,
        FiveQuadStratsSumType'.CommonFace(commonfacerule),
    )

    hits == 2 && return convert(
        sumtype,
        FiveQuadStratsSumType'.CommonEdge(commonedgerule),
    )

    hits == 1 && return convert(
        sumtype,
        FiveQuadStratsSumType'.CommonVertex(commonvertexrule),
    )

    h2 = volume(σ)
    xtol2 = 0.2 * 0.2
    k2 = abs2(gamma(op))
    if max(dmin2 * k2, dmin2 / 16h2) < xtol2
        return convert(sumtype, FiveQuadStratsSumType'.WiltonSE(wiltonrule))
    end

    return convert(sumtype, FiveQuadStratsSumType'.DoubleQuad(doublequadrule))
end
