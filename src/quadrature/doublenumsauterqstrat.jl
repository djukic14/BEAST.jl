function quaddata(op::IntegralOperator,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_charts, trial_charts, qs::DoubleNumSauterQstrat)

    T = coordtype(test_charts[1])

    tqd = quadpoints(test_local_space, test_charts, (qs.outer_rule,))
    bqd = quadpoints(trial_local_space, trial_charts, (qs.inner_rule,))

    leg = (
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_vert, 0, 1)),
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_edge, 0, 1)),
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_face, 0, 1)),
        convert.(NTuple{2,T}, _legendre(qs.sauter_schwab_common_tetr, 0, 1)),
    )

    return (tpoints=tqd, bpoints=bqd, gausslegendre=leg)
end

"""
    _hitsanddmin2(τ, σ, distancetol=1.0e3)

Compute the number of hits and the minimum squared distance between two elements.

# Arguments
- `τ`: First element.
- `σ`: Second element.
- `distancetol`: The distance tolerance used to determine if two vertices are considered close.

# Returns
- `hits::Int`: The number of hits, i.e., the number of pairs of vertices that are considered close.
- `dmin2::Real`: The minimum squared distance between the vertices of the two elements.
"""
function _hitsanddmin2(τ, σ; distancetol=1.0e3)
    T = eltype(eltype(τ.vertices))
    hits = 0
    dtol = distancetol * eps(T)
    dmin2 = floatmax(T)
    for t in τ.vertices
        for s in σ.vertices
            d2 = LinearAlgebra.norm_sqr(t - s)
            dmin2 = min(dmin2, d2)
            hits += (d2 < dtol)
        end
    end
    return hits, dmin2
end

"""
    quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd, qs::DoubleNumSauterQstrat)

Compute the quadrature rule for the given integral operator, reference spaces, indices, and quadrature strategy.

# Arguments
- `op::IntegralOperator`: The integral operator.
- `g::RefSpace`: The reference space for the test function.
- `f::RefSpace`: The reference space for the trial function.
- `i`: The index of the test function.
- `τ`: Test element.
- `j`: The index of the trial function.
- `σ`: Trial element.
- `qd`: The quadrature data.
- `qs::DoubleNumSauterQstrat`: The double numerical quadrature strategy.

# Returns
- The computed quadrature rule as a `FourQuadStratsSumType`.

"""
function quadrule(op::IntegralOperator, g::RefSpace, f::RefSpace, i, τ, j, σ, qd,
    qs::DoubleNumSauterQstrat)

    hits, _ = _hitsanddmin2(τ, σ; distancetol=1.0e3)

    commonfacerule = SauterSchwabQuadrature.CommonFace(qd.gausslegendre[3])
    commonedgerule = SauterSchwabQuadrature.CommonEdge(qd.gausslegendre[2])
    commonvertexrule = SauterSchwabQuadrature.CommonVertex(qd.gausslegendre[1])
    doublequadrule = DoubleQuadRule(qd.tpoints[1, i], qd.bpoints[1, j])

    return _sumtypequadrule(
        hits, commonfacerule, commonedgerule, commonvertexrule,
        doublequadrule,
    )
end

"""
    _sumtypequadrule(hits, commonfacerule, commonedgerule, commonvertexrule, doublequadrule)

This function returns a `FourQuadStratsSumType` based on the number of `hits` provided.
The `FourQuadStratsSumType` is determined by the type of `commonfacerule`, `commonedgerule`,
`commonvertexrule`, and `doublequadrule` arguments.

## Arguments
- `hits`: The number of hits.
- `commonfacerule`: The rule for common face.
- `commonedgerule`: The rule for common edge.
- `commonvertexrule`: The rule for common vertex.
- `doublequadrule`: The rule for double quad.

## Returns
- `FourQuadStratsSumType`: The determined `FourQuadStratsSumType` based on the number of hits.

"""
@inline function _sumtypequadrule(
    hits, commonfacerule, commonedgerule, commonvertexrule, doublequadrule,
)

    sumtype = FourQuadStratsSumType{typeof(commonfacerule),typeof(commonedgerule),
        typeof(commonvertexrule),typeof(doublequadrule)}

    hits == 3 && return sumtype(FourQuadStratsSumType'.CommonFace(commonfacerule))
    hits == 2 && return sumtype(FourQuadStratsSumType'.CommonEdge(commonedgerule))
    hits == 1 && return sumtype(FourQuadStratsSumType'.CommonVertex(commonvertexrule))

    return sumtype(FourQuadStratsSumType'.DoubleQuad(doublequadrule))
end
