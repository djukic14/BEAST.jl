

using Test
using BEAST
using CompScienceMeshes
using SumTypes
λ = 1.0

m = meshsphere(λ / 4, λ / 10)
X = raviartthomas(m)
Y = buffachristiansen(m)

testshapes = refspace(X)
trialshapes = refspace(Y)
function anyuninit(q::BEAST.DoubleNumWiltonSauterQStratSumType{F,E,V,W,D}) where {F,E,V,W,D}
    (F == SumTypes.Uninit) && return true
    (E == SumTypes.Uninit) && return true
    (V == SumTypes.Uninit) && return true
    (W == SumTypes.Uninit) && return true
    (D == SumTypes.Uninit) && return true

    return false
end

@test anyuninit(BEAST.DoubleNumWiltonSauterQStratSumType'.WiltonSE(1))

for op in [Maxwell3D.singlelayer(; wavenumber=1.0), Maxwell3D.doublelayer()]
    qs = BEAST.defaultquadstrat(op, X, Y)
    @test qs isa BEAST.DoubleNumWiltonSauterQStrat

    testelements, testassemblydata, trialelements, trialassemblydata, quaddata, _ = BEAST.assembleblock_primer(
        op, X, Y; quadstrat=qs
    )

    for p in eachindex(testelements)
        tcell = testelements[p]
        for q in eachindex(trialelements)
            bcell = trialelements[q]

            @test BEAST.quadrule(
                op, testshapes, trialshapes, p, tcell, q, bcell, quaddata, qs
            ) isa BEAST.DoubleNumWiltonSauterQStratSumType
            @test !anyuninit(
                BEAST.quadrule(
                    op, testshapes, trialshapes, p, tcell, q, bcell, quaddata, qs
                ),
            )
            @inferred BEAST.quadrule(
                op, testshapes, trialshapes, p, tcell, q, bcell, quaddata, qs
            )
        end
    end
end
