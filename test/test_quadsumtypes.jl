

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
function anyuninit(q::BEAST.FiveQuadStratsSumType{F,E,V,W,D}) where {F,E,V,W,D}
    (F == SumTypes.Uninit) && return true
    (E == SumTypes.Uninit) && return true
    (V == SumTypes.Uninit) && return true
    (W == SumTypes.Uninit) && return true
    (D == SumTypes.Uninit) && return true

    return false
end

function anyuninit(q::BEAST.FourQuadStratsSumType{F,E,V,D}) where {F,E,V,D}
    (F == SumTypes.Uninit) && return true
    (E == SumTypes.Uninit) && return true
    (V == SumTypes.Uninit) && return true
    (D == SumTypes.Uninit) && return true

    return false
end

@test anyuninit(BEAST.FiveQuadStratsSumType'.WiltonSE(1))
for qs in [BEAST.DoubleNumSauterQstrat(2, 3, 5, 5, 4, 3), BEAST.DoubleNumWiltonSauterQStrat(2, 3, 6, 7, 5, 5, 4, 3)]
    @show qs
    for op in [Maxwell3D.singlelayer(; wavenumber=1.0), Maxwell3D.doublelayer()]

        testelements, testassemblydata, trialelements, trialassemblydata, quaddata, _ = BEAST.assembleblock_primer(
            op, X, Y; quadstrat=qs
        )

        for p in eachindex(testelements)
            tcell = testelements[p]
            for q in eachindex(trialelements)
                bcell = trialelements[q]

                qr = BEAST.quadrule(
                    op, testshapes, trialshapes, p, tcell, q, bcell, quaddata, qs
                )
                @test qr isa BEAST.FiveQuadStratsSumType || qr isa BEAST.FourQuadStratsSumType
                @test !anyuninit(qr)
                @inferred BEAST.quadrule(
                    op, testshapes, trialshapes, p, tcell, q, bcell, quaddata, qs
                )
            end
        end
    end
end
