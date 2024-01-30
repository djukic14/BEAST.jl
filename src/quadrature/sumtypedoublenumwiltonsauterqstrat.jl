@sum_type DoubleNumWiltonSauterQStratSumType{F,E,V,W,D} :hidden begin
    SauterSchwabCommonFace{F}(::F)
    SauterSchwabCommonEdge{E}(::E)
    SauterSchwabCommonVertex{V}(::V)
    WiltonSE{W}(::W)
    DoubleQuad{D}(::D)
end

# executes momintegrals! of the configured quadstrat in DoubleNumWiltonSauterQStratSumType
function momintegrals!(
    biop, tshs, bshs, tcell, bcell, z, strat::D
) where {D<:DoubleNumWiltonSauterQStratSumType}
    @cases strat begin
        [
            SauterSchwabCommonFace,
            SauterSchwabCommonEdge,
            SauterSchwabCommonVertex,
            WiltonSE,
            DoubleQuad,
        ](
            quadstrat
        ) => momintegrals!(biop, tshs, bshs, tcell, bcell, z, quadstrat)
    end
end

# executes momintegrals_test_refines_trial! of the configured quadstrat in
# DoubleNumWiltonSauterQStratSumType
function momintegrals_test_refines_trial!(
    out,
    op,
    test_functions,
    test_cell,
    test_chart,
    trial_functions,
    trial_cell,
    trial_chart,
    quadrule::D,
    quadstrat,
) where {D<:DoubleNumWiltonSauterQStratSumType}
    @cases quadrule begin
        [
            SauterSchwabCommonFace,
            SauterSchwabCommonEdge,
            SauterSchwabCommonVertex,
            WiltonSE,
            DoubleQuad,
        ](
            qr
        ) => momintegrals_test_refines_trial!(
            out,
            op,
            test_functions,
            test_cell,
            test_chart,
            trial_functions,
            trial_cell,
            trial_chart,
            qr,
            quadstrat,
        )
    end
end

# executes momintegrals_trial_refines_test! of the configured quadstrat in
# DoubleNumWiltonSauterQStratSumType
function momintegrals_trial_refines_test!(
    out,
    op,
    test_functions,
    test_cell,
    test_chart,
    trial_functions,
    trial_cell,
    trial_chart,
    qr::D,
    quadstrat,
) where {D<:DoubleNumWiltonSauterQStratSumType}
    @cases qr begin
        [
            SauterSchwabCommonFace,
            SauterSchwabCommonEdge,
            SauterSchwabCommonVertex,
            WiltonSE,
            DoubleQuad,
        ](
            qr
        ) => momintegrals_trial_refines_test!(
            out,
            op,
            test_functions,
            test_cell,
            test_chart,
            trial_functions,
            trial_cell,
            trial_chart,
            qr,
            quadstrat,
        )
    end
end
