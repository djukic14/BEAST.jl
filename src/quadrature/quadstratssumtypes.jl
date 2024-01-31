@sum_type FiveQuadStratsSumType{F,E,V,W,D} :hidden begin
    CommonFace{F}(::F)
    CommonEdge{E}(::E)
    CommonVertex{V}(::V)
    WiltonSE{W}(::W) #TODO this has to be renamed along with W
    DoubleQuad{D}(::D) #TODO this has to be renamed along wiht D
end
#TODO: consider integrating BogaertInt

@sum_type FourQuadStratsSumType{F,E,V,D} :hidden begin
    CommonFace{F}(::F)
    CommonEdge{E}(::E)
    CommonVertex{V}(::V)
    DoubleQuad{D}(::D) #TODO this has to be renamed alon with D
end

#TODO: find smarter way to do this
# executes momintegrals! of the configured quadstrat in FiveQuadStratsSumType
function momintegrals!(
    biop, tshs, bshs, tcell, bcell, z, strat::D
) where {D<:FiveQuadStratsSumType}
    @cases strat begin
        [CommonFace, CommonEdge, CommonVertex, WiltonSE, DoubleQuad](
            quadstrat
        ) => momintegrals!(
            biop, tshs, bshs, tcell, bcell, z, quadstrat
        )
    end
end

function momintegrals!(
    biop, tshs, bshs, tcell, bcell, z, strat::D
) where {D<:FourQuadStratsSumType}
    @cases strat begin
        [CommonFace, CommonEdge, CommonVertex, DoubleQuad](
            quadstrat
        ) => momintegrals!(
            biop, tshs, bshs, tcell, bcell, z, quadstrat
        )
    end
end

# executes momintegrals_test_refines_trial! of the configured quadstrat in
# FiveQuadStratsSumType
function momintegrals_test_refines_trial!(
    out, op, test_functions, test_cell, test_chart, trial_functions, trial_cell,
    trial_chart, quadrule::D, quadstrat,
) where {D<:FiveQuadStratsSumType}
    @cases quadrule begin
        [CommonFace, CommonEdge, CommonVertex, WiltonSE, DoubleQuad](
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

# executes momintegrals_test_refines_trial! of the configured quadstrat in
# FourQuadStratsSumType
function momintegrals_test_refines_trial!(
    out, op, test_functions, test_cell, test_chart, trial_functions, trial_cell, trial_chart,
    quadrule::D, quadstrat,
) where {D<:FourQuadStratsSumType}
    @cases quadrule begin
        [CommonFace, CommonEdge, CommonVertex, DoubleQuad](
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
# FiveQuadStratsSumType
function momintegrals_trial_refines_test!(
    out, op, test_functions, test_cell, test_chart, trial_functions, trial_cell, trial_chart,
    qr::D, quadstrat,
) where {D<:FiveQuadStratsSumType}
    @cases qr begin
        [CommonFace, CommonEdge, CommonVertex, WiltonSE, DoubleQuad](
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

# executes momintegrals_trial_refines_test! of the configured quadstrat in
# FourQuadStratsSumType
function momintegrals_trial_refines_test!(
    out, op, test_functions, test_cell, test_chart, trial_functions, trial_cell, trial_chart,
    qr::D, quadstrat,
) where {D<:FourQuadStratsSumType}
    @cases qr begin
        [CommonFace, CommonEdge, CommonVertex, DoubleQuad](
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
