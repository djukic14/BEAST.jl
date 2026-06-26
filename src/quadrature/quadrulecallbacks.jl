struct ApplyMomintegrals{Z,Op,TS,TP,TC,BS,BP,BC,Buf} <: QuadruleCallback
    zlocal::Z
    biop::Op
    test_space::TS
    tptr::TP
    tcell::TC
    trial_space::BS
    bptr::BP
    bcell::BC
    qbuffer::Buf
end

function (f::ApplyMomintegrals)(qrule)
    return momintegrals!(f.zlocal, f.biop,
        f.test_space, f.tptr, f.tcell,
        f.trial_space, f.bptr, f.bcell,
        qrule, f.qbuffer)
end

struct ApplyLocalMomintegrals{Op,TS,BS,TC,BC,Z,Buf} <: QuadruleCallback
    op::Op
    test_local_space::TS
    trial_local_space::BS
    test_chart::TC
    trial_chart::BC
    out::Z
    qbuffer::Buf
end

function (f::ApplyLocalMomintegrals)(qrule)
    return momintegrals!(f.op,
        f.test_local_space, f.trial_local_space,
        f.test_chart, f.trial_chart,
        f.out, qrule, f.qbuffer)
end
