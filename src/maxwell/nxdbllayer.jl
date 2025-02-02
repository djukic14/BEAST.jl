

mutable struct DoubleLayerRotatedMW3D{T,K} <: MaxwellOperator3D{T,K}
    "im times the wavenumber"
    alpha::T
    gamma::K
end

# defaultquadstrat(op::DoubleLayerRotatedMW3D, tfs::Space, bfs::Space) = DoubleNumSauterQstrat(6,7,5,5,4,3)
defaultquadstrat(op::DoubleLayerRotatedMW3D, tfs::Space, bfs::Space) = DoubleNumQStrat(6,7)

LinearAlgebra.cross(::NormalVector, a::MWDoubleLayer3D) = DoubleLayerRotatedMW3D(a.alpha, a.gamma)

function (igd::Integrand{<:DoubleLayerRotatedMW3D})(x,y,f,g)

    nx = normal(x)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    γ = gamma(igd.operator)
    G = exp(-γ*R)/(4π*R)
    K = -(γ + iR) * G * (iR * r)

    fvalue = getvalue(f)
    gvalue = getvalue(g)

    Kg = cross.(Ref(K), gvalue)
    nxKg = cross.(Ref(nx), Kg)
    return _krondot(fvalue, nxKg)
end
