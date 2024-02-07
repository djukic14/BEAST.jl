

mutable struct DoubleLayerRotatedMW3D{T,K} <: MaxwellOperator3D{T,K}
    "im times the wavenumber"
    alpha::T
    gamma::K
end

# defaultquadstrat(op::DoubleLayerRotatedMW3D, tfs::Space, bfs::Space) = DoubleNumSauterQstrat(6,7,5,5,4,3)
defaultquadstrat(op::DoubleLayerRotatedMW3D, tfs::Space, bfs::Space) = DoubleNumQStrat(6, 7)

LinearAlgebra.cross(::NormalVector, a::MWDoubleLayer3D) = DoubleLayerRotatedMW3D(a.alpha, a.gamma)

function (igd::Integrand{<:DoubleLayerRotatedMW3D})(x, y, f, g)

    nx = normal(x)

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    γ = gamma(igd.operator)
    G = exp(-γ * R) / (4π * R)
    K = -(γ + iR) * G * (iR * r)

    fvalue = getvalue(f)
    gvalue = getvalue(g)

    Kg = cross.(Ref(K), gvalue)
    nxKg = cross.(Ref(nx), Kg)
    return _krondot(fvalue, nxKg)
end

struct MWDoubleLayerRotatedFarField3D{K,U} <: MWFarField
    gamma::K
    waveimpedance::U
end


function MWDoubleLayerRotatedFarField3D(;
    gamma=nothing,
    wavenumber=nothing,
    waveimpedance=nothing
)
    if (gamma === nothing) && (wavenumber === nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if (gamma !== nothing) && (wavenumber !== nothing)
        error("Supply one of (not both) gamma or wavenumber")
    end

    if gamma === nothing
        if iszero(real(wavenumber))
            gamma = -imag(wavenumber)
        else
            gamma = im * wavenumber
        end
    end

    @assert gamma !== nothing

    waveimpedance === nothing && (waveimpedance = 1.0)

    return MWDoubleLayerRotatedFarField3D(gamma, waveimpedance)
end

MWDoubleLayerRotatedFarField3D(op::DoubleLayerRotatedMW3D{T,U}) where {T,U} = MWDoubleLayerRotatedFarField3D(op.alpha, 1.0)

function integrand(op::MWDoubleLayerRotatedFarField3D, krn, y, f, p)
    op.waveimpedance * (y × ((krn * f[1]) × normal(p)))
end
