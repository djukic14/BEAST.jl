using BEAST
using Test

@test_throws ErrorException BEAST.gamma_wavenumber_handler(1.0, 2.0)

gamma, wavenumber = BEAST.gamma_wavenumber_handler(1.0, nothing)
@test gamma == 1.0

gamma, wavenumber = BEAST.gamma_wavenumber_handler(1.0f0, nothing)
@test gamma == 1.0f0

gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, 2.0)
@test gamma == 2.0*im

gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, 2.0f0)
@test gamma == 2.0f0*im

gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, nothing)
@test BEAST.isstatic(gamma) 
@test gamma == Val(0)

gamma, wavenumber = BEAST.gamma_wavenumber_handler(0.0, nothing)
@test BEAST.isstatic(gamma)
@test gamma == Val(0)

gamma, wavenumber = BEAST.gamma_wavenumber_handler(nothing, 0.0)
@test BEAST.isstatic(gamma)
@test gamma == Val(0)




