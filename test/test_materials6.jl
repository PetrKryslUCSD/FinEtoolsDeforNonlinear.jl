module mcssf1
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNonlinearModule: estimatesoundspeed
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using BenchmarkTools
using Test

function test()
    mr = DeforModelRed3D
    mass_density, E, nu = 8000*phun("kg/m^3"), 200000.0*phun("MPa"), 0.3
    m = FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean(mr, mass_density, E, nu)
    ss = estimatesoundspeed(m)
    @test abs(ss - 5801.193511153214) / 5801.193511153214 < 1.0e-6
    true
end
end
using .mcssf1
mcssf1.test()
