module mcmatopt1
using FinEtools
using FinEtoolsDeforNonlinear
using LinearAlgebra: Transpose, I, diagm, mul!, norm
using BenchmarkTools
using Test

function test()
	F = rand(3, 3)
	b = fill(0.0, 3, 3)
	mul!(b, F, Transpose(F))
	@test norm(b - F*F') < 1.0e-6
	true
end
end
using .mcmatopt1
mcmatopt1.test()

module mneotest1
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNonlinearModule: estimatesoundspeed
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.MatDeforNeohookeanNaiveModule: MatDeforNeohookeanNaive
using BenchmarkTools
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3

    m1 = MatDeforNeohookean(mr, E, nu)
    m2 = MatDeforNeohookeanNaive(mr, E, nu)
    # @show m
    update! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!
    
    stress = fill(0.0, 6)
    D = fill(0.0, 6, 6)
    output = FFlt[]
    statev = FFlt[]
    tn = 0.0
    dtn = 0.0
    loc = [0.0 0.0 0.0]
    label = 0
    quantity=:nothing
    Fn1 = [1.1 0 0; 0 1.0 0.02; 0.1 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]

    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    @show stress
    update!(m2, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    @show stress
    true
end
end
using .mneotest1
mneotest1.test()