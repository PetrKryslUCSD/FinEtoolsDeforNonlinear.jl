"""
Test of the neohookean material model
"""
module m1test1
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using LinearAlgebra: norm
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean(mr, E, nu)
    # @show m
    update! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!

    stress = fill(0.0, 6)
    output = FFlt[]
    statev = FFlt[]
    tn = 0.0
    dtn = 0.0
    loc = [0.0 0.0 0.0]
    label = 0
    quantity=:nothing

    Fn1 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress) == 0.0

    Fn1 = [1.01 0 0; 0 1.0 0; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[93365.8, 39786.2, 39786.2, 0.0, 0.0, 0.0]) / E < 1.0e-7

    Fn1 = [1.0 0 0; 0 1.01 0; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[39786.2, 93365.8, 39786.2, 0.0, 0.0, 0.0]) / E < 1.0e-7

    Fn1 = [1.0 0 0; 0 1.0 0.01; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[0.0, 269.231, 0.0, 0.0, 0.0, 26923.1]) / E < 1.0e-7
end
end
using .m1test1
m1test1.test()
