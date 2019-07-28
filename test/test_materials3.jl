
module m3test1a
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using LinearAlgebra: norm
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3

    m1 = FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean(mr, E, nu)
    m2 = FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule.MatDeforNeohookeanAD(mr, E, nu)
    # @show m
    update! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!
    tangentmoduli! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.tangentmoduli!

    stress = fill(0.0, 6)
    D = fill(0.0, 6, 6)
    Dtrue = [7.866603574863347e6 3.6713286713286713e6 3.6713286713286713e6 0.0 0.0 0.0; 3.6713286713286713e6  7.866603574863347e6 3.671328671328671e6 0.0 0.0 0.0;  3.6713286713286713e6 3.671328671328671e6 7.866603574863347e6 0.0 0.0 0.0; 0.0 0.0 0.0 2.0976374517673384e6 0.0 0.0; 0.0 0.0 0.0 0.0 2.0976374517673384e6 0.0; 0.0 0.0 0.0 0.0 0.0 2.097637451767338e6]
    output = FFlt[]
    statev = FFlt[]
    tn = 0.0
    dtn = 0.0
    loc = [0.0 0.0 0.0]
    label = 0
    quantity=:nothing

    Fn1 = [1.1 0 0; 0 1.0 0; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[863901.0097711232, 349914.99578510894, 349914.99578510894, 0.0, 0.0, 0.0]) / E < 1.0e-7
    tangentmoduli!(m1, D, statev, Fn1, Fn, tn, dtn, loc, label)
    @test norm(D-Dtrue) / E < 1.0e-7

    Fn1 = [1.1 0 0; 0 1.0 0; 0 0 1.0]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m2, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[863901.0097711232, 349914.99578510894, 349914.99578510894, 0.0, 0.0, 0.0] ) / E < 1.0e-7
    tangentmoduli!(m2, D, statev, Fn1, Fn, tn, dtn, loc, label)
    @test norm(D-Dtrue) / E < 1.0e-7

    # Fn1 = [1.0 0 0; 0 1.01 0; 0 0 1.0]
    # Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    # update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # # @show stress
    # @test norm(stress-[40184.6915460777, 95648.94230769236, 40184.69154607771, 0.0, 0.0, 0.0] ) / E < 1.0e-7

    # Fn1 = [1.0 0 0; 0 1.0 0.01; 0 0 1.0]
    # Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    # update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # # @show stress
    # @test norm(stress-[201.92307692305477, 740.4317307692086, 471.1538461537944, 0.0, 0.0, 26927.78846153846]) / E < 1.0e-7
end
end
using .m3test1a
m3test1a.test()
