module m3test2a
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforLinear: stress3x3tto6v!
using FinEtoolsDeforNonlinear
using LinearAlgebra
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    c1, c2, K = 0.55, 0.3, 55.0
    m1 = FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule.MatDeforI1RivlinAD(mr, c1, c2, K)
    # @show m
    update! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!
    tangentmoduli! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.tangentmoduli!

    stress = fill(0.0, 6)
    D = fill(0.0, 6, 6)
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
    J = det(Fn1)
    b = Fn1 * Fn1'
    cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
    cauv = fill( 0.0 , 6)
    stress3x3tto6v!(cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    

    Fn1 = [1.1 0 0.07; 0.001 0.97 0; -0.01 0 1.03]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    J = det(Fn1)
    b = Fn1 * Fn1'
    cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
    cauv = fill( 0.0 , 6)
    stress3x3tto6v!(cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    
    Fn1 = [1.1 -0.0333 0.07; 0.001 0.97 0; -0.01 -0.05 1.03]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    J = det(Fn1)
    b = Fn1 * Fn1'
    cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
    cauv = fill( 0.0 , 6)
    stress3x3tto6v!(cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    
end
end
using .m3test2a
m3test2a.test()



