
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

module testelasticitytensor1
using FinEtools
using FinEtoolsDeforNonlinear.MatDeforNonlinearModule: totlag2currsymm!
using FinEtoolsDeforLinear.MatDeforModule: tens4symmto6x6t!
using LinearAlgebra
using Test

function test()
	delta = (I, J) -> I == J ? 1.0 : 0.0
	lambda = 3.3
	mu = 0.156
	F = [1.02 0.03 -0.04; 0.01 0.99 -0.03; -0.01 0.02 0.95]
	C = fill(0.0, 3, 3, 3, 3)
	for I in 1:3, J in 1:3, K in 1:3, L in 1:3
		C[I, J, K, L] = lambda * delta(I, J) * delta(K, L) + 
		mu * (delta(I, K) * delta(J, L) + delta(I, L) * delta(J, K))
	end
	Cm = fill(0.0, 6, 6)
	tens4symmto6x6t!(Cm, C)
	c = fill(0.0, 6, 6)
	totlag2currsymm!(c, Cm, F)
	c2 = fill(0.0, 3, 3, 3, 3)
	for i in 1:3, j in 1:3, k in 1:3, l in 1:3
		c2[i, j, k, l] = 0.0
		for I in 1:3, J in 1:3, K in 1:3, L in 1:3
			c2[i, j, k, l] += C[I, J, K, L] / det(F) * F[i, I] * F[j, J] * F[k, K] * F[l, L]
		end
	end
	c2m = fill(0.0, 6, 6)
	tens4symmto6x6t!(c2m, c2)
	@test norm(c - c2m) <= 1.0e-12
end
end
using .testelasticitytensor1
testelasticitytensor1.test()

module m3test2a
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforLinear: stressttov!
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
    stressttov!(mr, cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    

    Fn1 = [1.1 0 0.07; 0.001 0.97 0; -0.01 0 1.03]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    J = det(Fn1)
    b = Fn1 * Fn1'
    cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
    cauv = fill( 0.0 , 6)
    stressttov!(mr, cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    
    Fn1 = [1.1 -0.0333 0.07; 0.001 0.97 0; -0.01 -0.05 1.03]
    Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
    update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    J = det(Fn1)
    b = Fn1 * Fn1'
    cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
    cauv = fill( 0.0 , 6)
    stressttov!(mr, cauv, cau)
    @test norm(cauv - stress) < 1.0e-10
    
end
end
using .m3test2a
m3test2a.test()

