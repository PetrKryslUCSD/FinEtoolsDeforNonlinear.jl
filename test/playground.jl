module m7test13a
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
    D1 = deepcopy(D)

    update!(m2, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
    # @show stress
    @test norm(stress-[863901.0097711232, 349914.99578510894, 349914.99578510894, 0.0, 0.0, 0.0] ) / E < 1.0e-7
    tangentmoduli!(m2, D, statev, Fn1, Fn, tn, dtn, loc, label)
    D2 = deepcopy(D)
    @test norm(D1-D2) / E < 1.0e-7

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
using .m7test13a
m7test13a.test()

struct DuN<:Number
	n::Tuple{Float64,Float64}
end

Base.:+(a::DuN, b::DuN) = DuN((a.n[1] + b.n[1], a.n[2] + b.n[2]))
Base.:-(a::DuN, b::DuN) = DuN((a.n[1] - b.n[1], a.n[2] - b.n[2]))
Base.:*(a::DuN, b::DuN) = DuN((a.n[1] * b.n[1], a.n[2]*b.n[1] + a.n[1]*b.n[2]))
Base.:/(a::DuN, b::DuN) = DuN((a.n[1] / b.n[1], (a.n[2]*b.n[1] - a.n[1]*b.n[2])/(b.n[1]^2)))
Base.promote_rule(::Type{DuN}, ::Type{<:Number})  = DuN
Base.convert(::Type{DuN}, x::Real) = DuN((x, zero(x)))

a = DuN((sqrt(2.0), 1.0))
b = DuN((3.0, 0.0))


convert(DuN, 13.0)
a + 1.0
1.0 + a


f(x) = x*x - x

f(a)

g(x) = 1.0/x

g(a)

# module m3test2a
# using FinEtools
# using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
# using FinEtoolsDeforLinear: stress3x3tto6v!
# using FinEtoolsDeforNonlinear
# using LinearAlgebra
# using Test
# function test()
#     mr = DeforModelRed3D
#     E, nu = 7.0*phun("MPa"), 0.3
#     c1, c2, K = 0.55, 0.3, 55.0
#     m1 = FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule.MatDeforI1RivlinAD(mr, c1, c2, K)
#     # @show m
#     update! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.update!
#     tangentmoduli! = FinEtoolsDeforNonlinear.MatDeforNonlinearModule.tangentmoduli!

#     stress = fill(0.0, 6)
#     D = fill(0.0, 6, 6)
#     output = FFlt[]
#     statev = FFlt[]
#     tn = 0.0
#     dtn = 0.0
#     loc = [0.0 0.0 0.0]
#     label = 0
#     quantity=:nothing

#     Fn1 = [1.1 0 0; 0 1.0 0; 0 0 1.0]
#     Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
#     update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#     # @show stress
#     J = det(Fn1)
#     b = Fn1 * Fn1'
#     cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
#     cauv = fill( 0.0 , 6)
#     stress3x3tto6v!(cauv, cau)
#     @test norm(cauv - stress) < 1.0e-10
    

#     Fn1 = [1.1 0 0.07; 0.001 0.97 0; -0.01 0 1.03]
#     Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
#     update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#     # @show stress
#     J = det(Fn1)
#     b = Fn1 * Fn1'
#     cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
#     cauv = fill( 0.0 , 6)
#     stress3x3tto6v!(cauv, cau)
#     @test norm(cauv - stress) < 1.0e-10
    
#     Fn1 = [1.1 -0.0333 0.07; 0.001 0.97 0; -0.01 -0.05 1.03]
#     Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
#     update!(m1, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#     # @show stress
#     J = det(Fn1)
#     b = Fn1 * Fn1'
#     cau = 2/J * (c1 + 2*c2*(tr(b) - 3)) * b + (K*(J - 1) - 2*c1/J) * I
#     cauv = fill( 0.0 , 6)
#     stress3x3tto6v!(cauv, cau)
#     @test norm(cauv - stress) < 1.0e-10
    
# end
# end
# using .m3test2a
# m3test2a.test()



