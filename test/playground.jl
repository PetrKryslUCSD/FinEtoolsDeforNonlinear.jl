# module m7test13bc
# using FinEtools
# using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
# using FinEtoolsDeforNonlinear
# using LinearAlgebra: norm
# using BenchmarkTools
# using Test
# function test()
#     mr = DeforModelRed3D
#     E, nu = 7.0*phun("MPa"), 0.3

#     m1 = FinEtoolsDeforNonlinear.MatDeforNeohookeanModule.MatDeforNeohookean(mr, E, nu)
#     m2 = FinEtoolsDeforNonlinear.MatDeforNeohookeanADCModule.MatDeforNeohookeanADC(mr, E, nu)
#     m3 = FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule.MatDeforNeohookeanAD(mr, E, nu)
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
#     stress1 = deepcopy(stress)
#    	tangentmoduli!(m1, D, statev, Fn1, Fn, tn, dtn, loc, label)
#    	D1 = deepcopy(D)

#    	update!(m2, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#    	stress2 = deepcopy(stress)
#    	tangentmoduli!(m2, D, statev, Fn1, Fn, tn, dtn, loc, label)
#    	D2 = deepcopy(D)

#    	@test norm(stress1-stress2) / E < 1.0e-7
#    	@test norm(D1-D2) / E < 1.0e-7

#    	update!(m3, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#    	stress3 = deepcopy(stress)
#    	tangentmoduli!(m3, D, statev, Fn1, Fn, tn, dtn, loc, label)
#    	D3 = deepcopy(D)

#     # Fn1 = [1.0 0 0; 0 1.01 0; 0 0 1.0]
#     # Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
#     # update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#     # # @show stress
#     # @test norm(stress-[40184.6915460777, 95648.94230769236, 40184.69154607771, 0.0, 0.0, 0.0] ) / E < 1.0e-7

#     # Fn1 = [1.0 0 0; 0 1.0 0.01; 0 0 1.0]
#     # Fn = [1.0 0 0; 0 1.0 0; 0 0 1.0]
#     # update!(m, statev, stress, output, Fn1, Fn, tn, dtn, loc, label, quantity)
#     # # @show stress
#     # @test norm(stress-[201.92307692305477, 740.4317307692086, 471.1538461537944, 0.0, 0.0, 26927.78846153846]) / E < 1.0e-7
# end
# end
# using .m7test13bc
# m7test13bc.test()

module mctransf1
using FinEtools
using FinEtoolsDeforNonlinear.MatDeforNonlinearModule: totalLagrangean2current!
using LinearAlgebra: det, I, diagm, diag, norm
using BenchmarkTools
using Test
symmetrize!(a) = begin
    @inbounds for c in 1:size(a, 2)
    	@inbounds for r in c:size(a, 1)
    		a[c, r] =  a[r, c] += a[c, r]
    	end
    end
    a
end
totalLagrangean2current2!(c, C, F) = begin
    n = size(F, 1)
    c .= 0.0
    for i in 1:n
    	for j in 1:n
    		for k in 1:n
    			for l in 1:n
    				for I in 1:n
    					for J in 1:n
    						for K in 1:n
    							for L in 1:n
    								# c[i, j, k, l] += F[I, i] * F[J, j] * F[K, k] * F[L, l] * C[I, J, K, L]
    								c[i, j, k, l] += F[i, I] * F[j, J] * F[k, K] * F[l, L] * C[I, J, K, L]
    							end
    						end
    					end
    				end
    			end
    		end
    	end
    end
    c ./= det(F)
    return c
end
t3x3x3x3tom6x6!(m, t) = begin
	map = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3]
    for i in 1:size(m, 1)
    	for j in 1:size(m, 2)
    		m[i, j] = t[map[i,1], map[i,2], map[j,1], map[j,2]]
    	end
    end
    return m
end
m6x6tot3x3x3x3!(t, m) = begin
	# map = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3]
	map = [1 4 5; 4 2 6; 5 6 3]
	n = 3
	for i in 1:n
		for j in 1:n
			for k in 1:n
				for l in 1:n
					t[i, j, k, l] = m[map[i, j], map[k, l]]
				end
			end
		end
	end
    return t
end
checksymmetryt3x3x3x3(C4th) = begin
	for I in 1:3
		for J in 1:3
			for K in 1:3
				for L in 1:3
					@assert C4th[I, J, K, L] == C4th[K, L, I, J]
					@assert C4th[I, J, K, L] == C4th[J, I, K, L]
					@assert C4th[I, J, K, L] == C4th[I, J, L, K]
				end
			end
		end
	end
end
checksymmetrym6x6(m) = begin
    for i in 1:size(m, 1)
    	for j in 1:size(m, 2)
    		@assert m[i, j] == m[j, i]
    	end
    end
end
function test()
	F = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.1]
	F = rand(3, 3)
	C = rand(6,6)
	symmetrize!(C)
	C = diagm(0 => diag(C))
	C4th = fill(0.0, 3, 3, 3, 3)
	C4th = m6x6tot3x3x3x3!(C4th, C)
	checksymmetryt3x3x3x3(C4th)
	# C2 = fill(0.0, 6, 6)
	# C2 = t3x3x3x3tom6x6!(C2, C4th)
	# @show C - C2
	c1 = fill(0.0, 6, 6)
	# @btime totalLagrangean2current!($c1, $C, $F)
	totalLagrangean2current!(c1, C, F)
	c14th = fill(0.0, 3, 3, 3, 3)
	c14th = m6x6tot3x3x3x3!(c14th, c1)
	c4th = fill(0.0, 3, 3, 3, 3)
	# @btime totalLagrangean2current2!($c4th, $C4th, $F)
	totalLagrangean2current2!(c4th, C4th, F)
	# @show c14th - c4th
	c2 = fill(0.0, 6, 6)
	c2 = t3x3x3x3tom6x6!(c2, c4th)
	# @show c2 - c2'
	# checksymmetrym6x6(c2)
	@show norm(c1 - c2)
end
end
using .mctransf1
mctransf1.test()


# module tensioncompression2
# using FinEtools
# using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
# using FinEtoolsDeforNonlinear
# using FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule: MatDeforNeohookeanAD
# using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
# using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
# using LinearAlgebra: norm
# using Statistics: mean
# using SparseArrays
# using DelimitedFiles
# using Interpolations
# using UnicodePlots
# using Test
# function test()
# 	mr = DeforModelRed3D
# 	E, nu = 7.0*phun("MPa"), 0.3
# 	m = MatDeforNeohookeanAD(mr, E, nu)
# 	L= 6/2*phun("mm");
# 	H = 2/2*phun("mm");
# 	W = 2/2*phun("mm");
# 	umag = 1.5*phun("mm");# Magnitude of the displacement
# 	nincr = 48
# 	utol = 10e-7;
# 	graphics = ~true;
# 	maxdu_tol = W/1e7;
# 	maxiter = 5
# 	tolerance = W / 1000

# 	fens, fes = H8block(L, W, H, 2, 1, 1)

# 	l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
# 	e1 = FDataDict("node_list"=>l1, "component"=>1)
# 	l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
# 	e2 = FDataDict("node_list"=>l2, "component"=>2)
# 	l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
# 	e3 = FDataDict("node_list"=>l3, "component"=>3)

# 	movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
# 	table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
# 	move(x, lambda) = table(lambda);
# 	e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

# 	# Rmout = fill(0.0, 3, 3)
# 	# rv = vec([-0.56 -0.1361 0.35])
# 	# rotmat3!(Rmout, rv)
# 	# mcsys = CSys(Rmout)
# 	# femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
# 	femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

# 	region1 = FDataDict("femm"=>femm)
# 	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

# 	modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
# 	modeldata["maxdu_tol"] = maxdu_tol;
# 	modeldata["maxiter"] = maxiter;
# 	modeldata["line_search"]  = true;
# 	modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
# 	    # @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
# 	end
# 	Ux = FFlt[]; Rx = FFlt[]
# 	modeldata["increment_observer"] = (lambda, modeldata) -> begin
# 	println("lambda = $(lambda)")
# 	@show modeldata["un1"]
# 	    push!(Ux, mean(modeldata["un1"].values[movel1,1]));
# 	    push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
# 	end

# 	modeldata = nonlinearstatics(modeldata);
# 	# @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

# 	@test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈
# 	 [-0.0015, -6.547464398217384, 0.0014999999999999994, 2.6479612347602295]    
# 	# pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
# 	# display(pl)
# 	return true
# end
# end
# using .tensioncompression2
# tensioncompression2.test()
