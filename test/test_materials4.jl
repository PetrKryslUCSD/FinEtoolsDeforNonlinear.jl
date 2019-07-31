module tensioncompression1
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots
using Test
function test()
	mr = DeforModelRed3D
	E, nu = 7.0*phun("MPa"), 0.3
	m = MatDeforNeohookean(mr, E, nu)
	L= 6/2*phun("mm");
	H = 2/2*phun("mm");
	W = 2/2*phun("mm");
	umag = 1.5*phun("mm");# Magnitude of the displacement
	nincr = 48
	utol = 10e-7;
	graphics = ~true;
	maxdu_tol = W/1e7;
	maxiter = 5
	tolerance = W / 1000

	fens, fes = H8block(L, W, H, 2, 1, 1)

	l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	e1 = FDataDict("node_list"=>l1, "component"=>1)
	l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
	e2 = FDataDict("node_list"=>l2, "component"=>2)
	l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
	e3 = FDataDict("node_list"=>l3, "component"=>3)

	movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
	move(x, lambda) = table(lambda);
	e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

	# Rmout = fill(0.0, 3, 3)
	# rv = vec([-0.56 -0.1361 0.35])
	# rotmat3!(Rmout, rv)
	# mcsys = CSys(Rmout)
	# femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
	femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

	region1 = FDataDict("femm"=>femm)
	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

	modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
	modeldata["maxdu_tol"] = maxdu_tol;
	modeldata["maxiter"] = maxiter;
	modeldata["line_search"]  = true;
	modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
		# @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
	end
	Ux = FFlt[]; Rx = FFlt[]
	modeldata["increment_observer"] = (lambda, modeldata) -> begin
		push!(Ux, mean(modeldata["un1"].values[movel1,1]));
		push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
	end

	modeldata = nonlinearstatics(modeldata);
	# @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

	@test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈
	 [-0.0015, -6.547464398217384, 0.0014999999999999994, 2.6479612347602295]    
	# pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
	# display(pl)
	return true
end
end
using .tensioncompression1
tensioncompression1.test()

module tensioncompression2
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule: MatDeforNeohookeanAD
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots
using Test
function test()
	mr = DeforModelRed3D
	E, nu = 7.0*phun("MPa"), 0.3
	m = MatDeforNeohookeanAD(mr, E, nu)
	L= 6/2*phun("mm");
	H = 2/2*phun("mm");
	W = 2/2*phun("mm");
	umag = 1.5*phun("mm");# Magnitude of the displacement
	nincr = 48
	utol = 10e-7;
	graphics = ~true;
	maxdu_tol = W/1e7;
	maxiter = 5
	tolerance = W / 1000

	fens, fes = H8block(L, W, H, 2, 1, 1)

	l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	e1 = FDataDict("node_list"=>l1, "component"=>1)
	l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
	e2 = FDataDict("node_list"=>l2, "component"=>2)
	l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
	e3 = FDataDict("node_list"=>l3, "component"=>3)

	movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
	move(x, lambda) = table(lambda);
	e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

	# Rmout = fill(0.0, 3, 3)
	# rv = vec([-0.56 -0.1361 0.35])
	# rotmat3!(Rmout, rv)
	# mcsys = CSys(Rmout)
	# femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
	femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

	region1 = FDataDict("femm"=>femm)
	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

	modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
	modeldata["maxdu_tol"] = maxdu_tol;
	modeldata["maxiter"] = maxiter;
	modeldata["line_search"]  = true;
	modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
	    # @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
	end
	Ux = FFlt[]; Rx = FFlt[]
	modeldata["increment_observer"] = (lambda, modeldata) -> begin
	    push!(Ux, mean(modeldata["un1"].values[movel1,1]));
	    push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
	end

	modeldata = nonlinearstatics(modeldata);
	# @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

	@test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈
	 [-0.0015, -6.547464398217384, 0.0014999999999999994, 2.6479612347602295]    
	# pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
	# display(pl)
	return true
end
end
using .tensioncompression2
tensioncompression2.test()

module tensioncompression3
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforStVKModule: MatDeforStVK
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots
using Test
function test()
	mr = DeforModelRed3D
	E, nu = 7.0*phun("MPa"), 0.3
	m = MatDeforStVK(mr, E, nu)
	L= 6/2*phun("mm");
	H = 2/2*phun("mm");
	W = 2/2*phun("mm");
	umag = 1.0*phun("mm");# Magnitude of the displacement
	nincr = 48
	utol = 10e-7;
	graphics = ~true;
	maxdu_tol = W/1e7;
	maxiter = 5
	tolerance = W / 1000

	fens, fes = H8block(L, W, H, 2, 1, 1)

	l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	e1 = FDataDict("node_list"=>l1, "component"=>1)
	l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
	e2 = FDataDict("node_list"=>l2, "component"=>2)
	l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
	e3 = FDataDict("node_list"=>l3, "component"=>3)

	movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
	move(x, lambda) = table(lambda);
	e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

	# Rmout = fill(0.0, 3, 3)
	# rv = vec([-0.56 -0.1361 0.35])
	# rotmat3!(Rmout, rv)
	# mcsys = CSys(Rmout)
	# femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
	femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

	region1 = FDataDict("femm"=>femm)
	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

	modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
	modeldata["maxdu_tol"] = maxdu_tol;
	modeldata["maxiter"] = maxiter;
	modeldata["line_search"]  = true;
	modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
	    # @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
	end
	Ux = FFlt[]; Rx = FFlt[]
	modeldata["increment_observer"] = (lambda, modeldata) -> begin
	    push!(Ux, mean(modeldata["un1"].values[movel1,1]));
	    push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
	end

	modeldata = nonlinearstatics(modeldata);
	# @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

	@test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈
	 [-0.001, -1.2962962962962954, 0.0009999999999999996, 3.629629629629627]    
	# pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
	# display(pl)
	return true
end
end
using .tensioncompression3
tensioncompression3.test()

module tensioncompression4
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforStVKADModule: MatDeforStVKAD
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots
using Test
function test()
	mr = DeforModelRed3D
	E, nu = 7.0*phun("MPa"), 0.3
	m = MatDeforStVKAD(mr, E, nu)
	L= 6/2*phun("mm");
	H = 2/2*phun("mm");
	W = 2/2*phun("mm");
	umag = 1.0*phun("mm");# Magnitude of the displacement
	nincr = 48
	utol = 10e-7;
	graphics = ~true;
	maxdu_tol = W/1e7;
	maxiter = 5
	tolerance = W / 1000

	fens, fes = H8block(L, W, H, 2, 1, 1)

	l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	e1 = FDataDict("node_list"=>l1, "component"=>1)
	l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
	e2 = FDataDict("node_list"=>l2, "component"=>2)
	l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
	e3 = FDataDict("node_list"=>l3, "component"=>3)

	movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
	table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
	move(x, lambda) = table(lambda);
	e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

	# Rmout = fill(0.0, 3, 3)
	# rv = vec([-0.56 -0.1361 0.35])
	# rotmat3!(Rmout, rv)
	# mcsys = CSys(Rmout)
	# femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
	femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

	region1 = FDataDict("femm"=>femm)
	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

	modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
	modeldata["maxdu_tol"] = maxdu_tol;
	modeldata["maxiter"] = maxiter;
	modeldata["line_search"]  = true;
	modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
	    # @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
	end
	Ux = FFlt[]; Rx = FFlt[]
	modeldata["increment_observer"] = (lambda, modeldata) -> begin
	    push!(Ux, mean(modeldata["un1"].values[movel1,1]));
	    push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
	end

	modeldata = nonlinearstatics(modeldata);
	# @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

	@test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈
	 [-0.001, -1.2962962962962954, 0.0009999999999999996, 3.629629629629627]    
	# pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
	# display(pl)
	return true
end
end
using .tensioncompression4
tensioncompression4.test()


module tensioncompression5
using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule: MatDeforI1RivlinAD
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots
using Test

function test()
	mr = DeforModelRed3D
	E, nu = 7.0*phun("MPa"), 0.3
	mu = E / 2 / (1 + nu)
	c1 = mu / 2
	c2 = c1 / 2
	K = E / 3 / (1 - 2 * nu)
	m = MatDeforI1RivlinAD(mr, c1, c2, K)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    umag = 1.5*phun("mm");# Magnitude of the displacement
    nincr = 48
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 5
    tolerance = W / 1000

    fens, fes = H8block(L, W, H, 2, 1, 1)

    l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1)
    l2  = selectnode(fens; box = [-Inf,Inf,0,0,-Inf,Inf], inflate  =  tolerance)
    e2 = FDataDict("node_list"=>l2, "component"=>2)
    l3  = selectnode(fens; box = [-Inf,Inf,-Inf,Inf,0,0], inflate  =  tolerance)
    e3 = FDataDict("node_list"=>l3, "component"=>3)

    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
    move(x, lambda) = table(lambda);
    e4 = FDataDict("node_list"=>movel1, "component"=>1, "displacement"=>move)

    # Rmout = fill(0.0, 3, 3)
    # rv = vec([-0.56 -0.1361 0.35])
    # rotmat3!(Rmout, rv)
    # mcsys = CSys(Rmout)
    # femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    region1 = FDataDict("femm"=>femm)
    modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

    modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
    modeldata["maxdu_tol"] = maxdu_tol;
    modeldata["maxiter"] = maxiter;
    modeldata["line_search"]  = true;
    modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
        # @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,1]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
    end

    modeldata = nonlinearstatics(modeldata);
    @test [minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)] ≈ [-0.0015, -6.43875417356319, 0.0014999999999999994, 4.473481547305163]

    # pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    # display(pl)
end # function stvkad_h8

end
using .tensioncompression5
tensioncompression5.test()
