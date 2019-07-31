module tension_compression_examples

using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule: MatDeforNeohookeanAD
using FinEtoolsDeforNonlinear.MatDeforStVKModule: MatDeforStVK
using FinEtoolsDeforNonlinear.MatDeforStVKADModule: MatDeforStVKAD
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots

function neohookean_h8()
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
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,1]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
    end

    modeldata = nonlinearstatics(modeldata);
    @show Ux, Rx

    pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function neohookean_h8

function neohookeanad_h8()
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
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,1]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
    end

    modeldata = nonlinearstatics(modeldata);
    @show Ux, Rx

    pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function neohookeanad_h8

function stvk_h8()
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
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,1]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
    end

    modeldata = nonlinearstatics(modeldata);
    @show Ux, Rx

    pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function stvk_h8

function stvkad_h8()
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
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,1]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,1]));
    end

    modeldata = nonlinearstatics(modeldata);
    @show minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)

    pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function stvkad_h8

function allrun()
    println("#####################################################")
    println("# neohookean_h8 ")
    neohookean_h8()
    println("#####################################################")
    println("# neohookeanad_h8 ")
    neohookeanad_h8()
    println("#####################################################")
    println("# stvk_h8 ")
    stvk_h8()
    println("#####################################################")
    println("# stvkad_h8 ")
    stvkad_h8()
    return true
end # function allrun

end # module tension_compression_examples
