module cantilever_examples

using FinEtools
using FinEtools.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforLinear: FEMMDeforLinear
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

function neohookean_h8()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    tmag = 0.2*phun("MPa");# Magnitude of the traction
    nincr = 8
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 9
    tolerance = W / 1000
    traction_vector = [0.0, 0.0, -tmag]

    fens, fes = H8block(L, W, H, 16, 9, 9)

    l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1)
    e2 = FDataDict("node_list"=>l1, "component"=>2)
    e3 = FDataDict("node_list"=>l1, "component"=>3)

    bfes = meshboundary(fes)
    el1 = selectelem(fens, bfes, box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
        v .= time .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    t1 = FDataDict("femm"=>FEMMBase(IntegDomain(subset(bfes, el1), GaussRule(2, 2))), "traction_vector"=>fi)

    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    # Rmout = fill(0.0, 3, 3)
    # rv = vec([-0.56 -0.1361 0.35])
    # rotmat3!(Rmout, rv)
    # mcsys = CSys(Rmout)
    # femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    region1 = FDataDict("femm"=>femm)
    modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1], "traction_bcs"=>[t1], "essential_bcs"=>[e1, e2, e3])

    modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
    modeldata["maxdu_tol"] = maxdu_tol;
    modeldata["maxiter"] = maxiter;
    modeldata["line_search"]  = true;
    modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; lambdas = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,3]));
        push!(lambdas, lambda);
    end

    modeldata = nonlinearstatics(modeldata);
    @show Ux / phun("mm"), lambdas

    pl = lineplot(Ux / phun("mm"), lambdas)
    display(pl)

    vtkexportmesh("neohookean_h8.vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
end # function neohookean_h8

function neohookean_h20()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    tmag = 0.2*phun("MPa");# Magnitude of the traction
    nincr = 8
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 9
    tolerance = W / 1000
    traction_vector = [0.0, 0.0, -tmag]

    fens, fes = H20block(L, W, H, 6, 3, 3)

    l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1)
    e2 = FDataDict("node_list"=>l1, "component"=>2)
    e3 = FDataDict("node_list"=>l1, "component"=>3)

    bfes = meshboundary(fes)
    el1 = selectelem(fens, bfes, box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
        v .= time .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    t1 = FDataDict("femm"=>FEMMBase(IntegDomain(subset(bfes, el1), GaussRule(2, 2))), "traction_vector"=>fi)

    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    # Rmout = fill(0.0, 3, 3)
    # rv = vec([-0.56 -0.1361 0.35])
    # rotmat3!(Rmout, rv)
    # mcsys = CSys(Rmout)
    # femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), mcsys, m)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    region1 = FDataDict("femm"=>femm)
    modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1], "traction_bcs"=>[t1], "essential_bcs"=>[e1, e2, e3])

    modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
    modeldata["maxdu_tol"] = maxdu_tol;
    modeldata["maxiter"] = maxiter;
    modeldata["line_search"]  = true;
    modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; lambdas = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
        push!(Ux, mean(modeldata["un1"].values[movel1,3]));
        push!(lambdas, lambda);
    end

    modeldata = nonlinearstatics(modeldata);
    @show Ux / phun("mm"), lambdas

    pl = lineplot(Ux / phun("mm"), lambdas)
    display(pl)

    vtkexportmesh("neohookean_h20.vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
end # function neohookean_h8

function allrun()
    println("#####################################################")
    println("# neohookean_h8 ")
    neohookean_h8()
    println("#####################################################")
    println("# neohookean_h20 ")
    neohookean_h20()
    return true
end # function allrun

end # module cantilever_examples
