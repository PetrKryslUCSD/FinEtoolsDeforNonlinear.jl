module bertoldi_compression_examples

using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshImportModule: import_ABAQUS 
using FinEtools.DeforModelRedModule: DeforModelRed2DStrain
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.MatDeforNeohookeanADModule: MatDeforNeohookeanAD
using FinEtoolsDeforNonlinear.MatDeforStVKModule: MatDeforStVK
using FinEtoolsDeforNonlinear.MatDeforStVKADModule: MatDeforStVKAD
using FinEtoolsDeforNonlinear.MatDeforI1RivlinADModule: MatDeforI1RivlinAD
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.AlgoDeforNonlinearModule: nonlinearstatics
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots

filetag(v) = Float64(Int(round(v*100000))) / 100000

function neohookeanad_q4()
    mr = DeforModelRed2DStrain
    E, nu = 7.0*phun("MPa"), 0.4772727
    m = MatDeforNeohookeanAD(mr, E, nu)
    L= 101.0*phun("mm");
    W = 2/2*phun("mm");
    H = 2/2*phun("mm");
    umag = 10.1*phun("mm");# Magnitude of the displacement
    nincr = 48
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 15
    tolerance = 1.0 * phun("mm") / 10

    output = MeshImportModule.import_ABAQUS("./bertoldi-specimen-1-mesh-0.002.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz = fens.xyz[:, 1:2]
    box = boundingbox(fens.xyz)
    fens2, fes2 = mirrormesh(fens, fes, vec([-1.0 0.0]),
    	box[[1, 3]]; renumb = c -> c[end:-1:1])
    fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)
    fens2, fes2 = mirrormesh(fens, fes, vec([0.0 -1.0]),
    	box[[1, 3]]; renumb = c -> c[end:-1:1])
    fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)
    box = boundingbox(fens.xyz)
    
    l1  = selectnode(fens; box = [-Inf,Inf,box[3],box[3]], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1)
    e2 = FDataDict("node_list"=>l1, "component"=>2)
    l2  = selectnode(fens; box = [-Inf,Inf,box[4],box[4]], inflate  =  tolerance)
    e3 = FDataDict("node_list"=>l2, "component"=>1)
    table = LinearInterpolation([0, 1], umag*[0,-1])
    move(x, lambda) = table(lambda);
    movel1 = l2 
    e4 = FDataDict("node_list"=>l2, "component"=>2, "displacement"=>move)
    l3  = selectnode(fens; box = [box[1],box[1],-Inf,Inf], inflate  =  tolerance)
    l4  = selectnode(fens; box = [box[2],box[2],-Inf,Inf], inflate  =  tolerance)

	xyz = deepcopy(fens.xyz)
    vtkexportmesh("bertoldi_compression-boundary.vtk", fens, meshboundary(fes))

# Perturbed the locations of the nodes in the interior
    fens.xyz .+= 0.1 .* phun("mm") .* 2.0 .* (rand(size(fens.xyz)...) .- 0.5)
    # Restore the coordinates of the points on the outer boundary
    fens.xyz[l1, :] = xyz[l1, :]
    fens.xyz[l2, :] = xyz[l2, :]
    fens.xyz[l3, :] = xyz[l3, :]
    fens.xyz[l4, :] = xyz[l4, :]

    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(2, 2), H), m)

    region1 = FDataDict("femm"=>femm)
    modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

    modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
    modeldata["maxdu_tol"] = maxdu_tol;
    modeldata["maxiter"] = maxiter;
    modeldata["line_search"]  = true;
    modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
        @show iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
    	@show lambda
        push!(Ux, mean(modeldata["un1"].values[movel1,2]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,2]));
        vtkexportmesh("bertoldi_compression-$(filetag(lambda)).vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
    end

    modeldata = nonlinearstatics(modeldata);
    # @show Ux, Rx
    @show (minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)) 

    pl = lineplot(Ux ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function neohookeanad_q4

function neohookeanad_q8()
    mr = DeforModelRed2DStrain
    E, nu = 7.0*phun("MPa"), 0.4772727
    m = MatDeforNeohookeanAD(mr, E, nu)
    L= 101.0*phun("mm");
    W = 2/2*phun("mm");
    H = 2/2*phun("mm");
    umag = 10.1*phun("mm");# Magnitude of the displacement
    nincr = 48
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 15
    tolerance = 1.0 * phun("mm") / 10

    output = MeshImportModule.import_ABAQUS("./bertoldi-specimen-1-mesh-0.002.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz = fens.xyz[:, 1:2]
    box = boundingbox(fens.xyz)
    fens2, fes2 = mirrormesh(fens, fes, vec([-1.0 0.0]),
    	box[[1, 3]]; renumb = c -> c[end:-1:1])
    fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)
    fens2, fes2 = mirrormesh(fens, fes, vec([0.0 -1.0]),
    	box[[1, 3]]; renumb = c -> c[end:-1:1])
    fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, tolerance)
    fes = cat(fes1, fes2)
    box = boundingbox(fens.xyz)

    fens, fes = Q4toQ8(fens, fes)
    
    l1  = selectnode(fens; box = [-Inf,Inf,box[3],box[3]], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1)
    e2 = FDataDict("node_list"=>l1, "component"=>2)
    l2  = selectnode(fens; box = [-Inf,Inf,box[4],box[4]], inflate  =  tolerance)
    e3 = FDataDict("node_list"=>l2, "component"=>1)
    table = LinearInterpolation([0, 1], umag*[0,-1])
    move(x, lambda) = table(lambda);
    movel1 = l2 
    e4 = FDataDict("node_list"=>l2, "component"=>2, "displacement"=>move)
    l3  = selectnode(fens; box = [box[1],box[1],-Inf,Inf], inflate  =  tolerance)
    l4  = selectnode(fens; box = [box[2],box[2],-Inf,Inf], inflate  =  tolerance)

	xyz = deepcopy(fens.xyz)
    vtkexportmesh("bertoldi_compression-boundary.vtk", fens, meshboundary(fes))

    # Perturb the locations of the nodes in the interior
    fens.xyz .+= 0.05 .* phun("mm") .* 2.0 .* (rand(size(fens.xyz)...) .- 0.5)
    # Restore the coordinates of the points on the outer boundary
    fens.xyz[l1, :] = xyz[l1, :]
    fens.xyz[l2, :] = xyz[l2, :]
    fens.xyz[l3, :] = xyz[l3, :]
    fens.xyz[l4, :] = xyz[l4, :]
    # The mid-edge nodes should be at the midpoint
    for i in 1:count(fes)
    	k = fes.conn[i][5]
    	p = fes.conn[i][1]
    	q = fes.conn[i][2]
    	fens.xyz[k, :] .= (fens.xyz[p, :] + fens.xyz[q, :]) / 2
    	k = fes.conn[i][6]
    	p = fes.conn[i][2]
    	q = fes.conn[i][3]
    	fens.xyz[k, :] .= (fens.xyz[p, :] + fens.xyz[q, :]) / 2
    	k = fes.conn[i][7]
    	p = fes.conn[i][3]
    	q = fes.conn[i][4]
    	fens.xyz[k, :] .= (fens.xyz[p, :] + fens.xyz[q, :]) / 2
    	k = fes.conn[i][8]
    	p = fes.conn[i][4]
    	q = fes.conn[i][1]
    	fens.xyz[k, :] .= (fens.xyz[p, :] + fens.xyz[q, :]) / 2
    end

    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(2, 2), H), m)

    region1 = FDataDict("femm"=>femm)
    modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],  "essential_bcs"=>[e1, e2, e3, e4])

    modeldata["load_multipliers"] = (1:nincr)./nincr*1.0;
    modeldata["maxdu_tol"] = maxdu_tol;
    modeldata["maxiter"] = maxiter;
    modeldata["line_search"]  = true;
    modeldata["iteration_observer"] = (lambda, iter, du, modeldata) -> begin
        @show iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    Ux = FFlt[]; Rx = FFlt[]
    modeldata["increment_observer"] = (lambda, modeldata) -> begin
    	@show lambda
        push!(Ux, mean(modeldata["un1"].values[movel1,2]));
        push!(Rx, sum(modeldata["reactions"].values[movel1,2]));
        vtkexportmesh("bertoldi_compression-$(filetag(lambda)).vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
    end

    modeldata = nonlinearstatics(modeldata);
    # @show Ux, Rx
    @show (minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)) 

    pl = lineplot(Ux ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function neohookeanad_q8

function allrun()
    println("#####################################################")
    println("# neohookeanad_q4 ")
    neohookeanad_q4()
    println("#####################################################")
    println("# neohookeanad_q8 ")
    neohookeanad_q8()
    return true
end # function allrun

end # module tension_compression_examples
