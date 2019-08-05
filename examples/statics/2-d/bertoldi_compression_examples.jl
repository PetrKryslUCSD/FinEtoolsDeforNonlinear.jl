module bertoldi_compression_examples

using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshImportModule: import_ABAQUS 
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed2DStrain
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

function neohookeanad_q4()
    mr = DeforModelRed2DStrain
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookeanAD(mr, E, nu)
    L= 6/2*phun("mm");
    W = 2/2*phun("mm");
    H = 2/2*phun("mm");
    umag = 34.5*phun("mm");# Magnitude of the displacement
    nincr = 48
    utol = 10e-7;
    graphics = ~true;
    maxdu_tol = W/1e7;
    maxiter = 15
    tolerance = 0.1 * phun("mm") / 1000

    output = MeshImportModule.import_ABAQUS("./bertoldi-specimen-1-mesh-0.001.inp")
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

vtkexportmesh("junk.vtk", fens, fes)
fens.xyz .+= 0.1 .* phun("mm") .* 2.0 .* (rand(size(fens.xyz)...) .- 0.5)

    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(2, 2), H), m)

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
        push!(Rx, sum(modeldata["reactions"].values[movel1,2]));
        vtkexportmesh("junk$(lambda).vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
    end

    modeldata = nonlinearstatics(modeldata);
    # @show Ux, Rx
    @show (minimum(Ux), minimum(Rx), maximum(Ux), maximum(Rx)) 

    pl = lineplot((L .+ Ux) ./ L, Rx, canvas = DotCanvas)
    display(pl)
end # function neohookeanad_q4

function allrun()
    println("#####################################################")
    println("# neohookeanad_q4 ")
    neohookeanad_q4()
    return true
end # function allrun

end # module tension_compression_examples
