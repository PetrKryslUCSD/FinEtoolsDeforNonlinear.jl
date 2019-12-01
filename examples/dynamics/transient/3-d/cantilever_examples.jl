module cantilever_examples

using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
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
    tend = 1.0

    fens, fes = H8block(L, W, H, 16, 9, 9)
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3))

    l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    setebc!(u, l1, true, 1:3, 0.0)
    
    applyebc!(u)
    numberdofs!(u)
    
    bfes = meshboundary(fes)
    el1 = selectelem(fens, bfes, box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
        v .= time .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    eL1femm =  FEMMBase(IntegDomain(subset(bdryfes, el1), GaussRule(2, 2)))
    FL = distribloads(eL1femm, geom, u, fi, 2);
    
    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)
    femm = associategeometry!(femm, geom)
    
    
    iteration_observer = (lambda, iter, du, modeldata) -> begin
        @show lambda, iter, modeldata["maxdu"], modeldata["maxbal"]
    end
    
    Ux = FFlt[]; ts = FFlt[]
    increment_observer = (t, un1) -> begin
        push!(Ux, mean(un1.values[movel1,3]));
        push!(ts, t);
    end

    un1 = clone(u)
    un = clone(u)
    vn1 = clone(u)
    
    stabldt = estimatestablestep(femm, geom);

    t = 0.0
    # Initial displacement, velocity, and acceleration.
    U0 = gathersysvec(un1);
    V0 = gathersysvec(vn1);
    # The acceleration will be computed from the initial loads.
    A0 = deepcopy(V0);
    F0 = deepcopy(V0);
    # Temporary vectors
    U1 = deepcopy(U1);
    V1 = deepcopy(V0);
    A1 = deepcopy(A0);
    F1 = deepcopy(F0); 
    step = 0;
    while t < tend
        step = step + 1      # Step counter
        fill!(F0, 0.0) # Zero out the load
        F0 .= F0 .+ FL # Add on the time-independent load vector
        # If this is the first step compute the initial acceleration.
        if step == 1
            A0 = M\F0;
        end
        # Update the displacements
        @. U1 = U0 + dt*V0 + ((dt^2)/2)*A0;# displacement update
        scattersysvec!(un1, U1);
        # Add the restoring forces, starting from the time-independent load.
        F0 .+= restoringforce(femm, geom, un1, un, t, dt, savestate = true)
        # Compute the new acceleration.
        A1 = M\F0;
        # Update the velocity
        @. V1 = V0 + (dt/2) * (A0+A1);
        # Bring the the displacement and velocity fields up to date
        scattersysvec!(vn1, V1);
        # Switch the temporary vectors for the next step.
        (U0, U1) = (U1, U0);
        (V0, V1) = (V1, V0);
        (A0, A1) = (A1, A0);
        (un, un1) = (un1, un);
        if (t == tend)   # are we at the end?
            break;
        end
        if (t+dt > tend) # If the next step would take us beyond the end,
            dt = tend-t; # hit the end precisely.
        end
        t = t+dt;
        increment_observer(t,model_data);
        end
    end
    
    @show Ux / phun("mm"), lambdas

    pl = lineplot(Ux / phun("mm"), lambdas)
    display(pl)

    vtkexportmesh("neohookean_h8.vtk", fens, fes; vectors = [("u", modeldata["un1"].values)])
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
