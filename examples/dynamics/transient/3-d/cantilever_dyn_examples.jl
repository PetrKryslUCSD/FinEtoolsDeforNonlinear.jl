module cantilever_dyn_examples

using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforLinear: FEMMDeforLinear, lumpedmass
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness, geostiffness, nzebcloads, restoringforce, estimatestablestep
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots

function neohookean_h8()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    mass_density = 1000.0*phun("kg/m^3")
    m = MatDeforNeohookean(mr, mass_density, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    tmag = 0.02*phun("MPa");# Magnitude of the traction
    tolerance = W / 1000
    traction_vector = [0.0, 0.0, -tmag]
    tend = 0.25e-3

    fens, fes = H8block(L, W, H, 16, 9, 9)
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3))

    l1  = selectnode(fens; box = [0,0,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    for i in 1:3
        setebc!(u, l1, true, i, 0.0)
    end
    applyebc!(u)
    numberdofs!(u)
    
    bfes = meshboundary(fes)
    el1 = selectelem(fens, bfes, box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)
    setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
        v .= 1.0 .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    eL1femm =  FEMMBase(IntegDomain(subset(bfes, el1), GaussRule(2, 2)))
    FL = distribloads(eL1femm, geom, u, fi, 2);
    
    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)
    femm = associategeometry!(femm, geom)

    Ux = FFlt[]; ts = FFlt[]
    increment_observer = (step, t, un1) -> begin
    	if rem(step, 100) == 0
    		println("step = $(step)")
    	    push!(Ux, mean(un1.values[movel1,3]));
    	    push!(ts, t);
    	end
    end

    un1 = deepcopy(u)
    un = deepcopy(u)
    vn1 = deepcopy(u)
    
    @show stabldt = estimatestablestep(femm, geom);
    M = lumpedmass(femm, geom, un1)
    invMv = [1.0 / M[idx, idx] for idx in 1:size(M, 1)] 
    
    dtn = 0.8 * stabldt
    tn = 0.0
    # Initial displacement, velocity, and acceleration.
    Un = gathersysvec(un1);
    Vn = gathersysvec(vn1);
    # The acceleration will be computed from the initial loads.
    An = deepcopy(Vn);
    Fn = deepcopy(Vn);
    # Temporary vectors
    Un1 = deepcopy(Un);
    Vn1 = deepcopy(Vn);
    An1 = deepcopy(An);
    step = 0;
    while tn < tend
        step = step + 1      # Step counter
        fill!(Fn, 0.0) # Zero out the load
        Fn .= Fn .+ FL # Add on the time-independent load vector
        # If this is the first step compute the initial acceleration.
        if step == 1
            An .= invMv .* Fn;
            # An .= M \ Fn;
        end
        # Update the displacements
        @. Un1 = Un + dtn*Vn + ((dtn^2)/2)*An;# displacement update
        scattersysvec!(un1, Un1);
        # Add the restoring forces, starting from the time-independent load.
        Fn .+= restoringforce(femm, geom, un1, un, tn, dtn, true)
        # Compute the new acceleration.
        An1 .= invMv .* Fn;
        # An1 .= M \ Fn;
        # Update the velocity
        @. Vn1 = Vn + (dtn/2) * (An+An1);
        # Bring the the displacement and velocity fields up to date
        # scattersysvec!(vn1, Vn1);
        # Switch the temporaries for the next step.
        (Un, Un1) = (Un1, Un);
        (Vn, Vn1) = (Vn1, Vn);
        (An, An1) = (An1, An);
        (un, un1) = (un1, un);
        if (tn == tend)   # are we at the end?
            break;
        end
        if (tn+dtn > tend) # If the next step would take us beyond the end,
            dtn = tend-tn; # hit the end precisely.
        end
        tn = tn+dtn;
        increment_observer(step, tn, un1);
    end

    # @show Ux / phun("mm"), ts

    pl = lineplot(ts, Ux / phun("mm"), canvas = DotCanvas)
    display(pl)

    vtkexportmesh("neohookean_h8.vtk", fens, fes; vectors = [("u", un1.values)])
    true
end # function neohookean_h8

function allrun()
    println("#####################################################")
    println("# neohookean_h8 ")
    neohookean_h8()             # cantilever_dyn_examples.neohookean_h8()
    return true
end # function allrun (setq term-suppress-hard-newline nil)

end # module cantilever_dyn_examples
