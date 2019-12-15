module cantilever_dyn_examples

using FinEtools
using FinEtoolsDeforLinear.DeforModelRedModule: DeforModelRed3D
using FinEtoolsDeforLinear: FEMMDeforLinear, lumpedmass
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.MatDeforNeohookeanNaiveModule: MatDeforNeohookeanNaive

using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness, geostiffness, nzebcloads, restoringforce, estimatestablestep
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearExplModule: FEMMDeforNonlinearExpl
using FinEtoolsDeforNonlinear.AssemblyModule: SysvecAssemblerOpt
using LinearAlgebra: norm
using Statistics: mean
using SparseArrays
using DelimitedFiles
using Interpolations
using UnicodePlots

function neohookean_h8()
		timing = time()

	    mr = DeforModelRed3D
	    E, nu = 7.0*phun("MPa"), 0.3
	    mass_density = 1000.0*phun("kg/m^3")
	    m = MatDeforNeohookean(mr, mass_density, E, nu)
	    L= 6/2*phun("mm");
	    H = 2/2*phun("mm");
	    W = 2/2*phun("mm");
	    tmag = 0.1*phun("MPa");# Magnitude of the traction
	    tolerance = W / 1000
	    traction_vector = [0.0, 0.0, -tmag]
	    tend = 0.00075e-3

	    fens, fes = H8block(L, W, H, 160, 80, 80)
    @info count(fens), count(fes)
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
    function setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
        v .= 1.0 .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    eL1femm =  FEMMBase(IntegDomain(subset(bfes, el1), GaussRule(2, 2)))
    FL = distribloads(eL1femm, geom, u, fi, 2);

    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    femm = FEMMDeforNonlinearExpl(mr, IntegDomain(fes, GaussRule(3, 2)), m)
    femm = associategeometry!(femm, geom)

    Ux = FFlt[]; ts = FFlt[]
    function increment_observer(step, t, un1)
        if step == 1 || rem(step, 10) == 0
            println("$(step) steps")
            push!(Ux, mean(un1.values[movel1,3]));
            push!(ts, t);
        end
    end

    un1 = deepcopy(u)
    un = deepcopy(u)
    vn1 = deepcopy(u)

    stabldt = estimatestablestep(femm, geom);
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
    # Global vector assembler
    assembler = SysvecAssemblerOpt(deepcopy(Un1), length(Un1))
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
        tim = time()
        Fn .+= restoringforce(femm, assembler, geom, un1, un, tn, dtn, true)
        println("$(time() - tim)")
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

	timing = time() - timing
    
    pl = lineplot(ts, Ux / phun("mm"), canvas = DotCanvas)
    display(pl)

    vtkexportmesh("neohookean_h8.vtk", fens, fes; vectors = [("u", un1.values)])

    println("$step steps done in $(timing) seconds")
    mus = timing / count(fes) / step * 1.0e6
    println("$(round(100 * mus) / 100) microseconds per element in one time step")

    true
end # function neohookean_h8

# For threaded execution we will need a per-thread buffer:
# Each thread has its own machine, its own assembler, and its own assembled
# vector. This is to avoid memory contention.
struct ThreadBuffer
	femm # finite element machine, with private copy of finite elements and the material
	assembler # assembler: each thread needs to have its own in order not to write into one the same
end

function neohookean_h8_thr()
	timing = time()

    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    mass_density = 1000.0*phun("kg/m^3")
    m = MatDeforNeohookean(mr, mass_density, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    tmag = 0.1*phun("MPa");# Magnitude of the traction
    tolerance = W / 1000
    traction_vector = [0.0, 0.0, -tmag]
    tend = 0.00075e-3

    fens, fes = H8block(L, W, H, 80, 40, 40)
    @info "Mesh of $(count(fens)) nodes, $(count(fes)) elements"
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
    function setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0)
        v .= 1.0 .* traction_vector
        return v
    end
    fi = ForceIntensity(FFlt, length(traction_vector), setvector!, 0.0);
    eL1femm =  FEMMBase(IntegDomain(subset(bfes, el1), GaussRule(2, 2)))
    FL = distribloads(eL1femm, geom, u, fi, 2);

    movel1  = selectnode(fens; box = [L,L,-Inf,Inf,-Inf,Inf], inflate  =  tolerance)

    # Now we prepare  the assembly for threaded execution.
    @show nth = 1 #Base.Threads.nthreads()
    chunk = Int(floor(count(fes) / nth))
    threadbuffs = ThreadBuffer[];
    for th in 1:nth
    	# This thread will work on a subset of the mesh.
    	meshrange = th < nth ? (chunk*(th-1)+1:chunk*(th)+1-1) : (chunk*(th-1)+1:count(fes))
    	feschnk = deepcopy(subset(fes, collect(meshrange)))
    	material = deepcopy(m)
    	femm = FEMMDeforNonlinearExpl(mr, IntegDomain(feschnk, GaussRule(3, 2)), material)
    	femm = associategeometry!(femm, geom)
    	# It will use its own assembler in order to avoid memory contention.
    	assembler = SysvecAssemblerOpt(fill(0.0, u.nfreedofs), u.nfreedofs)
    	push!(threadbuffs, ThreadBuffer(femm, assembler));
    end

    Ux = FFlt[]; ts = FFlt[]
    function increment_observer(step, t, un1)
        if step == 1 || rem(step, 10) == 0
            println("$(step) steps")
            push!(Ux, mean(un1.values[movel1,3]));
            push!(ts, t);
        end
    end

    # Fields for the current, last displacements
    un1 = deepcopy(u)
    un = deepcopy(u)

    # Assemble the mass matrix for the entire mesh
    femm = FEMMDeforNonlinearExpl(mr, IntegDomain(fes, GaussRule(3, 2)), m)
    stabldt = estimatestablestep(femm, geom);
    M = lumpedmass(femm, geom, un1)
    invMv = [1.0 / M[idx, idx] for idx in 1:size(M, 1)] 
    femm = nothing

    # Prepare the time stepping
    dtn = 0.8 * stabldt
    tn = 0.0
    # Initial displacement, velocity, and acceleration.
    Un = gathersysvec(un1);
    Vn = gathersysvec(un1);
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
        end
        # Update the displacements
        @. Un1 = Un + dtn*Vn + ((dtn^2)/2)*An;# displacement update
        scattersysvec!(un1, Un1);
        # Add the restoring forces, starting from the time-independent load.
        tim = time()
        tim1  = time()
        tasks = [];
        for th in 1:length(threadbuffs)
        	push!(tasks, Threads.@spawn begin 
        		tim2  = time()
        		fill!(threadbuffs[th].assembler.F_buffer, 0.0)
        		# now add the restoring force from the subset of the mesh handled by this thread
        		restoringforce(threadbuffs[th].femm, threadbuffs[th].assembler, geom, un1, un, tn, dtn, true)
        		println("Thread $(Threads.threadid()): $(time() - tim2)")
        	end);
        	
        end
        println("Farm out work: $(time() - tim1)")# Wait for the threads to finish, and then add the force from the thread to the global force vector
        tim1  = time()
        for th in 1:length(tasks)
        	Threads.wait(tasks[th]);
        	Fn .+= threadbuffs[th].assembler.F_buffer
        end
        println("Collect results: $(time() - tim1)")
        println("Total: $(time() - tim)")
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

	timing = time() - timing
    
    pl = lineplot(ts, Ux / phun("mm"), canvas = DotCanvas)
    display(pl)

    vtkexportmesh("neohookean_h8_thr.vtk", fens, fes; vectors = [("u", un1.values)])

    println("$step steps done in $(timing) seconds")
    mus = timing / count(fes) / step * 1.0e6
    println("$(round(100 * mus) / 100) microseconds per element in one time step")

    true
end # function neohookean_h8_thr

function allrun()
    println("#####################################################")
    println("# neohookean_h8 ")
    neohookean_h8()             # cantilever_dyn_examples.neohookean_h8()
    println("#####################################################")
    println("# neohookean_h8_thr ")
    neohookean_h8_thr()             # cantilever_dyn_examples.neohookean_h8_thr()
    return true
end # function allrun 

end # module cantilever_dyn_examples
