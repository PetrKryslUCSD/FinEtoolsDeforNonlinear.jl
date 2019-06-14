module m1testsingle4
using FinEtools
using FinEtoolsDeforNonlinear
using FinEtoolsDeforNonlinear.MatDeforNeohookeanModule: MatDeforNeohookean
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearModule: FEMMDeforNonlinear
using FinEtoolsDeforNonlinear.FEMMDeforNonlinearBaseModule: stiffness
using LinearAlgebra: norm, cross
using SparseArrays
using DelimitedFiles
# using Debugger
using Interpolations
using Test
function test()
    mr = DeforModelRed3D
    E, nu = 7.0*phun("MPa"), 0.3
    m = MatDeforNeohookean(mr, E, nu)
    L= 6/2*phun("mm");
    H = 2/2*phun("mm");
    W = 2/2*phun("mm");
    umag=-1.5*phun("mm");# Magnitude of the displacement
    nincr = 48

    fens, fes = H8hexahedron([0 0 0; 1.2 0 -0.1; 1.0 0.9 0.1; -0.05 1.1 0; 0 0 1.03; 1.2 0 1.1; 1.0 0.95 0.81; -0.05 1.01 0.9], 2, 1, 1)
    femm = FEMMDeforNonlinear(mr, IntegDomain(fes, GaussRule(3, 2)), m)

    geom = NodalField(fens.xyz)
    un1 = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    un = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # Clamped cross-section
    node_list = selectnode(fens, box = [0,0,-Inf,Inf,-Inf,Inf], inflate = W/1000);
    setebc!(un, node_list, true, 1, 0.0)
    node_list = selectnode(fens, box = [-Inf,Inf,0,0,-Inf,Inf], inflate = W/1000);
    setebc!(un, node_list, true, 2, 0.0)
    node_list = selectnode(fens, box = [-Inf,Inf,-Inf,Inf,0,0], inflate = W/1000);
    setebc!(un, node_list, true, 3, 0.0)
    table = LinearInterpolation([0, 0.25, 0.5, 0.75,1], umag*[0,-1,0,1,0])
    move(lambda) = table(lambda);
    movingl = selectnode(fens, box = [L,L,-Inf,Inf,-Inf,Inf], inflate = W/1000);
    setebc!(un, movingl, true, 1, 0.0)
    applyebc!(un)
    numberdofs!(un)
    copyto!(un1, un);

    tn, dtn = 0.0, 1.0
    load_multipliers = (0:1:nincr)/nincr*1.0;
    load_increments = diff(unique(sort(load_multipliers)));

    # For all  load increments
    for incr =1:length(load_increments) # Load-implementation loop

        lambda = sum(load_increments[1:incr];
        dlambda = load_increments[incr];

        # Initial value of the displacement  in the current step
        copyto!(un1, un);

        # Process load-factor-dependent essential boundary conditions:
        setebc!(un1, movingl, true, 1, move(lambda))
        applyebc!(un1)

        # Initialization
        applyebc!(un1) # Apply EBC

        du = deepcopy(un1); # Displacement increment
        du.values[:] = 0.0

        # Initialize the load vector
        F = fill(0.0, un1.nfreedofsyou);


    # If any boundary conditions are inhomogeneous, calculate  the force
    # vector due to the displacement increment. Then update the guess of
    # the new displacement.
    if (aany_nonzero_EBC)
        for i=1:length(model_data.region)
            F = F + nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, geom, un, unm1, un1-un, dlambda);
        end
        # Provided we got converged results  in the last step, we already
        # have a usable stiffness matrix.
        if (isempty(K)) # we don't have a stiffness matrix
            K=  sparse(un1.nfreedofs,un1.nfreedofs);
            for i=1:length(model_data.region)
                K = K + stiffness(model_data.region{i}.femm, model_data.region{i}.MA, geom, un, unm1, dlambda);
                K = K + stiffness_geo(model_data.region{i}.femm, model_data.region{i}.MA, geom, un, unm1, dlambda);
            end
        end
        un1 = un1 + scatter_sysvec(du, K\F);
    end


    # Iteration loop
    iter=1;
    while true #  Iteration loop

        # Initialize the load vector
        F =0*F;

        # Construct the system stiffness matrix.
        K=  sparse(un1.nfreedofs,un1.nfreedofs);
        for i=1:length(model_data.region)
            K = K + stiffness(model_data.region{i}.femm, model_data.region{i}.MA, geom, un1, un, dlambda);
            K = K + stiffness_geo(model_data.region{i}.femm, model_data.region{i}.MA, geom, un1, un, dlambda);
        end

        # Process the body load
        if (isfield(model_data, 'body_load' ))
            for j=1:length(model_data.body_load)
                body_load =model_data.body_load{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',body_load.fes,...
                    'integration_rule',body_load.integration_rule));
                fi= force_intensity(struct('magn',body_load.force));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 3);
            end
            clear body_load fi  femm
        end

        # Process the traction boundary condition
        if (isfield(model_data.boundary_conditions, 'traction' ))
            for j=1:length(model_data.boundary_conditions.traction)
                traction =model_data.boundary_conditions.traction{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',traction.fes,...
                    'integration_rule',traction.integration_rule));
                fi= force_intensity(struct('magn',traction.traction));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 2);
            end
            clear traction fi  femm
        end

        # Process the nodal force boundary condition
        if (isfield(model_data.boundary_conditions, 'nodal_force' ))
            for j=1:length(model_data.boundary_conditions.nodal_force)
                nodal_force =model_data.boundary_conditions.nodal_force{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
                    'integration_rule',point_rule));
                fi= force_intensity(struct('magn',nodal_force.force));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 0);
            end
            clear nodal_force fi femm
        end

        #  External loads vector
        FL = lambda * F;

        # The restoring  force vector
        FR =zeros(un1.nfreedofs,1);
        for i=1:length(model_data.region)
            FR = FR + restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1, un, dlambda);
        end

        # Solve the system of linear algebraic equations
        F = (FL + FR);
        dusol= K\F;

        #  Distribute the solution
        du = scatter_sysvec(du, dusol);

        # Do we need line search?
        eta= 1.0;
        if (line_search)
            R0 = dot(F,dusol);
            for  linesrch   =1:1 # How many  searches should we do?
                FR =zeros(un1.nfreedofs,1);
                for i=1:length(model_data.region)
                    FR = FR + restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1+eta*du, un, dlambda);
                end
                F = FL + FR;
                R1 = dot(F,dusol);
                a = R0/R1;
                if ( a<0 )
                    eta = a/2 +sqrt((a/2)^2 -a);
                else
                    eta =a/2;
                end
                eta=min( [eta, 1.0] );
            end
        end

        # Increment the displacements
        un1 = un1 + eta*du;

        # Compute the increment data
        maxdu  =eta*max(max(abs(du.values)));
        maxbal=max(abs(F));

        # Report  iteration results?
        if (~isempty(iteration_observer))
            model_data.un = un;
            model_data.un1 = un1;
            model_data.dt = dlambda;
            model_data.maxdu = maxdu;
            model_data.maxbal = maxbal;
            iteration_observer (lambda,iter,du,model_data);
        end

        # Converged?
        if (maxdu <= maxdu_tol)  || (maxbal <= maxbal_tol)
            break; # Breakout of the iteration loop
        end;

        iter=iter+1;

    end #  Iteration loop


    # Update the model data
    model_data.un = un;
    model_data.un1 = un1;
    model_data.dt = dlambda;
    #  Now the reactions
    un1Allfree.values=un1.values;# This field holds the current converged displacements;
    unAllfree=un1Allfree;
    unAllfree.values=un.values;# This field holds the current converged displacements;
    # since all the degrees of freedom are free, the restoring forces
    # can be used to compute the reactions. Note the negative sign: the
    # reactions are the opposite of the resisting forces of the
    # material.
    FR =zeros(un1Allfree.nfreedofs,1);
    for i=1:length(model_data.region)# Compute the resisting forces of the material
        FR = FR - restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1Allfree, unAllfree, dlambda);
    end
    model_data.reactions=scatter_sysvec(un1Allfree,FR);

    # Converged
    # Update the FEMMs
    for i=1:length(model_data.region)
       [~,model_data.region{i}.femm]  =...
           restoring_force(model_data.region{i}.femm, sysvec_assembler, ...
           geom, un1, un, dlambda);
    end

    # Report results
   if ~isempty(increment_observer)# report the progress
        increment_observer (lambda,model_data);
    end

    # Reset  the displacement field for the next load step
    unm1 =un;
    un = un1;
    F =0*F;

end # Load-incrementation loop

end

    # K = stiffness(femm, geom, un1, un, tn, dtn)

    # ijs = readdlm("m1testop4.csv", ',', Float64)
    # matK  = sparse(Int.(ijs[:, 1]), Int.(ijs[:, 2]), ijs[:, 3], size(K, 1), size(K, 2))
    # # @show norm(K - matK) / maximum(matK[:])
    # @test norm(K - matK) / maximum(matK[:]) < 1.0e-4
end
end
using .m1testsingle4
m1testsingle4.test()
