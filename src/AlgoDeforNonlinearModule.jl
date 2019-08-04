"""
    AlgoDeforNonlinearModule

Module for algorithms used in nonlinear deformation models.
"""
module AlgoDeforNonlinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.AlgoBaseModule: dcheck!
import Arpack: eigs
import SparseArrays: spzeros
import LinearAlgebra: mul!
my_A_mul_B!(C, A, B) = mul!(C, A, B)
import FinEtools.FieldModule: AbstractField, ndofs, setebc!, numberdofs!, applyebc!, scattersysvec!, wipe!, copyto!, incrscattersysvec!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.FEMMBaseModule: associategeometry!, distribloads, fieldfromintegpoints, elemfieldfromintegpoints
import ..FEMMDeforNonlinearBaseModule: stiffness, geostiffness, nzebcloads, restoringforce
import FinEtoolsDeforLinear.DeforModelRedModule: stresscomponentmap
import FinEtools.ForceIntensityModule: ForceIntensity, settime!
import FinEtools.MeshModificationModule: meshboundary
import FinEtools.MeshExportModule: vtkexportmesh
import LinearAlgebra: eigen, qr, dot, cholesky, sum
# using Debugger

"""
    AlgoDeforNonlinearModule.nonlinearstatics(modeldata::FDataDict)

Algorithm for static nonlinear deformation (stress) analysis.

The algorithm chooses steps from the array of load multipliers: the step
takes it precisely from the preceding step to the next step in one go.

# Argument
`modeldata` = dictionary with values for keys

* `"fens"`  = finite element node set
* `"regions"`  = array of region dictionaries
* `"essential_bcs"` = array of essential boundary condition dictionaries
* `"traction_bcs"` = array of traction boundary condition dictionaries
* `"temperature_change"` = dictionary of data for temperature change

For each region (connected piece of the domain made of a particular
material), mandatory, the  region dictionary  contains values for keys:
* `"femm"` = finite element model machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold
  + `"displacement"` = when this key is not present, the assumption is
        that the displacement is fixed at zero (0). Otherwise, this needs to
        be set to a function with signature `f(x, lambda)`, where `x` is the
        location of the node in the reference configuration, and `lambda` is the
        load factor. In other words, whenever this quantity is supplied, it is
        implied that the displacement depends on the load factor.
  + `"component"` = which component is prescribed  (1, 2, 3)?
  + `"node_list"` = list of nodes on the boundary to which the condition
        applies (mandatory)

For traction boundary conditions (optional) each dictionary
would hold
  + `"femm"` = finite element model machine (mandatory);
  + `"traction_vector"` = traction vector, a force-intensity
    (`ForceIntensity`) object.

Control parameters
The following attributes  may be supplied:
* `"load_multipliers"` = For what load multipliers should the solution be
          calculated? Array of monotonically increasing numbers.
* `"line_search"` = Should we use line search? Boolean.  Default = true.
* `"maxdu_tol"` = Tolerance on the magnitude  of the largest incremental
          displacement component.
* `"maxbal_tol"` = Tolerance on the magnitude of the largest out-of-balance
          force component.
* `"iteration_observer"` = observer function to be called
      after each iteration.  Default is to do nothing.
      The observer function has a signature
                `iteration_observer(lambda,iter,du,modeldata)`
      where `lambda` is the current load factor, `iter` is the iteration
      number, `du` is the nodal field of current displacement increments.
* `"increment_observer"` = observer function
      to be called after convergence is reached in each step (optional)
      The observer function has a signature
                `output(lambda, modeldata)`
      where `lambda` is the current load factor. Default is to do nothing.
      The increment observer can refer to the following key-value pairs in `modeldata`:
      - "un1": converged displacement in current step
      - "un": converged displacement in the last step
      - "t": current value of load factor
      - "dt": increment of the load factor

# Output
`modeldata` = the dictionary on input is augmented with the keys
* `"geom"` = the nodal field that is the geometry
* `"u"` = the nodal field that is the computed displacement
* `"reactions"` = computed reaction nodal field
* `"timing"` = dictionary with timing results
"""
function nonlinearstatics(modeldata::FDataDict)
    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = ["fens", "regions", "essential_bcs",  "traction_bcs", "factorize"]
    essential_bcs_recognized_keys = ["displacement", "node_list", "component"]
    traction_bcs_recognized_keys = ["femm", "traction_vector"]
    regions_recognized_keys = ["femm", "body_load"]

    # Should we use line search?
    line_search = get(modeldata, "line_search", false);

    # For what load multipliers should dissolution be calculated?
    load_multipliers = get(modeldata, "load_multipliers", 1.0);

    # Is something to be done after each iteration?
    iteration_observer = get(modeldata, "iteration_observer", nothing);

    increment_observer = get(modeldata, "increment_observer", (lambda, modeldata) -> println("Load multiplier = $(lambda)"));

    # Tolerance on the magnitude  of the largest incremental displacement component
    maxdu_tol = get(modeldata, "maxdu_tol", 0.0);

    # Tolerance on the magnitude  of the out-of-balance force
    maxbal_tol = get(modeldata, "maxbal_tol", 0.0);

    # Maximum number of iterations allowed
    maxiter = get(modeldata, "maxiter", 4);

    # Extract the nodes
    fens=get(()->error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the displacement field
    un1 = NodalField(zeros(nnodes(geom), ndofs(geom)))

    modeldata["timing"] = FDataDict()

    # Apply the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing)
        for j = 1:length(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(()->error("Must get node list!"), ebc, "node_list");
            displacement = get(ebc, "displacement", nothing);
            u_fixed = zeros(FFlt, length(fenids)); # default is  zero displacement
            if (displacement != nothing) # if it is nonzero, it is dependent on the load factor
                ebc["load_factor_dependent"] = true
            end
            component = get(ebc, "component", 0); # which component?
            setebc!(un1, fenids[:], true, component, u_fixed);
        end
        applyebc!(un1);
    end

    # Number the equations
    numberdofs!(un1)
    # Converged displacement at the beginning of the current step, i. e.  u_{n}
    un = deepcopy(un1);
    # Converged displacement one step back, i. e.  u_{n-1}
    unm1 = deepcopy(un);

    # Create the necessary FEMMs
    regions = get(()->error("Must get region list!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        femm = region["femm"];
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom);
    end

    # Initially the stiffness matrix is empty.  We can check it  in the
    # treatment of nonzero essential boundary conditions and construct it if
    # needed.
    K = nothing;

    # Increment magnitudes
    load_increments = diff(unique(sort(vcat([0], vec(load_multipliers)))));

    # For all  load increments
    for incr in 1:length(load_increments) # Load-incrementation loop

        lambda = sum(load_increments[1:incr]);
        dlambda = load_increments[incr];

        # Initial value of the displacement  in the current step
        copyto!(un1, un)

        # Process load-factor-dependent essential boundary conditions:
        aany_nonzero_EBC = false;
        essential_bcs = get(modeldata, "essential_bcs", nothing);
        if (essential_bcs != nothing)
            for j = 1:length(essential_bcs)
                ebc = essential_bcs[j]
                if haskey(ebc, "load_factor_dependent")
                    fenids = get(()->error("Must get node list!"), ebc, "node_list");
                    displacement = ebc["displacement"]
                    component = get(ebc, "component", 0); # which component?
                    for k in 1:length(fenids)
                        u_fixed = displacement(geom.values[fenids[k],:], lambda);
                        setebc!(un1, [fenids[k]], true, component, u_fixed[1]);
                        aany_nonzero_EBC = aany_nonzero_EBC || (u_fixed[1] != 0.0);
                    end
                    applyebc!(un1);
                    numberdofs!(un1)
                end
            end
        end

        # Initialization
        applyebc!(un1); # Apply EBC
        numberdofs!(un1)
        du = deepcopy(un1); # Displacement increment
        applyebc!(du);

        # Initialize the load vector
        F = fill(0.0, un1.nfreedofs);
        FL = similar(F);

        # If any boundary conditions are inhomogeneous, calculate  the force
        # vector due to the displacement increment. Then update the guess of
        # the new displacement.
        if (aany_nonzero_EBC)
            du.fixed_values .= un1.fixed_values .- un.fixed_values
            for i = 1:length(regions)
                femm = deepcopy(regions[i]["femm"])
                F .+= nzebcloads(femm, geom, un, unm1, du, lambda, dlambda);
            end
            # Provided we got converged results  in the last step, we already
            # have a usable stiffness matrix.
            if K == nothing # we don't have a stiffness matrix
                K = spzeros(un1.nfreedofs, un1.nfreedofs);
                for i = 1:length(regions)
                    femm = deepcopy(regions[i]["femm"])
                    K = K + stiffness(femm, geom, un, unm1, lambda, dlambda);
                    K = K + geostiffness(femm, geom, un, unm1, lambda, dlambda);
                end
            end
            incrscattersysvec!(un1, K\F);
        end

        # Iteration loop
        iterationsuccessful = false
        iter = 1;
        while true #  Iteration loop

            # Initialize the load vector
            fill!(F, 0.0)

            # Construct the system stiffness matrix.
            K =  spzeros(un1.nfreedofs,un1.nfreedofs);
            for i = 1:length(regions)
                femm = regions[i]["femm"];
                K = K + stiffness(femm, geom, un1, un, lambda, dlambda);
                K = K + geostiffness(femm, geom, un1, un, lambda, dlambda);
            end

            # # Process the body load
            # if (isfield(model_data, 'body_load' ))
            #     for j=1:length(model_data.body_load)
            #         body_load =model_data.body_load{j};
            #         femm = femm_deformation (struct ('material',[],...
            #             'fes',body_load.fes,...
            #             'integration_rule',body_load.integration_rule));
            #         fi= force_intensity(struct('magn',body_load.force));
            #         F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 3);
            #     end
            #     clear body_load fi  femm
            # end

            # Process the traction boundary condition
            traction_bcs = get(modeldata, "traction_bcs", nothing);
            if (traction_bcs != nothing)
                for j=1:length(traction_bcs)
                    tractionbc = traction_bcs[j]
                    dcheck!(tractionbc, traction_bcs_recognized_keys)
                    femm = tractionbc["femm"]
                    fi = tractionbc["traction_vector"];
                    settime!(fi, lambda)
                    F .= F .+ distribloads(femm, geom, un1, fi, 2);
                end
            end

            # # Process the nodal force boundary condition
            # if (isfield(model_data.boundary_conditions, 'nodal_force' ))
            #     for j=1:length(model_data.boundary_conditions.nodal_force)
            #         nodal_force =model_data.boundary_conditions.nodal_force{j};
            #         femm = femm_deformation (struct ('material',[],...
            #             'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
            #             'integration_rule',point_rule));
            #         fi= force_intensity(struct('magn',nodal_force.force));
            #         F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 0);
            #     end
            #     clear nodal_force fi femm
            # end

            #  External loads vector
            copyto!(FL, F);

            # The restoring  force vector
            FR = fill(0.0, un1.nfreedofs);
            for i = 1:length(regions)
                femm = regions[i]["femm"];
                FR .+= restoringforce(femm, geom, un1, un, lambda, dlambda);
            end

            # Solve the system of linear algebraic equations
            F .= (FL .+ FR);
            dusol = K\F;

            #  Distribute the solution
            fill!(du.values, 0.0)
            scattersysvec!(du, dusol);

            # Do we need line search?
            eta = 1.0;
            if (line_search)
                un1p = deepcopy(un1)
                R0 = dot(F, dusol);
                for  linesrch in 1:1 # How many  searches should we do?
                    fill!(FR, 0.0);
                    un1p.values .= un1.values .+ eta.*du.values
                    for i = 1:length(regions)
                        femm = regions[i]["femm"];
                        FR .+= restoringforce(femm, geom, un1p, un, lambda, dlambda);
                    end
                    F .= FL .+ FR;
                    R1 = dot(F, dusol);
                    a = R0/R1;
                    if ( a < 0 )
                        eta = a/2 +sqrt((a/2)^2 -a);
                    else
                        eta = a/2;
                    end
                    eta = min(eta, 1.0);
                end
            end

            # Increment the displacements
            un1.values .= un1.values .+ eta*du.values;

            # Compute the increment data
            maxdu  = eta * maximum(abs.(du.values));
            maxbal = maximum(abs.(F));

            # Report  iteration results?
            if iteration_observer != nothing
                setindex!(modeldata, geom, "geom");
                setindex!(modeldata, un1, "un1");
                setindex!(modeldata, un, "un");
                setindex!(modeldata, lambda, "t");
                setindex!(modeldata, dlambda, "dt");
                setindex!(modeldata, maxdu, "maxdu");
                setindex!(modeldata, maxbal, "maxbal");
                iteration_observer(lambda, iter, du, modeldata);
            end

            # Converged?
            if (maxdu <= maxdu_tol)  || (maxbal <= maxbal_tol)
                iterationsuccessful = true
                break; # Breakout of the iteration loop
            end

            iter = iter+1;
            if iter > maxiter
                @error "Maximum number of iterations reached"
                break
            end

        end #  Iteration loop

        if !iterationsuccessful
            break
        end

        # Update the model data
        setindex!(modeldata, un1, "un1");
        setindex!(modeldata, un, "un");
        setindex!(modeldata, lambda, "t");
        setindex!(modeldata, dlambda, "dt");

        # Calculate the reactions
        # This field holds the current converged displacements;
        un1Allfree = deepcopy(un1);
        fill!(un1Allfree.is_fixed, false);# All nodes, all dofs are free.
        numberdofs!(un1Allfree);
        copyto!(un1Allfree.values, un1.values);
        unAllfree = deepcopy(un1Allfree);
        copyto!(unAllfree.values, un.values);
        # Since all the degrees of freedom are free, the restoring force method
        # can be used to compute the reactions. Note the negative sign: the
        # reactions are the opposite of the resisting forces of the
        # material.
        FR = fill(0.0, un1Allfree.nfreedofs);
        for i = 1:length(regions)
            femm = deepcopy(regions[i]["femm"])
            FR .-= restoringforce(femm, geom, un1Allfree, unAllfree, lambda, dlambda);
        end
        reactions = deepcopy(un1Allfree);
        scattersysvec!(reactions, FR);
        setindex!(modeldata, reactions, "reactions");

        # Converged
        # Update the FEMMs: save the material state vectors
        for i = 1:length(regions)
            femm = regions[i]["femm"];
            restoringforce(femm, geom, un1, un, lambda, dlambda, true);
        end

        # Report results
        increment_observer(lambda, modeldata);

        # Reset  the displacement field for the next load step
        copyto!(unm1, un);
        copyto!(un, un1);
        fill!(F, 0.0);

    end # Load-incrementation loop

    # Return some computed quantities: we stored most of them inside the
    # incrementation loop, now we add a reference to the final displacement
    setindex!(modeldata, modeldata["un1"], "u");

    return modeldata
end

end
