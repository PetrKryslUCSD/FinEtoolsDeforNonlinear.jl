"""
    FEMMDeforNonlinearBaseModule

Base module for operations on interiors of domains to construct system matrices and
system vectors for nonlinear deformation models.
"""
module FEMMDeforNonlinearBaseModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, locjac!, add_nnt_ut_only!, add_gkgt_ut_only!
import FinEtools.MatModule: massdensity
import FinEtoolsDeforLinear.DeforModelRedModule: nstressstrain, nthermstrain, Blmat!, divmat, vgradmat
import FinEtoolsDeforLinear.MatDeforModule: rotstressvec!, stressvtot!
import FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, tangentmoduli!, update!
import FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
import LinearAlgebra: Transpose, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, dot, det

rotF!(Fr, F, RmTF, Rm) = begin
    At_mul_B!(RmTF, Rm, F); A_mul_B!(Fr, RmTF, Rm)
end

"""
    AbstractFEMMDeforNonlinear <: AbstractFEMMDeforLinear

Abstract type of FEMM for nonlinear deformation.
"""
abstract type AbstractFEMMDeforNonlinear <: AbstractFEMMDeforLinear end

function _buff1(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    dofnums = zeros(FInt, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradXN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    gradxmN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    return dofnums, loc, J, csmatTJ, gradXN, gradxmN
end

function _buff2(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    return D, B, DB, elmat
end

function _buff3(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    X = fill(zero(FFlt), nne, sdim)
    xn = fill(zero(FFlt), nne, sdim)
    xn1 = fill(zero(FFlt), nne, sdim)
    Un = fill(zero(FFlt), nne, sdim)
    Un1 = fill(zero(FFlt), nne, sdim)
    Fn = fill(zero(FFlt), sdim, mdim);
    Fn1 = fill(zero(FFlt), sdim, mdim);
    Fnm = fill(zero(FFlt), sdim, mdim);
    Fn1m = fill(zero(FFlt), sdim, mdim);
    RmTF = fill(zero(FFlt), sdim, mdim);
    return X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF
end

function _buff4(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elvec = fill(zero(FFlt), elmatdim);      # element matrix -- buffer
    pu = fill(zero(FFlt), elmatdim);      # element matrix -- buffer
    return B, elvec, pu
end

function _buff5(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    cauchy = fill(zero(FFlt), nstrs);      # element matrix -- buffer
    output = fill(zero(FFlt), nstrs);      # element matrix -- buffer
    sigma = fill(zero(FFlt), sdim, sdim);      # element matrix -- buffer
    return cauchy, output, sigma
end

function _buff6(self::AbstractFEMMDeforNonlinear, geom::NodalField, u::NodalField)
    fes = self.integdomain.fes
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    nstrs = nstressstrain(self.mr);  # number of stresses
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    #  Indexing vector
    idx = [(1:ndn:(nne-1)*ndn+1).+(i-1) for i in 1:sdim]; # indexing vector -- buffer
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
    # to compute c1 = gradxmN*sigma*gradxmN';
    c1 = fill(zero(FFlt), nne, nne); # strain-displacement matrix -- buffer
    sg = fill(zero(FFlt), sdim, nne); # strain-displacement matrix -- buffer
    return idx, elmat, c1, sg
end

"""
    stiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble  stiffness matrix.

# Arguments
- `assembler` = matrix assembler,
- `geom` = geometry field: reference coordinates of the nodes,
- `un1` = displacement field at time `tn1 = tn + dtn`,
- `un` = displacement field at time `tn`,
- `tn` = time in step `n`
- `dtn` = increment of time
"""
function stiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradXN, gradxmN = _buff1(self, geom, un1)
    D, B, DB, elmat = _buff2(self, geom, un1)
    X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF = _buff3(self, geom, un1)
    statev = deepcopy(self.statev) # work with copies of material state
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), un1.nfreedofs, un1.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, X, fes.conn[i]);
        gathervalues_asmat!(un, Un, fes.conn[i]);
        gathervalues_asmat!(un1, Un1, fes.conn[i]);
        @. xn = X + Un; # previous known coordinates
        @. xn1 = X + Un1; # current coordinates
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]); Rm = self.mcsys.csmat
            # At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradXN, gradNparams[j], J); # wrt global c.s.
            At_mul_B!(Fn, xn, gradXN) # Previous deformation gradient
            At_mul_B!(Fn1, xn1, gradXN) # Current deformation gradient
            rotF!(Fnm, Fn, RmTF, Rm) # wrt local c.s.
            rotF!(Fn1m, Fn1, RmTF, Rm) # wrt local c.s.
            tangentmoduli!(self.material, D, statev[i][j], Fn1m, Fnm, tn, dtn, loc, fes.label[i])
            A_mul_B!(gradxmN, (gradXN / Fn1), Rm)
            Blmat!(self.mr, B, Ns[j], gradxmN, loc, Rm); # local strain-global disp
            add_btdb_ut_only!(elmat, B, Jac*w[j]*det(Fn1), D, DB)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(un1, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::AbstractFEMMDeforNonlinear, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom, un1, un, tn, dtn);
end

"""
    nzebcloads(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, du::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:SysvecAssembler, T<:Number}

Compute and assemble load vector due to prescribed increment of displacements.

# Arguments
- `assembler` = matrix assembler,
- `geom` = geometry field: reference coordinates of the nodes,
- `un1` = displacement field at time `tn1 = tn + dtn`,
- `un` = displacement field at time `tn`,
- `du` = field of displacement increment prescribed at time `tn1 = tn + dtn`,
    The increment is stored in the fixed values of this field.
- `tn` = time in step `n`
- `dtn` = increment of time
"""
function nzebcloads(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, du::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:SysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradXN, gradxmN = _buff1(self, geom, un1)
    D, B, DB, elmat = _buff2(self, geom, un1)
    X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF = _buff3(self, geom, un1)
    B, elvec, pu = _buff4(self, geom, un1)
    statev = deepcopy(self.statev) # work with copies of material state
    startassembly!(assembler, un1.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(du, pu, fes.conn[i]);
        if norm(pu) != 0.0
            gathervalues_asmat!(geom, X, fes.conn[i]);
            gathervalues_asmat!(un, Un, fes.conn[i]);
            gathervalues_asmat!(un1, Un1, fes.conn[i]);
            @. xn = X + Un; # previous known coordinates
            @. xn1 = X + Un1; # current coordinates
            fill!(elmat,  0.0); # Initialize element matrix
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                updatecsmat!(self.mcsys, loc, J, fes.label[i]); Rm = self.mcsys.csmat
                # At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
                gradN!(fes, gradXN, gradNparams[j], J); # wrt global c.s.
                At_mul_B!(Fn, xn, gradXN) # Previous deformation gradient
                At_mul_B!(Fn1, xn1, gradXN) # Current deformation gradient
                rotF!(Fnm, Fn, RmTF, Rm) # wrt local c.s.
                rotF!(Fn1m, Fn1, RmTF, Rm) # wrt local c.s.
                tangentmoduli!(self.material, D, statev[i][j], Fn1m, Fnm, tn, dtn, loc, fes.label[i])
                A_mul_B!(gradxmN, (gradXN / Fn1), Rm)
                Blmat!(self.mr, B, Ns[j], gradxmN, loc, Rm); # local strain-global disp
                add_btdb_ut_only!(elmat, B, Jac*w[j]*det(Fn1), D, DB)
            end # Loop over quadrature points
            complete_lt!(elmat)
            elvec .= -elmat * pu
            gatherdofnums!(un1, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler, elvec, dofnums); # assemble load vector
        end
    end # Loop over elements
    return makevector!(assembler);
end

function nzebcloads(self::AbstractFEMMDeforNonlinear, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, du::NodalField{T}, tn::FFlt, dtn::FFlt) where {T<:Number}
    assembler = SysvecAssembler();
    return nzebcloads(self, assembler, geom, un1, un, du, tn, dtn);
end

"""
    restoringforce(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysvecAssembler, T<:Number}

Compute the restoring force vector.

Note: This method *UPDATES* the state of the FEMM object.  In
particular, the material state gets updated.

# Arguments
- `assembler` = vector assembler,
- `geom` = geometry field: reference coordinates of the nodes,
- `un1` = displacement field at time `tn1 = tn + dtn`,
- `un` = displacement field at time `tn`,
- `tn` = time in step `n`
- `dtn` = increment of time
"""
function restoringforce(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {A<:AbstractSysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradXN, gradxmN = _buff1(self, geom, un1)
    B, elvec, pu = _buff4(self, geom, un1)
    X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF = _buff3(self, geom, un1)
    cauchy, output, sigma = _buff5(self, geom, un1)
    statev = deepcopy(self.statev) # work with copies of material state
    if savestate
        statev = self.statev # until we are sure the state should be saved
    end
    startassembly!(assembler, un1.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, X, fes.conn[i]);
        gathervalues_asmat!(un, Un, fes.conn[i]);
        gathervalues_asmat!(un1, Un1, fes.conn[i]);
        @. xn = X + Un; # previous known coordinates
        @. xn1 = X + Un1; # current coordinates
        fill!(elvec,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]); Rm = self.mcsys.csmat
            # At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradXN, gradNparams[j], J); # wrt global c.s.
            At_mul_B!(Fn, xn, gradXN) # Previous deformation gradient
            At_mul_B!(Fn1, xn1, gradXN) # Current deformation gradient
            rotF!(Fnm, Fn, RmTF, Rm) # wrt local c.s.
            rotF!(Fn1m, Fn1, RmTF, Rm) # wrt local c.s.
            update!(self.material, statev[i][j], cauchy, output, Fn1m, Fnm, tn, dtn, loc, fes.label[i])
            A_mul_B!(gradxmN, (gradXN / Fn1), Rm)
            Blmat!(self.mr, B, Ns[j], gradxmN, loc, Rm); # local strain-global disp
            elvec = elvec - B'* (cauchy * (Jac * w[j] * det(Fn1))); # note the sign
        end # Loop over quadrature points
        gatherdofnums!(un1, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elvec, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makevector!(assembler);
end

function restoringforce(self::AbstractFEMMDeforNonlinear, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {T<:Number}
    assembler = SysvecAssembler();
    return restoringforce(self, assembler, geom, un1, un, tn, dtn, savestate);
end

"""
    geostiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}

Compute and assemble geometric stiffness matrix.

# Arguments
- `assembler` = matrix assembler,
- `geom` = geometry field: reference coordinates of the nodes,
- `un1` = displacement field at time `tn1 = tn + dtn`,
- `un` = displacement field at time `tn`,
- `tn` = time in step `n`
- `dtn` = increment of time
"""
function geostiffness(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {A<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    dofnums, loc, J, csmatTJ, gradXN, gradxmN = _buff1(self, geom, un1)
    X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF = _buff3(self, geom, un1)
    cauchy, output, sigma = _buff5(self, geom, un1)
    idx, elmat, c1, sg = _buff6(self, geom, un1)
    statev = deepcopy(self.statev) # work with copies of material state
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), un1.nfreedofs, un1.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, X, fes.conn[i]);
        gathervalues_asmat!(un, Un, fes.conn[i]);
        gathervalues_asmat!(un1, Un1, fes.conn[i]);
        @. xn = X + Un; # previous known coordinates
        @. xn1 = X + Un1; # current coordinates
        fill!(elmat,  0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]); Rm = self.mcsys.csmat
            # At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradXN, gradNparams[j], J); # wrt global c.s.
            At_mul_B!(Fn, xn, gradXN) # Previous deformation gradient
            At_mul_B!(Fn1, xn1, gradXN) # Current deformation gradient
            rotF!(Fnm, Fn, RmTF, Rm) # wrt local c.s.
            rotF!(Fn1m, Fn1, RmTF, Rm) # wrt local c.s.
            update!(self.material, statev[i][j], cauchy, output, Fn1m, Fnm, tn, dtn, loc, fes.label[i])
            A_mul_B!(gradxmN, (gradXN / Fn1), Rm)
            stressvtot!(self.material.mr, sigma, cauchy);
            fill!(c1, 0.0)
            add_gkgt_ut_only!(c1, gradxmN, Jac*w[j]*det(Fn1), sigma, sg)
            complete_lt!(c1)
            for d in 1:length(idx)
            	elmat[idx[d],idx[d]]  .+= c1;
            end
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(un1, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function geostiffness(self::AbstractFEMMDeforNonlinear, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return geostiffness(self, assembler, geom, un1, un, tn, dtn);
end

#
# """
#     inspectintegpoints(self::AbstractFEMMDeforNonlinear,
#       geom::NodalField{FFlt},  u::NodalField{T},
#       dT::NodalField{FFlt},
#       felist::FIntVec,
#       inspector::F,  idat, quantity=:Cauchy;
#       context...) where {T<:Number, F<:Function}
#
# Inspect integration point quantities.
#
# - `geom` - reference geometry field
# - `u` - displacement field
# - `dT` - temperature difference field
# - `felist` - indexes of the finite elements that are to be inspected:
#      The fes to be included are: `fes[felist]`.
# - `context`    - structure: see the update!() method of the material.
# - `inspector` - functionwith the signature
#         idat = inspector(idat, j, conn, x, out, loc);
#    where
#     `idat` - a structure or an array that the inspector may
#            use to maintain some state,  for instance minimum or maximum of
#            stress, `j` is the element number, `conn` is the element connectivity,
#            `out` is the output of the update!() method,  `loc` is the location
#            of the integration point in the *reference* configuration.
# ### Return
# The updated inspector data is returned.
# """
# function inspectintegpoints(self::FEMM, geom::NodalField{FFlt},  u::NodalField{T}, dT::NodalField{FFlt}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {FEMM<:AbstractFEMMDeforNonlinear, T<:Number, F<:Function}
#     fes = self.integdomain.fes
#     npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
#     dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix = _buffers(self, geom, u)
#     # Sort out  the output requirements
#     outputcsys = self.mcsys; # default: report the stresses in the material coord system
#     for apair in pairs(context)
#         sy, val = apair
#         if sy == :outputcsys
#             outputcsys = val
#         end
#     end
#     t= 0.0
#     dt = 0.0
#     dTe = fill(zero(FFlt), nodesperelem(fes)) # nodal temperatures -- buffer
#     ue = fill(zero(FFlt), size(elmat, 1)); # array of node displacements -- buffer
#     nne = nodesperelem(fes); # number of nodes for element
#     sdim = ndofs(geom);            # number of space dimensions
#     xe = fill(zero(FFlt), nne, sdim); # array of node coordinates -- buffer
#     qpdT = 0.0; # node temperature increment
#     qpstrain = fill(zero(FFlt), nstressstrain(self.mr), 1); # total strain -- buffer
#     qpthstrain = fill(zero(FFlt), nthermstrain(self.mr)); # thermal strain -- buffer
#     qpstress = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
#     out1 = fill(zero(FFlt), nstressstrain(self.mr)); # stress -- buffer
#     out =  fill(zero(FFlt), nstressstrain(self.mr));# output -- buffer
#     # Loop over  all the elements and all the quadrature points within them
#     for ilist = 1:length(felist) # Loop over elements
#         i = felist[ilist];
#         gathervalues_asmat!(geom, xe, fes.conn[i]);# retrieve element coords
#         gathervalues_asvec!(u, ue, fes.conn[i]);# retrieve element displacements
#         gathervalues_asvec!(dT, dTe, fes.conn[i]);# retrieve element temp. increments
#         for j = 1:npts # Loop over quadrature points
#             locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
#             Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
#             updatecsmat!(self.mcsys, loc, J, fes.label[i]);
#             At_mul_B!(csmatTJ,  self.mcsys.csmat,  J); # local Jacobian matrix
#             gradN!(fes, gradN, gradNparams[j], csmatTJ);
#             Blmat!(self.mr, B, Ns[j], gradN, loc, self.mcsys.csmat);
#             updatecsmat!(outputcsys, loc, J, fes.label[i]);
#             # Quadrature point quantities
#             A_mul_B!(qpstrain, B, ue); # strain in material coordinates
#             qpdT = dot(vec(dTe), vec(Ns[j]));# Quadrature point temperature increment
#             thermalstrain!(self.material, qpthstrain, qpdT)
#             # Material updates the state and returns the output
#             out = update!(self.material, qpstress, out, vec(qpstrain), qpthstrain, t, dt, loc, fes.label[i], quantity)
#             if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
#                 (length(out1) >= length(out)) || (out1 = zeros(length(out)))
#                 rotstressvec(self.mr, out1, out, transpose(self.mcsys.csmat))# To global coord sys
#                 rotstressvec(self.mr, out, out1, outputcsys.csmat)# To output coord sys
#             end
#             # Call the inspector
#             idat = inspector(idat, i, fes.conn[i], xe, out, loc);
#         end # Loop over quadrature points
#     end # Loop over elements
#     return idat; # return the updated inspector data
# end
#
# function inspectintegpoints(self::FEMM, geom::NodalField{FFlt},  u::NodalField{T}, felist::FIntVec, inspector::F, idat, quantity=:Cauchy; context...) where {FEMM<:AbstractFEMMDeforNonlinear, T<:Number, F<:Function}
#     dT = NodalField(fill(zero(FFlt), nnodes(geom), 1)) # zero difference in temperature
#     return inspectintegpoints(self, geom, u, dT, felist, inspector, idat, quantity; context...);
# end

end
