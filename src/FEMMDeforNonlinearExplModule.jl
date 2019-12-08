"""
    FEMMDeforNonlinearExplModule

Base module for operations on interiors of domains to construct 
system vectors for explicit-dynamics nonlinear deformation models.
"""
module FEMMDeforNonlinearExplModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import ..FEMMDeforNonlinearBaseModule: AbstractFEMMDeforNonlinear
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed2DAxisymm, nstressstrain, nthermstrain, Blmat!
import FinEtoolsDeforLinear.MatDeforModule: rotstressvec!
import FinEtools.MatModule: massdensity
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, tangentmoduli!, update!, newstate
import FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
import FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, locjac!, add_nnt_ut_only!, add_gkgt_ut_only!, add_btsigma!
import LinearAlgebra: Transpose, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, dot, det
import ..FEMMDeforNonlinearBaseModule: restoringforce

rotF!(Fr, F, RmTF, Rm) = begin
    At_mul_B!(RmTF, Rm, F); A_mul_B!(Fr, RmTF, Rm)
end

struct _Buffers
	initialized::Bool
	dofnums::FIntVec
	loc::FFltMat 
	J::FFltMat 
	csmatTJ::FFltMat 
	gradXN::FFltMat 
	gradxmN::FFltMat 
	gradxN::FFltMat 
	D::FFltMat 
	B::FFltMat 
	DB::FFltMat 
	elmat::FFltMat 
	elvec::FFltVec
	pu::FFltVec
	X::FFltMat 
	xn::FFltMat 
	xn1::FFltMat 
	Un::FFltMat 
	Un1::FFltMat 
	Fn::FFltMat 
	Fn1::FFltMat 
	Fnm::FFltMat 
	Fn1m::FFltMat 
	RmTF::FFltMat 
	cauchy::FFltVec
	output::FFltVec
	sigma::FFltMat
	idx::Array{StepRange{Int64,Int64},1}
	c1::FFltMat 
	sg::FFltMat 
end

function _makebuffers(initialized, ndn, nne, sdim, mdim, nstrs)
    elmatdim = ndn*nne;             # dimension of the element matrix
    # Prepare _buffers
    dofnums = zeros(FInt, elmatdim); # degree of freedom array -- buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim); # intermediate result -- buffer
    gradXN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    gradxmN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    gradxN = fill(zero(FFlt), nne, mdim); # intermediate result -- buffer
    # Prepare _buffers
    D = fill(zero(FFlt), nstrs, nstrs); # material stiffness matrix -- buffer
    B = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    DB = fill(zero(FFlt), nstrs, elmatdim); # strain-displacement matrix -- buffer
    elmat = fill(zero(FFlt), elmatdim, elmatdim);      # element matrix -- buffer
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
    elvec = fill(zero(FFlt), elmatdim);      # element matrix -- buffer
    pu = fill(zero(FFlt), elmatdim);      # element matrix -- buffer
    # Prepare _buffers
    cauchy = fill(zero(FFlt), nstrs);      # element matrix -- buffer
    output = fill(zero(FFlt), nstrs);      # element matrix -- buffer
    sigma = fill(zero(FFlt), sdim, sdim);      # element matrix -- buffer
    # Prepare _buffers
    idx = [(1:ndn:(nne-1)*ndn+1).+(i-1) for i in 1:sdim]; # indexing vector -- buffer
    # to compute c1 = gradxmN*sigma*gradxmN';
    c1 = fill(zero(FFlt), nne, nne); # strain-displacement matrix -- buffer
    sg = fill(zero(FFlt), sdim, nne); # strain-displacement matrix -- buffer
    return _buffers = _Buffers(initialized,
    	dofnums, loc, J, csmatTJ, gradXN, gradxN, gradxmN, 
    	D, B, DB, elmat, elvec, pu,
    	X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF,
    	cauchy, output, sigma,
    	idx, c1, sg
    	)
end

"""
    AbstractFEMMDeforNonlinear <: AbstractFEMMDeforLinear

Abstract type of FEMM for nonlinear deformation.
"""
mutable struct FEMMDeforNonlinearExpl{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear} <: AbstractFEMMDeforNonlinear
	mr::Type{MR} # model reduction type
	integdomain::IntegDomain{S, F} # integration domain data
	mcsys::CSys # updater of the material orientation matrix
	material::M # material object
	statev::Vector{Vector{FFltVec}} # vector of state vectors, one for each element, and then for each integration point
	_buffers::_Buffers
end


"""
    FEMMDeforNonlinearExpl(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}

Constructor of nonlinear deformation finite element modeling machine.

The material coordinate system is the default.
"""
function FEMMDeforNonlinearExpl(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdomain);
    statev = Vector{FFltVec}[]
    for i in 1:count(integdomain.fes)
        push!(statev, [newstate(material) for j in 1:npts])
    end
    return FEMMDeforNonlinearExpl(mr, integdomain, CSys(manifdim(integdomain.fes)), material, statev, _makebuffers(false, 1, 1, 1, 1, 1))
end

"""
    FEMMDeforNonlinearExpl(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}

Constructor of nonlinear deformation finite element modeling machine.
"""
function FEMMDeforNonlinearExpl(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdomain);
    statev = Vector{FFltVec}[]
    for i in 1:count(integdomain.fes)
        push!(statev, [newstate(material) for j in 1:npts])
    end
    return FEMMDeforNonlinearExpl(mr, integdomain, mcsys, material, statev, _makebuffers(false, 1, 1, 1, 1, 1))
end

function _buff1(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.dofnums, self._buffers.loc, self._buffers.J, self._buffers.csmatTJ, self._buffers.gradXN, self._buffers.gradxN, self._buffers.gradxmN
end

function _buff2(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.D, self._buffers.B, self._buffers.DB, self._buffers.elmat
end

function _buff3(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.X, self._buffers.xn, self._buffers.xn1, self._buffers.Un, self._buffers.Un1, self._buffers.Fn, self._buffers.Fn1, self._buffers.Fnm, self._buffers.Fn1m, self._buffers.RmTF
end

function _buff4(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.B, self._buffers.elvec, self._buffers.pu
end

function _buff5(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.cauchy, self._buffers.output, self._buffers.sigma
end

function _buff6(self::FEMMDeforNonlinearExpl, geom::NodalField, u::NodalField)
    return self._buffers.idx, self._buffers.elmat, self._buffers.c1, self._buffers.sg
end

function _makebuffers!(self::FEMMDeforNonlinearExpl,  geom::NodalField{FFlt}, u::NodalField)
	if !self._buffers.initialized
		fes = self.integdomain.fes
		ndn = ndofs(u); # number of degrees of freedom per node
		nne = nodesperelem(fes); # number of nodes for element
		sdim = ndofs(geom);            # number of space dimensions
		mdim = manifdim(fes); # manifold dimension of the element
		nstrs = nstressstrain(self.mr);  # number of stresses
		self._buffers = _makebuffers(true, ndn, nne, sdim, mdim, nstrs)
	end
    return self # default is no-op
end

"""
    restoringforce(self::AbstractFEMMDeforNonlinear, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {A<:AbstractSysvecAssembler, T<:Number}

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
- `savestate` = bool flag: should we modify the material states (`savestate = true`)? Otherwise work with a copy of the material state.
"""
function restoringforce(self::FEMMDeforNonlinearExpl, assembler::A, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {A<:AbstractSysvecAssembler, T<:Number}
    fes = self.integdomain.fes
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    _makebuffers!(self, geom, un1)
    dofnums, loc, J, csmatTJ, gradXN, gradxN, gradxmN = _buff1(self, geom, un1)
    B, elvec, pu = _buff4(self, geom, un1)
    X, xn, xn1, Un, Un1, Fn, Fn1, Fnm, Fn1m, RmTF = _buff3(self, geom, un1)
    cauchy, output, sigma = _buff5(self, geom, un1)
    statev = self.statev # Initially we are assuming that the state can be overridden
    if !savestate  # unless we are instructed to work with a copy
        statev = deepcopy(self.statev) # work with copies of material state
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
            locjac!(loc, J, X, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            updatecsmat!(self.mcsys, loc, J, fes.label[i]); Rm = self.mcsys.csmat
            # At_mul_B!(csmatTJ, self.mcsys.csmat, J); # local Jacobian matrix
            gradN!(fes, gradXN, gradNparams[j], J); # wrt global c.s.
            At_mul_B!(Fn, xn, gradXN) # Previous deformation gradient
            At_mul_B!(Fn1, xn1, gradXN) # Current deformation gradient
            rotF!(Fnm, Fn, RmTF, Rm) # wrt local c.s.
            rotF!(Fn1m, Fn1, RmTF, Rm) # wrt local c.s.
            update!(self.material, statev[i][j], cauchy, output, Fn1m, Fnm, tn, dtn, loc, fes.label[i])
            gradN!(fes, gradxN, gradXN, Fn1); # wrt current c.s.
            A_mul_B!(gradxmN, gradxN, Rm) # wrt material c.s.
            Blmat!(self.mr, B, Ns[j], gradxmN, loc, Rm); # local strain-global disp
            add_btsigma!(elvec, B, - Jac * w[j] * det(Fn1), cauchy)
        end # Loop over quadrature points
        gatherdofnums!(un1, dofnums, fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, elvec, dofnums); # assemble symmetric matrix
    end # Loop over elements
    return makevector!(assembler);
end

function restoringforce(self::FEMMDeforNonlinearExpl, geom::NodalField{FFlt}, un1::NodalField{T}, un::NodalField{T}, tn::FFlt, dtn::FFlt, savestate = false) where {T<:Number}
    assembler = SysvecAssembler();
    return restoringforce(self, assembler, geom, un1, un, tn, dtn, savestate);
end

end
