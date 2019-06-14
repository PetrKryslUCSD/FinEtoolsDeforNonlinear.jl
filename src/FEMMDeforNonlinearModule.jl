"""
    FEMMDeforNonlinearModule

Module for nonlinear deformation models. Standard isoparametric finite elements.
"""
module FEMMDeforNonlinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!, gathervalues_asmat!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import ..FEMMDeforNonlinearBaseModule: AbstractFEMMDeforNonlinear
import FinEtools.CSysModule: CSys, updatecsmat!
import FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed2DAxisymm, nstressstrain, nthermstrain, Blmat!, divmat, vgradmat
import FinEtools.MatrixUtilityModule: add_btdb_ut_only!, complete_lt!, add_btv!, locjac!, add_nnt_ut_only!
import FinEtools.MatDeforModule: rotstressvec!
import FinEtools.MatModule: massdensity
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, tangentmoduli!, update!, newstate
import FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
import LinearAlgebra: Transpose, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: norm, dot

"""
    FEMMDeforNonlinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear} <: AbstractFEMMDeforNonlinear

Class for nonlinear deformation finite element modeling machine.
"""
mutable struct FEMMDeforNonlinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear} <: AbstractFEMMDeforNonlinear
    mr::Type{MR} # model reduction type
    integdomain::IntegDomain{S, F} # integration domain data
    mcsys::CSys # updater of the material orientation matrix
    material::M # material object
    statev::Vector{Vector{FFltVec}} # vector of state vectors, one for each element, and then for each integration point
end

"""
    FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}

Constructor of nonlinear deformation finite element modeling machine.

The material coordinate system is the default.
"""
function FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdomain);
    statev = Vector{FFltVec}[]
    for i in 1:count(integdomain.fes)
        push!(statev, [newstate(material) for j in 1:npts])
    end
    return FEMMDeforNonlinear(mr, integdomain, CSys(manifdim(integdomain.fes)), material, statev)
end

"""
    FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}

Constructor of nonlinear deformation finite element modeling machine.
"""
function FEMMDeforNonlinear(mr::Type{MR}, integdomain::IntegDomain{S, F}, mcsys::CSys, material::M) where {MR<:AbstractDeforModelRed, S<:AbstractFESet, F<:Function, M<:AbstractMatDeforNonlinear}
    @assert mr == material.mr "Model reduction is mismatched"
    @assert (integdomain.axisymmetric) || (mr != DeforModelRed2DAxisymm) "Axially symmetric requires axisymmetric to be true"
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(integdomain);
    statev = Vector{FFltVec}[]
    for i in 1:count(integdomain.fes)
        push!(statev, [newstate(material) for j in 1:npts])
    end
    return FEMMDeforNonlinear(mr, integdomain, mcsys, material, statev)
end
end
