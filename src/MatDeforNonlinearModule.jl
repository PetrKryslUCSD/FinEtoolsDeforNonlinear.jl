module MatDeforNonlinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor

"""
    AbstractMatDeforNonlinear <: AbstractMatDefor

Abstract nonlinear material.

"""
abstract type AbstractMatDeforNonlinear <: AbstractMatDefor; end

"""
    tangentmoduli!(self::M, D::FFltMat, statev::FFltVec, Fn1::FFltMat,
        Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
        where {M<:AbstractMatDeforNonlinear}

Calculate the material stiffness matrix.

# Arguments
- `self` = material
- `D` = matrix of tangent moduli, supplied as a buffer and overwritten.
  Returned as output.
- `statev` = material state vector, the content of this vector must not change
    inside this function.
- `Fn1` = deformation gradient at time `tn1 = tn + dtn`,
- `Fn` = deformation gradient at time `tn`,
- `tn` = time in step `n`
- `dtn` = increment of time
- `loc` = location of the integration point in the reference coordinates (time `t0`),
- `label` = label of the element containing the integration point

!!! note
The deformation gradients and the matrix of the tangent moduli are expressed
with respect to the local material coordinate system.
"""
function tangentmoduli!(self::M, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt) where {M<:AbstractMatDeforNonlinear}
    return self.tangentmoduli!(self, D, statev, Fn1, Fn, tn, dtn, loc, label)
end

"""
    update!(self::M, statev::FFltVec, stress::FFltVec, output::FFltVec,
        Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt,
        loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
        where {M<:AbstractMatDeforNonlinear}

Update material state.

# Arguments
- `self` = material
- `statev` = state variables: array which is (if necessary) allocated  in an appropriate
     size, filled with the state variables, and returned. The contents of this
     vector may change as the state of the material is updated by the logic inside
     this function. If this change is to be saved, it must happen outside of this
     function.
- `cauchy` = Cauchy stress vector, allocated by the caller with a size of the
  number of stress and strain components, `nstressstrain`. The components of
  the stress vector are calculated and stored in the `stress` vector.
- `output` =  array which is (if necessary) allocated  in an appropriate size,
     filled with the output quantity, and returned.
- `Fn1` = deformation gradient at time `tn1 = tn + dtn`,
- `Fn` = deformation gradient at time `tn`,
- `tn` = time in step `n`
- `dtn` = increment of time
- `loc` = location of the integration point in the reference coordinates (time `t0`),
- `label` = label of the element containing the integration point

# Output
- `cauchy` = Cauchy stress vector
- `output` = output array

!!! note
The deformation gradients and the stress vector are expressed
with respect to the local material coordinate system.
"""
function update!(self::M, statev::FFltVec, cauchy::FFltVec, output::FFltVec,  Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing) where {M<:AbstractMatDeforNonlinear}
    return self.update!(self, statev, cauchy, output, Fn1, Fn, tn, dtn, loc, label, quantity)
end

"""
    newstate(self::M) where {M<:AbstractMatDeforNonlinear}

Create an initial material state at an integration point.
"""
function newstate(self::M) where {M<:AbstractMatDeforNonlinear}
    return FFlt[] # the default material state is an empty array: no state to store
end

end
