module MatDeforNonlinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor
import LinearAlgebra:  det

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

"""
    totlag2curr!(c, C, F)

Convert a total Lagrangean constitutive matrix to a current Lagrangean one
(sometimes known as "Eulerian").

- `C`    = Lagrangean constitutive matrix, 6x6, symmetric
- `F`    = current deformation gradient, F_iJ = partial x_i / partial X_J

The transformation is c_ijkl = 1/J C_IJKL F_iI F_jJ F_kK F_lL. In the present
case the fourth-order tensor is represented with a 6 x 6 matrix.
"""
function totlag2curr!(c, C, F) 
	@assert size(F) == (3, 3)
	@assert size(C) == (6, 6)
    t4 = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3));
    t2 = ((1, 4, 5), (4, 2, 6), (5, 6, 3)) 
    for mj in 1:size(c, 2)
    	k, l = t4[mj]
    	for mi in 1:size(c, 1)
    		i, j = t4[mi]
    		s = 0.0
    		@inbounds for I in 1:3
    			@inbounds for J in 1:3
    				mI = t2[I][J]
    				FiIFjJ = F[i, I] * F[j, J]
    				@inbounds for K in 1:3
    					@inbounds for L in 1:3
    						mJ = t2[K][L]
    						s += FiIFjJ * F[k, K] * F[l, L] * C[mI, mJ]
    					end
    				end
    			end
    		end
    		c[mi, mj] = s
    	end
    end
    c .*= (1.0 / det(F))
    return c
end

"""
    totlag2currsymm!(c, C, F)

Convert a total Lagrangean constitutive matrix to a current Lagrangean one
(sometimes known as "Eulerian").

- `C`    = Lagrangean constitutive matrix, 6x6, symmetric
- `F`    = current deformation gradient, F_iJ = partial x_i / partial X_J

The transformation is c_ijkl = 1/J C_IJKL F_iI F_jJ F_kK F_lL. In the present
case the fourth-order tensor is represented with a 6 x 6 matrix.

!!! note
The Lagrangean material stiffness matrices, both input and output, are
presumed symmetric.
"""
function totlag2currsymm!(c, C, F) 
	@assert size(F) == (3, 3)
	t4 = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3));
	t2 = ((1, 4, 5), (4, 2, 6), (5, 6, 3)) 
    for mj in 1:size(c, 2)
    	k, l = t4[mj]
    	for mi in 1:mj
    		i, j = t4[mi]
    		s = 0.0
    		@inbounds for I in 1:3
    			@inbounds for J in 1:3
    				mI = t2[I][J]
    				fiIjJ = F[i, I] * F[j, J]
    				@inbounds for K in 1:3
    					@inbounds for L in 1:3
    						mJ = t2[K][L]
    						s += fiIjJ * F[k, K] * F[l, L] * C[mI, mJ]
    					end
    				end
    			end
    		end
    		c[mj, mi] = c[mi, mj] = s
    	end
    end
    c .*= (1.0 / det(F))
    return c
end

"""
    totlag2curr4th!(c, C, F)

Convert a total Lagrangean constitutive matrix to a current Lagrangean one
(sometimes known as "Eulerian").

- `C`    = Lagrangean constitutive matrix, fourth-order tensor
- `F`    = current deformation gradient, F_iJ = partial x_i / partial X_J

The transformation is c_ijkl = 1/J C_IJKL F_iI F_jJ F_kK F_lL. Both the input
and the output are fourth-order tensors.
"""
function totlag2curr4th!(c, C, F)
	@assert size(F) == (3, 3)
	@assert size(c) == (3, 3, 3, 3)
	@assert size(C) == (3, 3, 3, 3)
    n = size(F, 1)
    c .= 0.0
    @inbounds for i in 1:n
    	@inbounds for j in 1:n
    		@inbounds for k in 1:n
    			@inbounds for l in 1:n
    				@inbounds for I in 1:n
    					@inbounds for J in 1:n
    						@inbounds for K in 1:n
    							@inbounds for L in 1:n
    								c[i, j, k, l] += F[i, I] * F[j, J] * F[k, K] * F[l, L] * C[I, J, K, L]
    							end
    						end
    					end
    				end
    			end
    		end
    	end
    end
    c ./= det(F)
    return c
end

end
