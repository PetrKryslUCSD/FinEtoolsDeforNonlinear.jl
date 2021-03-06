module MatDeforStVKModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!, stressttov!, strainttov!
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, totlag2currsymm!
import LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot, log, diagm, det


"""
	MatDeforStVK{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear

Type for triaxial St Venant-Kirchhoff hyperelastic material.

"""
struct  MatDeforStVK{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear
	mr::Type{MR} # model reduction type
	mass_density::FFlt # mass density
	E1::FFlt 
	E2::FFlt 
	E3::FFlt 
	G12::FFlt 
	G13::FFlt 
	G23::FFlt 
	nu12::FFlt 
	nu13::FFlt 
	nu23::FFlt 
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
	_D::FFltMat # cache the constant stiffness matrix of the material
end

################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforStVK(mr::Type{DeforModelRed3D}, mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt)

Create triaxial St Venant-Kirchhoff hyperelastic material.

In general, the material is assumed to be orthotropic. There is a
specialized constructor for an isotropic version.
"""
function MatDeforStVK(mr::Type{DeforModelRed3D}, mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt)
	compliance =[1/E1      -nu12/E1    -nu13/E1  0   0   0;
				-nu12/E1     1/E2      -nu23/E2  0   0   0;
				-nu13/E1   -nu23/E2       1/E3   0   0   0;
				0           0           0      1/G12 0   0;
				0           0           0        0 1/G13 0;
				0           0           0        0   0 1/G23];
	_D = inv(compliance);
	_I3 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
	function tangentmoduli3d!(self::MatDeforStVK, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		return totlag2currsymm!(D, self._D, Fn1)
	end
	function update3d!(self::MatDeforStVK, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		# Note: the code below allocates. This may not be a big deal when the
		# code gets called only a few times, but in explicit dynamics
		# situations  this would be very expensive. The material type needs to
		# provide these temporary arrays as buffers. 
		C = Fn1'*Fn1;# Cauchy-Green deformation tensor
		E = 1/2*(C-_I3)# Green-Lagrange strain
		Egl = fill(0.0, nstressstrain(self.mr))
		strainttov!(mr, Egl, E);
		S = self._D * Egl;
		St = fill(0.0, 3, 3)
		stressvtot!(mr, St, S)
		J = det(Fn1);
		cauchyt = Fn1*(St/J)*Fn1'; # Cauchy stress
		stressttov!(mr, cauchy, cauchyt)
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 6) || (output = zeros(6)) # make sure we can store it
			copyto!(output, cauchy);
		elseif quantity == :pressure || quantity == :Pressure
			output[1]  =  -sum(cauchy[1:3])/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			ep = eigen(cauchyt);
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			s1=cauchy[1]; s2=cauchy[2]; s3=cauchy[3];
			s4=cauchy[4]; s5=cauchy[5]; s6=cauchy[6];
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	return MatDeforStVK(mr, mass_density, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt,
		tangentmoduli3d!, update3d!, _D)
end

"""
    MatDeforStVK(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create St Venant-Kirchhoff isotropic elastic material.

The mass density is the default value.
"""
function MatDeforStVK(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
	mass_density = 1.0
	G = E / (2 * (1 + nu))
	return MatDeforStVK(mr, mass_density, E, E, E, G, G, G, nu, nu, nu)
end

"""
    MatDeforStVK(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create St Venant-Kirchhoff isotropic elastic material.

The mass density is the default value.
"""
function MatDeforStVK(mr::Type{MR}, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt) where {MR}
	mass_density = 1.0
	return MatDeforStVK(mr, mass_density, E1::FFlt, E2::FFlt, E3::FFlt, G12::FFlt, G13::FFlt, G23::FFlt, nu12::FFlt, nu13::FFlt, nu23::FFlt)
end

end
