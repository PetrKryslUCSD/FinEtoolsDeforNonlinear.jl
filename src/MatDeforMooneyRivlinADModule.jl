module MatDeforMooneyRivlinADModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!, stressttov!, strainttov!, strainvdet, strainvtr
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, totlag2currsymm!
using LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot, log, diagm, det
using ForwardDiff: gradient, hessian

"""
	MatDeforMooneyRivlinAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear

Type for triaxial Mooney-Rivlin hyperelastic material.

"""
struct  MatDeforMooneyRivlinAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear
	mr::Type{MR} # model reduction type
	mass_density::FFlt # mass density
	E::FFlt # Young's modulus
	nu::FFlt # Poisson ratio
	_lambda::FFlt
	_mu::FFlt
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
end

################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforMooneyRivlinAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)

Create triaxial Mooney-Rivlin hyperelastic material.
"""
function MatDeforMooneyRivlinAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)
	_mI = diagm(0=>[1, 1, 1, 0.5, 0.5, 0.5]);
	_m1 = vec(FFlt[1 1 1 0 0 0]) ;
	_m1m1 = _m1*_m1';
	_I3 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
	function strainenergy(Cv, MR, lambda, mu)
		trC = strain6vtr(MR, Cv)
		J = sqrt(strain6vdet(MR, Cv))
		lJ = log(J)
		return mu/2*(trC-3) - mu*lJ + lambda/2*(lJ)^2;
	end
	function tangentmoduli3d!(self::MatDeforMooneyRivlinAD, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strain3x3tto6v!(Cv, C)
		Dtotal = 4 .* hessian(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		return totlag2currsymm!(D, Dtotal, Fn1)
	end
	function update3d!(self::MatDeforMooneyRivlinAD, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		# Note: the code below allocates. This may not be a big deal when the
		# code gets called only a few times, but in explicit dynamics
		# situations  this would be very expensive. The material type needs to
		# provide these temporary arrays as buffers. 
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strain3x3tto6v!(Cv, C)
		Sv = 2 * gradient(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		S = fill(0.0, 3, 3)
		stress6vto3x3t!(S, Sv)
		cauchyt = Fn1*(S/det(Fn1))*Fn1'; # Cauchy stress
		stress3x3tto6v!(cauchy, cauchyt)
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 6) || (output = zeros(6)) # make sure we can store it
			copyto!(output, cauchy);
		elseif quantity == :pressure || quantity == :Pressure
			output[1]  =  -sum(cauchy[1:3])/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			ep = eigen(sigma);
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
	_lambda = E * nu / (1 + nu) / (1 - 2*(nu));
	_mu = E / 2. / (1+nu);
	return MatDeforMooneyRivlinAD(mr, mass_density, E, nu, _lambda, _mu,
		tangentmoduli3d!, update3d!)
end

"""
    MatDeforMooneyRivlinAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create Mooney-Rivlin isotropic elastic material.

The mass density is the default value.
"""
function MatDeforMooneyRivlinAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
	mass_density = 1.0
	return MatDeforMooneyRivlinAD(mr, mass_density, E, nu)
end



end
