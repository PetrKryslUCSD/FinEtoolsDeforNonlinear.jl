module MatDeforNeohookeanModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!, stressttov!
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear
import LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
import LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot, log, diagm, det


"""
	MatDeforNeohookean{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear

Type for triaxial neohookean hyperelastic material.

"""
struct  MatDeforNeohookean{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear
	mr::Type{MR} # model reduction type
	mass_density::FFlt # mass density
	E::FFlt # Young's modulus
	nu::FFlt # Poisson ratio
	_lambda::FFlt
	_mu::FFlt
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
end

# function _threedD(E::FFlt, nu::FFlt)
# 	mI = Matrix(Diagonal([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))
# 	m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
# 	_lambda = E * nu / (1 + nu) / (1 - 2*(nu));
# 	_mu = E / 2. / (1+nu);
# 	D = _lambda * m1 * m1' + 2. * _mu * mI;
# 	return D
# end


################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforNeohookean(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)

Create triaxial neohookean hyperelastic material.
"""
function MatDeforNeohookean(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)
	_mI = diagm(0=>[1, 1, 1, 0.5, 0.5, 0.5]);
	_m1 = vec(FFlt[1 1 1 0 0 0]) ;
	_m1m1 = _m1*_m1';
	_I3 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
	function tangentmoduli3d!(self::MatDeforNeohookean, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		J = det(Fn1);
		copyto!(D, (self._lambda / J) * _m1m1 + 2 * (self._mu - self._lambda*log(J))/J * _mI);
		return D
	end
	function update3d!(self::MatDeforNeohookean, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		# Finger deformation tensor
		b = Fn1*Fn1';
		J = det(Fn1); # Jacobian of the deformation gradient
		sigma = (self._mu/J) .* (b - _I3) + (self._lambda *log(J)/J) .* _I3;
		stressttov!(mr, cauchy, sigma)
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
	return MatDeforNeohookean(mr, mass_density, E, nu, _lambda, _mu,
		tangentmoduli3d!, update3d!)
end

"""
    MatDeforNeohookean(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create neohookean isotropic elastic material.

The mass density is the default value.
"""
function MatDeforNeohookean(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
	mass_density = 1.0
	return MatDeforNeohookean(mr, mass_density, E, nu)
end



end
