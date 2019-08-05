module MatDeforNeohookeanADModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stressvtot!, stressttov!, strainttov!, strainvdet, strainvtr
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, totalLagrangean2current!
using LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot, log, diagm, det
using ForwardDiff: gradient, hessian

"""
	MatDeforNeohookeanAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear

Type for triaxial neohookean hyperelastic material.

"""
struct  MatDeforNeohookeanAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear
	mr::Type{MR} # model reduction type
	mass_density::FFlt # mass density
	E::FFlt # Young's modulus
	nu::FFlt # Poisson ratio
	_lambda::FFlt
	_mu::FFlt
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
end

"""
    MatDeforNeohookeanAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create neohookean isotropic elastic material.

The mass density is the default value.
"""
function MatDeforNeohookeanAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
	mass_density = 1.0
	return MatDeforNeohookeanAD(mr, mass_density, E, nu)
end

################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforNeohookeanAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)

Create triaxial neohookean hyperelastic material.
"""
function MatDeforNeohookeanAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)
	function strainenergy(Cv, mr, lambda, mu)
		trC = strainvtr(mr, Cv)
		J = sqrt(strainvdet(mr, Cv))
		lJ = log(J)
		return mu/2*(trC-3) - mu*lJ + lambda/2*(lJ)^2;
	end
	function tangentmoduli3d!(self::MatDeforNeohookeanAD, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strainttov!(mr, Cv, C)
		Dtotal = 4 .* hessian(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		return totalLagrangean2current!(D, Dtotal, Fn1)
	end
	function update3d!(self::MatDeforNeohookeanAD, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strainttov!(mr, Cv, C)
		Sv = 2 * gradient(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		S = fill(0.0, 3, 3)
		stressvtot!(mr, S, Sv)
		cauchyt = Fn1*(S/det(Fn1))*Fn1'; # Cauchy stress
		stressttov!(mr, cauchy, cauchyt)
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
	return MatDeforNeohookeanAD(mr, mass_density, E, nu, _lambda, _mu,
		tangentmoduli3d!, update3d!)
end


################################################################################
# 2-D plane strain model
################################################################################

"""
    MatDeforNeohookeanAD(mr::Type{DeforModelRed2DStrain}, mass_density::FFlt, E::FFlt, nu::FFlt)

Create triaxial neohookean hyperelastic material.
"""
function MatDeforNeohookeanAD(mr::Type{DeforModelRed2DStrain}, mass_density::FFlt, E::FFlt, nu::FFlt)
	function strainenergy(Cv, mr, lambda, mu)
		trC = strainvtr(mr, Cv)
		J = sqrt(strainvdet(mr, Cv))
		lJ = log(J)
		return mu/2*(trC-3) - mu*lJ + lambda/2*(lJ)^2;
	end
	function tangentmoduli2dstrain!(self::MatDeforNeohookeanAD, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		Fn13d = fill(0.0, 3, 3)
		Fn13d[1:2, 1:2] = Fn1
		Fn13d[3, 3] = 1.0
		C = Fn13d'*Fn13d;
		Cv = fill(0.0, 6)
		strainttov!(mr, Cv, C)
		Dtotal = 4 .* hessian(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		Dcurrent = fill(0.0, size(Dtotal))
		totalLagrangean2current!(Dcurrent, Dtotal, Fn13d)
		fill!(D, 0.0)
		D[1:2, 1:2] = Dcurrent[1:2, 1:2]
		D[3, 3] = Dcurrent[4, 4]
		return D
	end
	function update2dstrain!(self::MatDeforNeohookeanAD, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		Fn13d = fill(0.0, 3, 3)
		Fn13d[1:2, 1:2] = Fn1
		Fn13d[3, 3] = 1.0
		C = Fn13d'*Fn13d;
		Cv = fill(0.0, 6)
		strainttov!(mr, Cv, C)
		Sv = 2 * gradient(Cv -> strainenergy(Cv, mr, self._lambda, self._mu), Cv);
		S = fill(0.0, 3, 3)
		stressvtot!(mr, S, Sv)
		cauchyt = Fn13d*(S/det(Fn1))*Fn13d'; # Cauchy stress
		stressttov!(mr, cauchy, cauchyt)
		if quantity == :nothing
			#Nothing to be copied to the output array
		elseif quantity == :cauchy || quantity == :Cauchy
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output, cauchy);
		elseif quantity == :pressure || quantity == :Pressure
			output[1]  =  -sum(tr(cauchyt))/3.
		elseif quantity == :princCauchy || quantity == :princcauchy
			ep = eigen(cauchyt);
			(length(output) >= 3) || (output = zeros(3)) # make sure we can store it
			copyto!(output,  sort(ep.values, rev=true));
		elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
			s1=cauchyt[1, 1]; s2=cauchyt[2, 2]; s3=cauchyt[3, 3];
			s4=cauchyt[1, 2]; s5=cauchyt[1, 3]; s6=cauchyt[2, 3];
			(length(output) >= 1) || (output = zeros(1)) # make sure we can store it
			output[1] = sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))
		end
		return output
	end
	_lambda = E * nu / (1 + nu) / (1 - 2*(nu));
	_mu = E / 2. / (1+nu);
	return MatDeforNeohookeanAD(mr, mass_density, E, nu, _lambda, _mu,
		tangentmoduli2dstrain!, update2dstrain!)
end


end
