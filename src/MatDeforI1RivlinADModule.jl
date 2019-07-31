module MatDeforI1RivlinADModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtoolsDeforLinear.DeforModelRedModule: AbstractDeforModelRed, DeforModelRed3D, DeforModelRed2DStrain, DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed1D, nstressstrain, nthermstrain
import FinEtoolsDeforLinear.MatDeforModule: AbstractMatDefor, stress6vto3x3t!, stress3vto2x2t!, stress4vto3x3t!, stress4vto3x3t!, stress3x3tto6v!, strain3x3tto6v!, strain6vto3x3t!, strain6vdet, strain6vtr
import ..MatDeforNonlinearModule: AbstractMatDeforNonlinear, totalLagrangean2current!
using LinearAlgebra: Transpose, Diagonal, mul!
At_mul_B!(C, A, B) = mul!(C, Transpose(A), B)
A_mul_B!(C, A, B) = mul!(C, A, B)
using LinearAlgebra: eigen, eigvals, norm, cholesky, cross, dot, log, diagm, det
using ForwardDiff: gradient, hessian

"""
	MatDeforI1RivlinAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear

Type for triaxial two-term I1-based Rivlin hyperelastic material. Described in
K Bertoldi and others, Journal of Mechanics and Physics of Solids 56 (2008)
2642-2668.

Material parameters suggested for one specific type of rubber in that paper
was specified by parameters
```
c1 = 0.55 MPa
c2 = 0.3 MPa
K = 50 * mu = 55 MPa
```
Initial Young's modulus 3.25 MPa, which corresponds to nu = 0.4772727
"""
struct  MatDeforI1RivlinAD{MR<:AbstractDeforModelRed, MTAN<:Function, MUPD<:Function} <: AbstractMatDeforNonlinear
	mr::Type{MR} # model reduction type
	mass_density::FFlt # mass density
	c1::FFlt # material constant: c1 = mu/2 (half the shear modulus).
	c2::FFlt # material constant: neo-hookean results if c2 = 0.
	K::FFlt # bulk modulus
	tangentmoduli!::MTAN # Function to return the tangent moduli matrix
	update!::MUPD # Function to update the material state
end

################################################################################
# 3-D solid model
################################################################################

"""
    MatDeforI1RivlinAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, E::FFlt, nu::FFlt)

Create triaxial two-term I1-based Rivlin hyperelastic material.
"""
function MatDeforI1RivlinAD(mr::Type{DeforModelRed3D}, mass_density::FFlt, c1::FFlt, c2::FFlt, K::FFlt)
	_mI = diagm(0=>[1, 1, 1, 0.5, 0.5, 0.5]);
	_m1 = vec(FFlt[1 1 1 0 0 0]) ;
	_m1m1 = _m1*_m1';
	_I3 = [1.0 0 0; 0 1.0 0; 0 0 1.0]
	# neo-hookean:  mu/2*(trC-3) - mu*lJ + lambda/2*(lJ)^2;
	# two-term I1-based Rivlin:
	function strainenergy(Cv, MR, c1, c2, K)
		trC = strain6vtr(MR, Cv)
		J = sqrt(strain6vdet(MR, Cv))
		lJ = log(J)
		return c1*(trC-3) + c2*(trC-3)^2 - 2*c1*(lJ) + K/2*(J - 1)^2;
	end
	function tangentmoduli3d!(self::MatDeforI1RivlinAD, D::FFltMat, statev::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat, label::FInt)
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strain3x3tto6v!(Cv, C)
		Dtotal = 4 .* hessian(Cv -> strainenergy(Cv, mr, self.c1, self.c2, self.K), Cv);
		return totalLagrangean2current!(D, Dtotal, Fn1)
	end
	function update3d!(self::MatDeforI1RivlinAD, statev::FFltVec, cauchy::FFltVec, output::FFltVec, Fn1::FFltMat, Fn::FFltMat, tn::FFlt, dtn::FFlt, loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
		@assert length(cauchy) == nstressstrain(self.mr)
		C = Fn1'*Fn1;
		Cv = fill(0.0, 6)
		strain3x3tto6v!(Cv, C)
		Sv = 2 * gradient(Cv -> strainenergy(Cv, mr, self.c1, self.c2, self.K), Cv);
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
	return MatDeforI1RivlinAD(mr, mass_density, c1, c2, K,
		tangentmoduli3d!, update3d!)
end

"""
    MatDeforI1RivlinAD(mr::Type{MR}, c1::FFlt, c2::FFlt, K::FFlt) where {MR}

Create two-term I1-based Rivlin hyperelastic material.

The mass density is the default value.
"""
function MatDeforI1RivlinAD(mr::Type{MR}, c1::FFlt, c2::FFlt, K::FFlt) where {MR}
	mass_density = 1.0
	return MatDeforI1RivlinAD(mr, mass_density, c1, c2, K)
end

"""
    MatDeforI1RivlinAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}

Create two-term I1-based Rivlin hyperelastic material equivalent to a simple neo-hookean material.

The mass density is the default value.
"""
function MatDeforI1RivlinAD(mr::Type{MR}, E::FFlt, nu::FFlt) where {MR}
	mass_density = 1.0
	c1 = E / 2 / (1 + nu) / 2 # half of the shear modulus
	c2 = 0.0
	K = E / 3 / (1 - 2*nu) # bulk modulus
	return MatDeforI1RivlinAD(mr, mass_density, c1, c2, K)
end


end
