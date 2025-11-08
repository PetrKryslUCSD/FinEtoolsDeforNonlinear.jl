"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.AssemblyModule: AbstractSysvecAssembler, startassembly!, assemble!, makevector!
import SparseArrays: sparse
import LinearAlgebra: diag

"""
    SysvecAssemblerOpt

Assembler for the system vector.
Optimized for multiple assemblies of the same-size vector.
"""
mutable struct SysvecAssemblerOpt{T<:Number} <: AbstractSysvecAssembler
    F_buffer::Vector{T};
    ndofs::FInt
end

"""
    SysvecAssemblerOpt(zero::T=0.0) where {T<:Number}

Construct blank system vector assembler. The vector entries are of type `T`.
"""
function SysvecAssemblerOpt(zero::T=0.0) where {T<:Number}
    return SysvecAssemblerOpt([zero], 1)
end

"""
    startassembly!(self::SysvecAssemblerOpt{T},
      ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.
"""
function startassembly!(self::SysvecAssemblerOpt{T},  ndofs_row::FInt) where {T<:Number}
	if ndofs_row != self.ndofs
		self.ndofs = ndofs_row
		self.F_buffer = zeros(T, self.ndofs);
	end
	fill!(self.F_buffer, zero(T))
	return self
end

"""
    assemble!(self::SysvecAssemblerOpt{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssemblerOpt{T}, vec::MV, dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}
	for i in 1:length(dofnums)
		gi = dofnums[i];
		if (0 < gi <= self.ndofs)
			self.F_buffer[gi] = self.F_buffer[gi] + vec[i];
		end
	end
	return self
end

"""
    makevector!(self::SysvecAssemblerOpt)

Make the global vector.
Beware: the buffer itself is returned as the vector, *not* a copy.
"""
function makevector!(self::SysvecAssemblerOpt)
	return self.F_buffer;
end



end
