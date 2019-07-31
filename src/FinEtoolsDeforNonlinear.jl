"""
FinEtoolsDeforNonlinear (C) 2017-2019, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for nonlinear stress analysis problems.
"""
module FinEtoolsDeforNonlinear

__precompile__(true)

include("MatDeforNonlinearModule.jl")
include("MatDeforNeohookeanModule.jl")
include("MatDeforNeohookeanADModule.jl")
include("MatDeforStVKModule.jl")
include("MatDeforStVKADModule.jl")
include("MatDeforI1RivlinADModule.jl")

include("FEMMDeforNonlinearBaseModule.jl")
include("FEMMDeforNonlinearModule.jl")

include("AlgoDeforNonlinearModule.jl")

# Exports follow:


end # module
