using Pkg
Pkg.activate(".")
Pkg.instantiate()
using FinEtoolsDeforNonlinear

include(".\\examples\\statics\\3-d\\cantilever_examples.jl");
cantilever_examples.neohookean_h8()

