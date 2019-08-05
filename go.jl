using Pkg
Pkg.activate(".")
# Pkg.instantiate()
using FinEtoolsDeforNonlinear

# include(".\\examples\\statics\\2-d\\tension_compression_examples.jl");
# tension_compression_examples.neohookeanad_q4()

Pkg.test()