using Pkg
Pkg.activate(".")  
Pkg.instantiate()
using FinEtoolsDeforNonlinear
using Profile

include(".\\examples\\dynamics\\transient\\3-d\\cantilever_dyn_examples.jl");                                  

cantilever_dyn_examples.neohookean_h8()   
# Profile.clear_malloc_data() 
cantilever_dyn_examples.neohookean_h8_thr()   

# exit()

# include("test/playground.jl")
# Pkg.test()

# include(".\\examples\\statics\\2-d\\tension_compression_examples.jl");
# tension_compression_examples.neohookeanad_q4()

# cd(".\\examples\\statics\\2-d")   
# include(".\\bertoldi_compression_examples.jl"); 
# bertoldi_compression_examples.neohookeanad_q4()
