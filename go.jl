using Pkg
Pkg.activate(".")  
Pkg.instantiate()
using FinEtoolsDeforNonlinear
using Profile

include("C:\\Users\\PK\\Documents\\work\\FinEtoolsDeforNonlinear.jl\\examples\\dynamics\\transient\\3-d\\cantilever_dyn_examples.jl");                                  

@time cantilever_dyn_examples.neohookean_h8()   
Profile.clear_malloc_data() 
@time cantilever_dyn_examples.neohookean_h8()   

exit()

# include("test/playground.jl")
# Pkg.test()

# include(".\\examples\\statics\\2-d\\tension_compression_examples.jl");
# tension_compression_examples.neohookeanad_q4()

# cd(".\\examples\\statics\\2-d")   
# include(".\\bertoldi_compression_examples.jl"); 
# bertoldi_compression_examples.neohookeanad_q4()
