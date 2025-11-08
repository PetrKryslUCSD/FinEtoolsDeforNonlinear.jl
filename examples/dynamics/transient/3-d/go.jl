# run as
# julia --project=.. -t 8 go.jl

include("cantilever_dyn_examples.jl");                                  

# Serial  execution
# cantilever_dyn_examples.neohookean_h8()   
# Parallel execution
for NTHREADS in [1 2 4 8]
	cantilever_dyn_examples.neohookean_h8_thr(NTHREADS)   
end


# exit()

# include("test/playground.jl")
# Pkg.test()

# include(".\\examples\\statics\\2-d\\tension_compression_examples.jl");
# tension_compression_examples.neohookeanad_q4()

# cd(".\\examples\\statics\\2-d")   
# include(".\\bertoldi_compression_examples.jl"); 
# bertoldi_compression_examples.neohookeanad_q4()
