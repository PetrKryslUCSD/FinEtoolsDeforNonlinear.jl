using Pkg
Pkg.activate(".")  
Pkg.instantiate()
using FinEtoolsDeforNonlinear
using Profile
using LinearAlgebra

include("./examples/dynamics/transient/3-d/weird_timing_examples.jl");                                  

# Serial  execution
# weird_timing_examples.neohookean_h8()   
# Parallel execution
LinearAlgebra.BLAS.set_num_threads(8) 
@show ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
for NTHREADS in [1 2 4 8]
	weird_timing_examples.neohookean_h8_thr(NTHREADS)   
	# weird_timing_examples.neohookean_h8_thr_2(NTHREADS)  
end

# (10.7 / 5.25) / 2
# (10.7 / 2.75) / 4
# (10.7 / 1.48) / 8

# (10.09 / 5.30) / 2
# (10.09 / 2.89) / 4
# (10.09 / 2.13) / 8
# (10.09 / 1.50) / 16

(23.6784 / 12.280) / 2
# exit()

# include("test/playground.jl")
# Pkg.test()

# include("./examples/statics/2-d/tension_compression_examples.jl");
# tension_compression_examples.neohookeanad_q4()

# cd("./examples/statics/2-d")   
# include("./bertoldi_compression_examples.jl"); 
# bertoldi_compression_examples.neohookeanad_q4()
