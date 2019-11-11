using Test

@time @testset "Materials" begin 
include("test_materials.jl") 
include("test_materials2.jl") 
include("test_materials3.jl") 
include("test_materials4.jl") 
include("test_materials5.jl") 
end

@time @testset "Operations" begin 
include("test_operations.jl") 
end

true
