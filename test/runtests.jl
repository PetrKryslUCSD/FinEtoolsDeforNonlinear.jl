

using Test
@time @testset "Materials" begin 
include("test_materials.jl") 
include("test_materials2.jl") 
end

@time @testset "Operations" begin include("test_operations.jl") end
true
