using NonadiabaticDynamicsBase
using Test, SafeTestsets

@time @safetestset "Atoms tests" begin include("atoms.jl") end
@time @safetestset "Cells tests" begin include("cells.jl") end
