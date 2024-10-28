
using Test
using NQCBase


rtol = 1e-7 # Relative tolerance for approximate equals - This was set for Julia 1.10 with ExtXYZ v0.2

atoms = Atoms([:H, :C, :O, :N])
cell = PeriodicCell(rand(3, 3) .* 10)
R = rand(3, 4) .* 10

@testset "to/from_extxyz_dict" begin
    dict = NQCBase.to_extxyz_dict(atoms, R, cell)
    @test dict["cell"] ≈ au_to_ang.(permutedims(cell.vectors, (2,1)))
    natoms, nR, ncell = NQCBase.from_extxyz_dict(dict)
    @test ncell.vectors ≈ cell.vectors
    @test ncell.inverse ≈ cell.inverse
end

@testset "Single frame" begin
    write_extxyz("output.xyz", atoms, R, cell)
    new_atoms, new_R, new_cell = read_extxyz("output.xyz")
    @test new_atoms == atoms
    @test isapprox(new_cell.vectors, cell.vectors; rtol=rtol)
    @test isapprox(new_cell.inverse, cell.inverse; rtol=rtol)
    @test isapprox(new_cell.periodicity, cell.periodicity; rtol=rtol)
    @test isapprox(new_R, [R]; rtol=rtol)
end

@testset "Multiple frames" begin
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell(rand(3, 3) .* 10)
    R = [rand(3, 4) .* 10 for _=1:100]
    write_extxyz("output.xyz", atoms, R, cell)
    new_atoms, new_R, new_cell = read_extxyz("output.xyz")
    @test new_atoms == atoms
    @test isapprox(new_cell.vectors, cell.vectors; rtol=rtol)
    @test isapprox(new_cell.inverse, cell.inverse; rtol=rtol)
    @test isapprox(new_cell.periodicity, cell.periodicity; rtol=rtol)
    @test isapprox(new_R, R; rtol=rtol)
end

rm("output.xyz")
