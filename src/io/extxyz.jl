
export write_extxyz
export read_extxyz

import ExtXYZ

write_extxyz(file, atoms, R::Matrix, cell) =
    ExtXYZ.write_frame(file, to_extxyz_dict(atoms, R, cell))

write_extxyz(file, atoms, R::Vector{<:Matrix}, cell) =
    ExtXYZ.write_frames(file, to_extxyz_dict.(Ref(atoms), R, Ref(cell)))

function read_extxyz(file)

    dict = ExtXYZ.read_frames(file, 1)[1]
    atoms, R, cell = from_extxyz_dict(dict)
    positions = [R]

    for dict in ExtXYZ.read_frames(file, Iterators.countfrom(2))
        _, R, _ = from_extxyz_dict(dict)
        push!(positions, R)
    end

    return atoms, positions, cell
end

function to_extxyz_dict(atoms::Atoms, R::Matrix, cell::PeriodicCell)

    dict = Dict{String, Any}()
    dict["N_atoms"] = length(atoms)

    dict["pbc"] = cell.periodicity

    dict["cell"] = au_to_ang.(permutedims(cell.vectors, (2,1)))

    dict["arrays"] = Dict{String,Any}()
    dict["arrays"]["species"] = String.(atoms.types) |> Array
    dict["arrays"]["pos"] = au_to_ang.(R)

    dict["info"] = Dict{String,Any}()

    return dict
end

function from_extxyz_dict(dict::Dict{String,Any})

    atoms = Atoms(Symbol.(dict["arrays"]["species"]))
    R = ang_to_au.(dict["arrays"]["pos"])
    cell = PeriodicCell{eltype(dict["cell"])}(ang_to_au.(dict["cell"]), dict["pbc"])

    return atoms, R, cell
end
