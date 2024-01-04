module NQCDynamicsBindings

using NQCBase,NQCDynamics

"""
    distance(config::Any, i1, i2)
    
Interatomic distance in Angstrom for DynamicsVariables. 
"""
function NQCBase.Structure.distance(config::Any, i1::int_or_index, i2::int_or_index)
    return NQCBase.Structure.distance(pos,i1,i2)
end

"""
    NQCBase.Structure.fractional_mass(sim::NQCDynamics.AbstractSimulation, index1::Int, index2::Int)

Returns m1/(m1+m2) and m2/(m1+m2) as a vector. 
"""
function NQCBase.Structure.fractional_mass(sim::NQCDynamics.AbstractSimulation, index1::Int, index2::Int)
	return NQCBase.Structure.fractional_mass(sim.atoms,index1,index2)
end

function NQCBase.Structure.reduced_mass(sim::NQCDynamics.AbstractSimulation, index1::Int, index2::Int)
	return NQCBase.Structure.reduced_mass(sim.atoms,index1,index2)
end

function NQCBase.Structure.minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation;cutoff::Int=50)
	return NQCBase.Structure.minimum_distance_translation(config,ind1,ind2,simulation.cell;cutoff=cutoff)
end

function NQCBase.Structure.minimum_distance_translation(config::Any, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation;cutoff::Int=50)
	return NQCBase.Structure.minimum_distance_translation(NQCDynamics.get_positions(config),ind1,ind2,simulation.cell;cutoff=cutoff)
end

function NQCBase.Structure.pbc_distance(config::Matrix, ind1, ind2, sim::NQCDynamics.AbstractSimulation; args...)
	return NQCBase.Structure.pbc_distance(config,ind1,ind2,sim.cell; args...)
end

function NQCBase.Structure.pbc_distance(config::Any, ind1, ind2, sim::NQCDynamics.AbstractSimulation; args...)
	return NQCBase.Structure.pbc_distance(NQCDynamics.get_positions(config),ind1,ind2,sim.cell; args...)
end

function NQCBase.Structure.pbc_center_of_mass(config::Matrix, ind1, ind2, sim::NQCDynamics.AbstractSimulation; args...)
	return NQCBase.Structure.pbc_center_of_mass(config,ind1,ind2,sim.cell, sim.atoms; args...)
end

function NQCBase.Structure.pbc_center_of_mass(config::Any, ind1, ind2, sim::NQCDynamics.AbstractSimulation; args...)
	return NQCBase.Structure.pbc_center_of_mass(NQCDynamics.get_positions(config),ind1,ind2,sim.cell, sim.atoms; args...)
end

function NQCBase.Structure.velocity_center_of_mass(config::Matrix, ind1, ind2, sim::NQCDynamics.AbstractSimulation)
	return NQCBase.Structure.velocity_center_of_mass(config,ind1,ind2, sim.atoms)
end

function NQCBase.Structure.velocity_center_of_mass(config::Any, ind1, ind2, sim::NQCDynamics.AbstractSimulation)
	return NQCBase.Structure.velocity_center_of_mass(NQCDynamics.get_velocities(config),ind1,ind2, sim.atoms)
end

	


end