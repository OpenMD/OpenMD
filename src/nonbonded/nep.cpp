#include "nonbonded/nep.hpp"
#include "utils/simError.h"

#include <fstream>
#include <cmath>
#include <unistd.h>

using namespace std;

namespace OpenMD {

  NEP::NEP(const std::string &model_filename, int N_atoms,
           std::vector<RealType> cell, std::vector<std::string> atom_symbols,
           std::vector<RealType> positions, std::vector<RealType> masses)
    : nep(model_filename), model_filename(model_filename) {
    /**
       @brief Wrapper class for NEP3 in nep.cpp.
       @details This class wraps the setup functionality of the NEP3 class in
       nep.cpp.
       @param model_filename       Path to the NEP model (/path/nep.txt).
       @param N_atoms              The number of atoms in the structure.
       @param cell                 The cell vectors for the structure.
       @param atom_symbols         The atomic symbol for each atom in the structure.
       @param positions            The position for each atom in the structure.
       @param masses               The mass for each atom in the structure.
    */
  atoms.N = N_atoms;
  atoms.position = positions;
  atoms.mass = masses;
  atoms.cell = cell;
  atoms.type.resize(atoms.N);

  std::vector<std::string> model_atom_symbols =
    _getAtomSymbols(model_filename); // load atom symbols used in model
  _convertAtomTypeNEPIndex(atoms.N, atom_symbols, model_atom_symbols, atoms.type);
  }
  
  std::vector<RealType> NEP::getDescriptors() {
    /**
       @brief Get NEP descriptors.
       @details Calculates NEP descriptors for a given structure and NEP model.
    */
    std::vector<RealType> descriptor(atoms.N * nep.annmb.dim);
    nep.find_descriptor(atoms.type, atoms.cell, atoms.position, descriptor);
    return descriptor;
  }

  std::vector<RealType> NEP::getLatentSpace() {
    /**
       @brief Get the NEP latent space.
       @details Calculates the latent space of a NEP model, i.e., the
       activations after the first layer. 
    */
    std::vector<RealType> latent(atoms.N * nep.annmb.num_neurons1);
    nep.find_latent_space(atoms.type, atoms.cell, atoms.position, latent);
    return latent;
  }

  std::vector<RealType> NEP::getDipole() {
    /**
       @brief Get dipole.
       @details Calculates the dipole (a vector of length 3) for the whole box.
    */
    std::vector<RealType> dipole(3);
    nep.find_dipole(atoms.type, atoms.cell, atoms.position, dipole);
    return dipole;
  }

  void NEP::_getCenterOfMass(std::vector<RealType> center_of_mass) {
    /**
       @brief Computes the center of mass for current atoms object.
       @details Computes the center of mass (COM) for the structure
       with positions defined by atoms.position. The COM will be
       written as a three vector, with in the order [x_component,
       y_component, z_component].
       @param center_of_mass      Vector to hold the center of mass.
    */

    RealType total_mass = 0.0;
    for (int i = 0; i < atoms.N; i++) {
      // positions are in order [x1, ..., xN, y1, ..., yN, z1, ..., zN]
      center_of_mass[0] += atoms.position[i] * atoms.mass[i];
      center_of_mass[1] += atoms.position[i + atoms.N] * atoms.mass[i];
      center_of_mass[2] += atoms.position[i + 2 * atoms.N] * atoms.mass[i];
      total_mass += atoms.mass[i];
    }
    center_of_mass[0] /= total_mass;
    center_of_mass[1] /= total_mass;
    center_of_mass[2] /= total_mass;
  }

  std::vector<RealType> NEP::getDipoleGradient(RealType displacement,
					       int method,
					       RealType charge) {
  /**
   @brief Get dipole gradient through finite differences.
   @details Calculates the dipole gradient, a (N_atoms, 3, 3) tensor
   for the gradients dµ_k/dr_ij, for atom i, Cartesian direction j (x,
   y, z) and dipole moment component k. Mode is either forward
   difference (method=0) or first or second order central difference
   (method=1 and method=2 respectively).  Before computing the
   gradient the dipoles are corrected using the center of mass and the
   total system charge, supplied via the parameter `charge`.   
   @param displacement        Displacement in Å.
   @param method              0 or 1, corresponding to forward or central differences.
   @param charge              Total system charge, used to correct dipoles.
  */
    const int N_cartesian = 3;
    const int N_components = 3;
    const int values_per_atom = N_cartesian * N_components;
    
    std::vector<RealType> dipole_gradient(atoms.N * 3 * 3);
    
    // Compute the total mass - need this for the corrections to COM
    RealType M = 0.0;
    for (int i = 0; i < atoms.N; i++) {
      M += atoms.mass[i];
    }

    if (method == 0) {
      // For forward differences we save the first dipole with no displacement
      std::vector<RealType> dipole(3);
      nep.find_dipole(atoms.type, atoms.cell, atoms.position, dipole);
      
      // Correct original dipole
      std::vector<RealType> center_of_mass(3);
      _getCenterOfMass(center_of_mass);
      for (int i = 0; i < 3; i++) {
	dipole[i] += charge * center_of_mass[i];
      }

      // dipole vectors are zeroed in find_dipole, can be allocated here
      std::vector<RealType> dipole_forward(3);
      
      // Positions are in order [x1, ..., xN, y1, ..., yN, ...]
      int index = 0;
      RealType old_position = 0.0;
      RealType old_center_of_mass = 0.0;
      const RealType displacement_over_M = displacement / M;
      const RealType one_over_displacement = 1.0 / displacement;

      for (int i = 0; i < N_cartesian; i++) {
	for (int j = 0; j < atoms.N; j++) {
	  index = atoms.N * i + j; // idx of position to change
	  old_position = atoms.position[index];
	  atoms.position[index] += displacement;
	  // center of mass gest moved by +h/N*m_j in the same direction
	  old_center_of_mass = center_of_mass[i];
	  center_of_mass[i] += displacement_over_M * atoms.mass[j];
	  
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position, dipole_forward);

	  for (int k = 0; k < N_components; k++) {
	    dipole_gradient[i * N_components + j * values_per_atom + k] =
              ((dipole_forward[k] + charge * center_of_mass[k]) - dipole[k]) *
              one_over_displacement;
	  }
	  // Restore positions
	  atoms.position[index] = old_position;
	  center_of_mass[i] = old_center_of_mass;
	}
      }
    } else if (method == 1) {
      // For central differences we need both forwards and backwards
      // displacements

      // dipole vectors are zeroed in find_dipole, can be allocated
      // here
      std::vector<RealType> dipole_forward(3);
      std::vector<RealType> dipole_backward(3);
    
      // use center of mass to correct for permanent dipole
      std::vector<RealType> center_of_mass_forward(3);
      _getCenterOfMass(center_of_mass_forward);
      std::vector<RealType> center_of_mass_backward(center_of_mass_forward);

      // Positions are in order [x1, ..., xN, y1, ..., yN, ...]
      int index = 0;
      RealType old_position = 0.0;
      RealType old_center_of_mass = 0.0;
      const RealType displacement_over_M = displacement / M;
      const RealType one_over_two_displacements = 0.5 / displacement;
      
      for (int i = 0; i < N_cartesian; i++) {
	for (int j = 0; j < atoms.N; j++) {
	  index = atoms.N * i + j; // idx of position to change
	  old_position = atoms.position[index];
	  old_center_of_mass = center_of_mass_forward[i];
	  
	  // Forward displacement
	  atoms.position[index] += displacement;
	  // center of mass gest moved by +h/N in the same direction
	  center_of_mass_forward[i] += displacement_over_M * atoms.mass[j];
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position, dipole_forward);
	  
	  // Backwards displacement
	  atoms.position[index] -= 2 * displacement; // +h - 2h = -h
	  center_of_mass_backward[i] -= displacement_over_M * atoms.mass[j];
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position, dipole_backward);
	  
	  for (int k = 0; k < N_components; k++) {
	    dipole_gradient[i * N_components + j * values_per_atom + k] =
	      ((dipole_forward[k] + charge * center_of_mass_forward[k]) -
	       (dipole_backward[k] + charge * center_of_mass_backward[k])) *
	      one_over_two_displacements;
	  }
	  // Restore positions
	  atoms.position[index] = old_position;
	  center_of_mass_forward[i] = old_center_of_mass;
	  center_of_mass_backward[i] = old_center_of_mass;
	}
      }
    } else if (method == 2) {
      // Second order central differences
      // Need to compute four dipoles for each structure, yielding an error O(h^4)
      // Coefficients are defined here:
      // https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference
      
      // dipole vectors are zeroed in find_dipole, can be allocated here
      std::vector<RealType> dipole_forward_one_h(3);
      std::vector<RealType> dipole_forward_two_h(3);
      std::vector<RealType> dipole_backward_one_h(3);
      std::vector<RealType> dipole_backward_two_h(3);

      // use center of mass to correct for permanent dipole
      std::vector<RealType> center_of_mass_forward_one_h(3);
      _getCenterOfMass(center_of_mass_forward_one_h);
      std::vector<RealType> center_of_mass_forward_two_h(center_of_mass_forward_one_h);
      std::vector<RealType> center_of_mass_backward_one_h(center_of_mass_forward_one_h);
      std::vector<RealType> center_of_mass_backward_two_h(center_of_mass_forward_one_h);

      // Positions are in order [x1, ..., xN, y1, ..., yN, ...]
      int index = 0;
      RealType old_position = 0.0;
      RealType old_center_of_mass = 0.0;
      const RealType displacement_over_M = displacement / M;
      const RealType one_over_displacement =
        1.0 / displacement; // coefficients are scaled properly
      
      const RealType c0 = -1.0 / 12.0; // coefficient for 2h
      const RealType c1 = 2.0 / 3.0;   // coefficient for h
      for (int i = 0; i < N_cartesian; i++) {
	for (int j = 0; j < atoms.N; j++) {
	  index = atoms.N * i + j; // idx of position to change
	  old_position = atoms.position[index];
	  old_center_of_mass = center_of_mass_forward_one_h[i];
	  
	  // --- Forward displacement
	  // Step one displacement forward
	  atoms.position[index] += displacement; // + h
	  center_of_mass_forward_one_h[i] +=
            displacement_over_M * atoms.mass[j]; // center of mass gets moved by
                                                 // +h/N in the same direction
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position,
			  dipole_forward_one_h);
	  // Step two displacements forward
	  atoms.position[index] += displacement; // + 2h total
	  center_of_mass_forward_two_h[i] +=
            2 * displacement_over_M *
            atoms.mass[j]; // center of mass gest moved by
                           // +2h/N in the same direction
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position,
			  dipole_forward_two_h);

	  // --- Backwards displacement
	  atoms.position[index] -= 3 * displacement; // 2h - 3h = -h
	  center_of_mass_backward_one_h[i] -= displacement_over_M * atoms.mass[j];
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position,
			  dipole_backward_one_h);

	  atoms.position[index] -= displacement; // -h - h = -2h
	  center_of_mass_backward_two_h[i] -=
            2 * displacement_over_M * atoms.mass[j];
	  nep.find_dipole(atoms.type, atoms.cell, atoms.position,
			  dipole_backward_two_h);
	  
	  for (int k = 0; k < N_components; k++) {
	    dipole_gradient[i * N_components + j * values_per_atom + k] =
              (c0 * (dipole_forward_two_h[k] +
                     charge * center_of_mass_forward_two_h[k]) +
               c1 * (dipole_forward_one_h[k] +
                     charge * center_of_mass_forward_one_h[k]) -
               c1 * (dipole_backward_one_h[k] +
                     charge * center_of_mass_backward_one_h[k]) -
               c0 * (dipole_backward_two_h[k] +
                     charge * center_of_mass_backward_two_h[k])) *
              one_over_displacement;
	  }
	  // Restore positions
	  atoms.position[index] = old_position;
	  center_of_mass_forward_one_h[i] = old_center_of_mass;
	  center_of_mass_forward_two_h[i] = old_center_of_mass;
	  center_of_mass_backward_one_h[i] = old_center_of_mass;
	  center_of_mass_backward_two_h[i] = old_center_of_mass;
	}
      }
    }
    // dipole gradient component d_x refers to cartesian direction x
    // x1 refers to x position of atom 1
    // order: [dx_x1, dy_x1, dz_x1,
    //         dx_y1, dy_y1, dz_y1,
    //         dx_z1, dy_z1, dz_z1,
    //         ...
    //         dx_zN, dy_zN, dz_zN]
    return dipole_gradient;
  }

  std::vector<RealType> NEP::getPolarizability() {

    /**
       @brief Get polarizability.
       @details Calculates the polarizability (a 2-tensor, represented as vector of
       length 9) for the whole box.
       Output order is pxx pxy pxz pyx pyy pyz pzx pzy pzz.
    */
    std::vector<RealType> polarizability(6);
    nep.find_polarizability(atoms.type, atoms.cell, atoms.position, polarizability);
    return polarizability;
}

  std::vector<RealType> NEP::getPolarizabilityGradient(RealType displacement, std::vector<int> components) {
    /**
       @brief Get polarizability gradient through finite differences.
       @details Calculates the polarizability gradient, a (N_atoms, 3, 6) tensor for the
       gradients dp_k/dr_ij, for atom i, Cartesian direction j (x, y, z) and polarizability
       component k, using second order finite differences.
       @param displacement        Displacement in Å.
    */
    const int N_cartesian = 3;
    const int N_components = 6;
    std::vector<RealType> polarizability_gradient(atoms.N * N_cartesian * N_components);
    
    // Second order central differences
    // Need to compute four dipoles for each structure, yielding an error O(h^4)
    // Coefficients are defined here:
    // https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference
    
    // Compute polarizability for forward displacement
    const int values_per_atom = N_cartesian * N_components;
    
    // polarizability vectors are zeroed in find_polarizability, can be allocated here
    std::vector<RealType> polarizability_forward_one_h(6);
    std::vector<RealType> polarizability_forward_two_h(6);
    std::vector<RealType> polarizability_backward_one_h(6);
    std::vector<RealType> polarizability_backward_two_h(6);
    
    // Positions are in order [x1, ..., xN, y1, ..., yN, ...]
    int index = 0;
    RealType old_position = 0.0;
    const RealType one_over_displacement =
      1.0 / displacement; // coefficients are scaled properly
    
    const RealType c0 = -1.0 / 12.0; // coefficient for 2h
    const RealType c1 = 2.0 / 3.0;   // coefficient for h
    for (int i = 0; i < N_cartesian; i++) {
      if (std::count(components.begin(), components.end(), i) == 0) {
	// if i is not present in components to calculate, skip this iteration
	continue;
      }
      for (int j = 0; j < atoms.N; j++) {
	// Positions have order x1, ..., xN, y1,...,yN, z1,...,zN
	index = atoms.N * i + j; // idx of position to change
	old_position = atoms.position[index];
	
	// --- Forward displacement
	// Step one displacement forward
	atoms.position[index] += displacement; // + h
	nep.find_polarizability(atoms.type, atoms.cell, atoms.position,
				polarizability_forward_one_h);
	// Step two displacements forward
	atoms.position[index] += displacement; // + 2h total
	nep.find_polarizability(atoms.type, atoms.cell, atoms.position,
				polarizability_forward_two_h);
	
	// --- Backwards displacement
	atoms.position[index] -= 3 * displacement; // 2h - 3h = -h
	nep.find_polarizability(atoms.type, atoms.cell, atoms.position,
				polarizability_backward_one_h);
	
	atoms.position[index] -= displacement; // -h - h = -2h
	nep.find_polarizability(atoms.type, atoms.cell, atoms.position,
				polarizability_backward_two_h);
	
	// --- Compute gradient ---
	for (int k = 0; k < N_components; k++) {
	  polarizability_gradient[i * N_components + j * values_per_atom + k] =
            (c0 * polarizability_forward_two_h[k] +
             c1 * polarizability_forward_one_h[k] -
             c1 * polarizability_backward_one_h[k] -
             c0 * polarizability_backward_two_h[k]) *
            one_over_displacement;
	}
	// Restore positions
	atoms.position[index] = old_position;
      }
    }
    // polarizability gradient component p_x refers to cartesian direction x
    // x1 refers to x position of atom 1
    // order: [pxx_x1, pyy_x1, pzz_x1, pxy_x1, pyz_x1, pzx_x1
    //         pxx_y2, pyy_y1, pzz_y1, pxy_y1, pyz_y1, pzx_y1
    //         pxx_z2, pyy_z1, pzz_z1, pxy_z1, pyz_z1, pzx_z1
    //         ...
    //         pxx_zN, pyy_zN, pzz_zN, pxy_zN, pyz_zN, pzx_zN]
    return polarizability_gradient;
  }
  
  std::tuple<std::vector<RealType>, std::vector<RealType>, std::vector<RealType>>
  NEP::getPotentialForcesAndVirials() {
    /**
       @brief Calculate potential, forces and virials.
       @details Calculates potential energy, forces and virials for a given
       structure and NEP model.
    */
    std::vector<RealType> potential(atoms.N);
    std::vector<RealType> force(atoms.N * 3);
    std::vector<RealType> virial(atoms.N * 9);
    nep.compute(atoms.type, atoms.cell, atoms.position, potential, force, virial);
    return std::make_tuple(potential, force, virial);
  }
  
  std::vector<std::string> NEP::_getAtomSymbols(std::string model_filename) {
    /**
       @brief Fetches atomic symbols
       @details This function fetches the atomic symbols from the header of a NEP
       model. These are later used to ensure consistent indices for the atom types.
       @param model_filename Path to the NEP model (/path/nep.txt).
   */
    std::ifstream input_potential(model_filename);
    if (!input_potential.is_open()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
	       "NEP: Cannot open input_potential file: %s\n",
	       model_filename.c_str());
      painCave.isFatal = 1;
      simError();
    }
    std::string potential_name;
    input_potential >> potential_name;
    int number_of_types;
    input_potential >> number_of_types;
    std::vector<std::string> atom_symbols(number_of_types);
    for (int n = 0; n < number_of_types; ++n) {
      input_potential >> atom_symbols[n];
    }
    input_potential.close();
    return atom_symbols;
  }
  
  void NEP::_convertAtomTypeNEPIndex(int N,
				     std::vector<std::string> atom_symbols,
				     std::vector<std::string> model_atom_symbols,
				     std::vector<int> &type) {
    /**
       @brief Converts atom species to NEP index.
       @details Converts atomic species to indicies, which are used in NEP.
       @param atom_symbols        List of atom symbols.
       @param model_atom_symbols  List of atom symbols used in the model.
       @param type                List of indices corresponding to atom type.
    */
    for (int n = 0; n < N; n++) {
      // Convert atom type to index for consistency with nep.txt
      // (model_filename)
      std::string atom_symbol = atom_symbols[n];
      bool is_allowed_element = false;
      for (int t = 0; (unsigned)t < model_atom_symbols.size(); ++t) {
	if (atom_symbol == model_atom_symbols[t]) {
	  type[n] = t;
	  is_allowed_element = true;
	}
      }
      if (!is_allowed_element) {
	snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
		 "NEP: Atom type: %s inot used in the given NEP model.\n",
		 atom_symbols[n].c_str());
	painCave.isFatal = 1;
	simError();
      }
    }
  }
  
  void NEP::setPositions(std::vector<RealType> positions) {
    /**
       @brief Sets positions.
       @details Sets the positions of the atoms object.
       Also updates the center of mass.
       @param positions           List of positions.
    */
    for (int i = 0; i < atoms.N * 3; i++) {
      atoms.position[i] = positions[i];
    }
  }
  
  void NEP::setCell(std::vector<RealType> cell) {
    /**
       @brief Sets cell.
       @details Sets the cell of the atoms object.
       @param Cell                Cell vectors.
    */
    for (int i = 0; i < 9; i++) {
      atoms.cell[i] = cell[i];
    }
  }
  
  void NEP::setMasses(std::vector<RealType> masses) {
    /**
       @brief Sets masses.
       @details Sets the masses of the atoms object.
       @param Cell                Atom masses.
    */
    for (int i = 0; i < atoms.N; i++) {
      atoms.mass[i] = masses[i];
    }
  }
  
  void NEP::setSymbols(std::vector<std::string> atom_symbols) {
    /**
       @brief Sets symbols.
       @details Sets the symbols of the atoms object from the ones used in the
       model.
       @param atom_symbols        List of symbols.
    */
    std::vector<std::string> model_atom_symbols =
      _getAtomSymbols(model_filename); // load atom symbols used in model
    _convertAtomTypeNEPIndex(atoms.N, atom_symbols, model_atom_symbols, atoms.type);
  }
}
