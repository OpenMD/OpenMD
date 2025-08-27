#ifndef BRAINS_NEP_HPP
#define BRAINS_NEP_HPP

#include <config.h>
#include "NEP_CPU/src/nep.h"

using namespace std;

namespace OpenMD {

  struct NEPAtoms {
    int N;
    std::vector<int> type;
    std::vector<RealType> cell, position, mass;
  };
  
  class NEP {
  public:
    
    NEP(const std::string &model_filename, int N_atoms,
	std::vector<RealType> box, std::vector<std::string> atom_symbols,
	std::vector<RealType> positions, std::vector<RealType> masses);
    
    std::vector<RealType> getDescriptors();
    std::vector<RealType> getDipole();
    std::vector<RealType> getDipoleGradient(RealType displacement, int method,
					    RealType charge);
    std::vector<RealType> getPolarizability();
    std::vector<RealType> getPolarizabilityGradient(RealType displacement,
						    std::vector<int> components);
    std::vector<RealType> getLatentSpace();
    std::tuple<std::vector<RealType>, std::vector<RealType>, std::vector<RealType>> getPotentialForcesAndVirials();
    std::vector<std::string> _getAtomSymbols(std::string model_filename);
    void _convertAtomTypeNEPIndex(int N, std::vector<std::string> atom_symbols,
				  std::vector<std::string> model_atom_symbols,
				  std::vector<int> &type);
    void _getCenterOfMass(std::vector<RealType> center_of_mass);
    void setPositions(std::vector<RealType> positions);
    void setCell(std::vector<RealType> cell);
    void setMasses(std::vector<RealType> masses);
    void setSymbols(std::vector<std::string> atom_symbols);

  private:
    NEP3 nep;
    struct NEPAtoms atoms;
    std::string model_filename;
  };
}
#endif
