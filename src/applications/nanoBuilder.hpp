#ifndef __NANOBUILDER_H__
#define __NANOBUILDER_H__

#include <vector>

#define DEFAULT 1
#define BUILD_CORE_SHELL 2
#define BUILD_CORE_SHELL_VACANCY 3
#define BUILD_RANDOM_PARTICLE 4
#define BUILD_NMOL_PARTICLE 5


namespace nano{
  struct Mols
  {
    double pos[3]; // Center of Molecule position
    MoleculeStamp* myStamp;
    int isCore;
    int isShell;
    int isVacancy;
  };
};



// returns 1 if successful, 0 otherwise
class nanoBuilder{

public:
  nanoBuilder(int &hasError);
  ~nanoBuilder();

  int buildNanoParticle(void);

  
private: 

  //Support Fncs.
  int sanityCheck( void );
  void buildVacancies(void);
  void orientationMunger(double A[3][3]);
  void placeRandom(int totalMol);
  // Builder Fncs
  void buildWithVacancies(double dist, double pos[3]);
  void buildRandomlyMixed(double dist, double pos[3]);
  void buildWithCoreShell(double dist,double pos[3]);
  void buildNmolParticle(double dist, double pos[3]);

  //Logicals
  int isRandom;
  int hasVacancies;
  int latticeType;
  int buildNmol;
  int buildType;
  int coreHasOrientation;
  int shellHasOrientation;

  // Int values
  int nVacancies;
  int nInterface;
  int nCells;
  int nMol;
  int nCoreModelAtoms;  // Number of atoms in Core model
  int nShellModelAtoms; // Number of atoms in Shell model
  int moleculeCount;
  int coreAtomCount;  // Count number of core atoms in system.
  int shellAtomCount; // Count number of shell atoms in system.
  int atomCount;
  int totalMolecules; // Total number of molecules
  int nCoreMolecules; // Total number of core molecules.
  int nShellMolecules; // Total number of shell molecules.
  int maxModelNatoms;
  double particleRadius;
  double coreRadius;
  double vacancyFraction;
  double vacancyRadius;
  double shellRadius;
  double latticeSpacing;
  double soluteX; //Mole fraction for randomly mixed nanoparticle.

  //Vector to store atoms while were building nanoparticle.
  std::vector<nano::Mols> moleculeVector;
  std::vector<int> vacancyInterface;

  
  nano::Mols myMol;
  
  MoleculeStamp* coreStamp;
  MoleculeStamp* shellStamp;


 

};



#endif // __NANOBUILDER_H__
