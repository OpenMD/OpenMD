#ifndef __SIMSETUP_H__
#define __SIMSETUP_H__
#include <string>
#include "StringUtils.hpp"
#include "MakeStamps.hpp"
#include "Globals.hpp"
#include "ForceFields.hpp"
#include "SimInfo.hpp"
#include "ReadWrite.hpp"
#include "AllIntegrator.hpp"

using namespace std;
// this routine is defined in BASS_interface.cpp
extern void set_interface_stamps( MakeStamps* ms, Globals* g );

string getPrefix(const string& str );

class SimSetup{

public:
  SimSetup();
  ~SimSetup();

  void setSimInfo( SimInfo* the_info ) { info = the_info; }
  void setSimInfo( SimInfo* the_info, int theNinfo );
  void suspendInit( void ) { initSuspend = true; }  
  void parseFile( char* fileName );
  void createSim( void );

  ForceFields* getForceField() {return the_ff;}
private:

#ifdef IS_MPI
  void receiveParse(void);
#endif

  void gatherInfo( void );
  void sysObjectsCreation( void );
  void finalInfoCheck( void );
  void initSystemCoords( void );
  void makeOutNames(void);
  void makeIntegrator(void);
  void initFortran(void);
  void makeMinimizer(void);

  void createFF( void );
  void compList( void );
  void calcSysValues( void );
  void makeSysArrays( void );

#ifdef IS_MPI
  void mpiMolDivide( void );

  int* mol2proc;
  int* molCompType;

#endif //is_mpi

  void initFromMetaDataFile( void );
  void makeMolecules( void );
  void makeElement( double x, double y, double z );

  int ensembleCase;
  int ffCase;

  MakeStamps* stamps;
  Globals* globals;
  char* inFileName;

  SimInfo* info;
  int isInfoArray;
  int nInfo;
  
  bool initSuspend;

  int n_components;
  int globalAtomCounter;
  int globalMolCounter;

  char force_field[100];
  char forcefield_variant[100];
  char ensemble[100];
  Component** the_components;

  int* components_nmol;
  MoleculeStamp** comp_stamps; //the stamps matching the components
  int tot_nmol;
  int tot_atoms;
  int tot_groups;
  int tot_rigid;
  int tot_bonds;
  int tot_bends;
  int tot_torsions;
  int tot_SRI;

  ForceFields* the_ff;

  // needed by makeElement

  int current_mol;
  int current_comp_mol;
  int current_comp;
  int current_atom_ndx;
  short int has_forcefield_variant;

  vector<int> globalAtomIndex;
  vector<int> globalGroupIndex;
  void setupZConstraint(SimInfo& theInfo);  //setup parameters for zconstraint method

};
#endif
