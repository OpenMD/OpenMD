#ifndef __MPISIMULATION__
#define __MPISIMULATION__

#include "brains/SimInfo.hpp"
#include "types/MakeStamps.hpp"
#define __C
#include "UseTheForce/mpiComponentPlan.h"

#include "UseTheForce/DarkSide/simParallel_interface.h"

class mpiSimulation{
public:

  mpiSimulation(SimInfo* the_entryPlug);
  ~mpiSimulation();
  
  void divideLabor();
  
  int  getMyNode( void )         { return parallelData->myNode; }
  int  getNProcessors( void )    { return parallelData->nProcessors; }
  int  getNMolLocal( void )      { return parallelData->nMolLocal; }
  int  getNMolGlobal( void )     { return parallelData->nMolGlobal; }
  int  getNAtomsLocal( void )    { return parallelData->nAtomsLocal; }
  int  getNAtomsGlobal( void )   { return parallelData->nAtomsGlobal; }
  int  getNGroupsLocal( void )   { return parallelData->nGroupsLocal; }
  int  getNGroupsGlobal( void )  { return parallelData->nGroupsGlobal; }
  int* getAtomToProcMap( void )  { return AtomToProcMap; }
  int* getGroupToProcMap( void ) { return GroupToProcMap; }
  int* getMolToProcMap( void )   { return MolToProcMap; }
  int* getMolComponentType(void) { return MolComponentType; }
  vector<int> getGlobalAtomIndex(void) {return globalAtomIndex; }
  vector<int> getGlobalGroupIndex(void) {return globalGroupIndex; }
  vector<int> getGlobalMolIndex(void) {return globalMolIndex;}
  int getGlobalToLocalMol(int globalIndex) {return globalToLocalMol[globalIndex];}
  int getGlobalToLocalAtom(int globalIndex) {return globalToLocalAtom[globalIndex];}
  int getGlobalToLocalGroup(int globalIndex) {return globalToLocalGroup[globalIndex];}

  // call at the begining and after load balancing
  
  void mpiRefresh( void );

protected:
  SimInfo* entryPlug;
  mpiSimData* parallelData;
  int *MolToProcMap;
  int *MolComponentType;
  int *AtomToProcMap;
  int *AtomType;
  int *GroupToProcMap;
  vector<int> globalAtomIndex;
  vector<int> globalMolIndex;
  vector<int> globalGroupIndex;

  vector<int> globalToLocalMol;
  vector<int> globalToLocalAtom;
  vector<int> globalToLocalGroup;

  // int *myIdents; // is needed by Cpp only. It tells the molecule which stamp it is.

};

/**
   The following pointer is the global declaration of the mpiSim
   object created when the mpiSimulation creation routine is
   called. Every one who includes the header file will then have
   access to all of the routines in mpiSimulation class. 
*/ 

extern mpiSimulation* mpiSim;

#endif
