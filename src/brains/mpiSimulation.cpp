#ifdef IS_MPI
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "brains/mpiSimulation.hpp"
#include "utils/simError.h"
#include "UseTheForce/fortranWrappers.hpp"
#include "math/randomSPRNG.hpp"

mpiSimulation* mpiSim;

mpiSimulation::mpiSimulation(SimInfo* the_entryPlug)
{
  entryPlug = the_entryPlug;
  parallelData = new mpiSimData;
  
  MPI_Comm_size(MPI_COMM_WORLD, &(parallelData->nProcessors) );
  parallelData->myNode = worldRank;

  MolToProcMap = new int[entryPlug->n_mol];
  MolComponentType = new int[entryPlug->n_mol];
  AtomToProcMap = new int[entryPlug->n_atoms];
  GroupToProcMap = new int[entryPlug->ngroup];

  mpiSim = this;
  wrapMeSimParallel( this );
}


mpiSimulation::~mpiSimulation(){
  
  delete[] MolToProcMap;
  delete[] MolComponentType;
  delete[] AtomToProcMap;
  delete[] GroupToProcMap;

  delete parallelData;
  // perhaps we should let fortran know the party is over.
  
}

void mpiSimulation::divideLabor( ){

  int nComponents;
  MoleculeStamp** compStamps;
  randomSPRNG *myRandom;
  int* componentsNmol;
  int* AtomsPerProc;
  int* GroupsPerProc;

  double numerator;
  double denominator;
  double precast;
  double x, y, a;
  int old_atoms, add_atoms, new_atoms;
  int old_groups, add_groups, new_groups;

  int nTarget;
  int molIndex, atomIndex, groupIndex;
  int done;
  int i, j, loops, which_proc;
  int nmol_global, nmol_local;
  int ngroups_global, ngroups_local;
  int natoms_global, natoms_local;
  int ncutoff_groups, nAtomsInGroups;
  int local_index;
  int baseSeed = entryPlug->getSeed();
  CutoffGroupStamp* cg;

  nComponents = entryPlug->nComponents;
  compStamps = entryPlug->compStamps;
  componentsNmol = entryPlug->componentsNmol;
  AtomsPerProc = new int[parallelData->nProcessors];
  GroupsPerProc = new int[parallelData->nProcessors];
  
  parallelData->nAtomsGlobal = entryPlug->n_atoms;
  parallelData->nBondsGlobal = entryPlug->n_bonds;
  parallelData->nBendsGlobal = entryPlug->n_bends;
  parallelData->nTorsionsGlobal = entryPlug->n_torsions;
  parallelData->nSRIGlobal = entryPlug->n_SRI;
  parallelData->nGroupsGlobal = entryPlug->ngroup;
  parallelData->nMolGlobal = entryPlug->n_mol;

  if (parallelData->nProcessors > parallelData->nMolGlobal) {
    sprintf( painCave.errMsg,
             "nProcessors (%d) > nMol (%d)\n"
             "\tThe number of processors is larger than\n"
             "\tthe number of molecules.  This will not result in a \n"
             "\tusable division of atoms for force decomposition.\n"
             "\tEither try a smaller number of processors, or run the\n"
             "\tsingle-processor version of OOPSE.\n",
             parallelData->nProcessors, parallelData->nMolGlobal );
    painCave.isFatal = 1;
    simError();
  }

  myRandom = new randomSPRNG( baseSeed );  


  a = 3.0 * (double)parallelData->nMolGlobal / (double)parallelData->nAtomsGlobal;

  // Initialize things that we'll send out later:
  for (i = 0; i < parallelData->nProcessors; i++ ) {
    AtomsPerProc[i] = 0;
    GroupsPerProc[i] = 0;
  }
  for (i = 0; i < parallelData->nMolGlobal; i++ ) {
    // default to an error condition:
    MolToProcMap[i] = -1;
    MolComponentType[i] = -1;
  }
  for (i = 0; i < parallelData->nAtomsGlobal; i++ ) {
    // default to an error condition:
    AtomToProcMap[i] = -1;
  }
  for (i = 0; i < parallelData->nGroupsGlobal; i++ ) {
    // default to an error condition:
    GroupToProcMap[i] = -1;
  }
    
  if (parallelData->myNode == 0) {
    numerator = (double) entryPlug->n_atoms;
    denominator = (double) parallelData->nProcessors;
    precast = numerator / denominator;
    nTarget = (int)( precast + 0.5 );

    // Build the array of molecule component types first
    molIndex = 0;
    for (i=0; i < nComponents; i++) {
      for (j=0; j < componentsNmol[i]; j++) {        
        MolComponentType[molIndex] = i;
        molIndex++;
      }
    }

    atomIndex = 0;
    groupIndex = 0;

    for (i = 0; i < molIndex; i++ ) {

      done = 0;
      loops = 0;

      while( !done ){
        loops++;
        
        // Pick a processor at random

        which_proc = (int) (myRandom->getRandom() * parallelData->nProcessors);

        // How many atoms does this processor have?
        
        old_atoms = AtomsPerProc[which_proc];
        add_atoms = compStamps[MolComponentType[i]]->getNAtoms();
        new_atoms = old_atoms + add_atoms;

        old_groups = GroupsPerProc[which_proc];
        ncutoff_groups = compStamps[MolComponentType[i]]->getNCutoffGroups(); 
        nAtomsInGroups = 0;
        for (j = 0; j < ncutoff_groups; j++) {
          cg = compStamps[MolComponentType[i]]->getCutoffGroup(j);
          nAtomsInGroups += cg->getNMembers();
        }
        add_groups = add_atoms - nAtomsInGroups + ncutoff_groups;        
        new_groups = old_groups + add_groups;

        // If we've been through this loop too many times, we need
        // to just give up and assign the molecule to this processor
        // and be done with it. 
        
        if (loops > 100) {          
          sprintf( painCave.errMsg,
                   "I've tried 100 times to assign molecule %d to a "
                   " processor, but can't find a good spot.\n"  
                   "I'm assigning it at random to processor %d.\n",
                   i, which_proc);
          painCave.isFatal = 0;
          simError();
          
          MolToProcMap[i] = which_proc;
          AtomsPerProc[which_proc] += add_atoms;
          for (j = 0 ; j < add_atoms; j++ ) {
	    AtomToProcMap[atomIndex] = which_proc;
	    atomIndex++;
          }
          GroupsPerProc[which_proc] += add_groups;
          for (j=0; j < add_groups; j++) {
            GroupToProcMap[groupIndex] = which_proc;
            groupIndex++;
          }
          done = 1;
          continue;
        }
    
        // If we can add this molecule to this processor without sending
        // it above nTarget, then go ahead and do it:
    
        if (new_atoms <= nTarget) {
          MolToProcMap[i] = which_proc;
          AtomsPerProc[which_proc] += add_atoms;
          for (j = 0 ; j < add_atoms; j++ ) {
	    AtomToProcMap[atomIndex] = which_proc;
	    atomIndex++;
          }
          GroupsPerProc[which_proc] += add_groups;
          for (j=0; j < add_groups; j++) {
            GroupToProcMap[groupIndex] = which_proc;
            groupIndex++;
          }
          done = 1;
          continue;
        }


        // The only situation left is when new_atoms > nTarget.  We
        // want to accept this with some probability that dies off the
        // farther we are from nTarget

        // roughly:  x = new_atoms - nTarget
        //           Pacc(x) = exp(- a * x)
        // where a = penalty / (average atoms per molecule)

        x = (double) (new_atoms - nTarget);
        y = myRandom->getRandom();
      
        if (y < exp(- a * x)) {
          MolToProcMap[i] = which_proc;
          AtomsPerProc[which_proc] += add_atoms;
          for (j = 0 ; j < add_atoms; j++ ) {
	    AtomToProcMap[atomIndex] = which_proc;
	    atomIndex++;
           }
          GroupsPerProc[which_proc] += add_groups;
          for (j=0; j < add_groups; j++) {
            GroupToProcMap[groupIndex] = which_proc;
            groupIndex++;
          }
          done = 1;
          continue;
        } else {
          continue;
        }       
        
      }
    }


    // Spray out this nonsense to all other processors:

    //std::cerr << "node 0 mol2proc = \n";
    //for (i = 0; i < parallelData->nMolGlobal; i++) 
    //  std::cerr << i << "\t" << MolToProcMap[i] << "\n";

    MPI_Bcast(MolToProcMap, parallelData->nMolGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(AtomToProcMap, parallelData->nAtomsGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(GroupToProcMap, parallelData->nGroupsGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(MolComponentType, parallelData->nMolGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(AtomsPerProc, parallelData->nProcessors,
	      MPI_INT, 0, MPI_COMM_WORLD);    

    MPI_Bcast(GroupsPerProc, parallelData->nProcessors,
	      MPI_INT, 0, MPI_COMM_WORLD);    
  } else {

    // Listen to your marching orders from processor 0:
    
    MPI_Bcast(MolToProcMap, parallelData->nMolGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(AtomToProcMap, parallelData->nAtomsGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(GroupToProcMap, parallelData->nGroupsGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(MolComponentType, parallelData->nMolGlobal, 
	      MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(AtomsPerProc, parallelData->nProcessors,
	      MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(GroupsPerProc, parallelData->nProcessors,
	      MPI_INT, 0, MPI_COMM_WORLD);


  }

  // Let's all check for sanity:

  nmol_local = 0;
  for (i = 0 ; i < parallelData->nMolGlobal; i++ ) {
    if (MolToProcMap[i] == parallelData->myNode) {
      nmol_local++;
    }
  }

  natoms_local = 0;
  for (i = 0; i < parallelData->nAtomsGlobal; i++) {
    if (AtomToProcMap[i] == parallelData->myNode) {
      natoms_local++;      
    }
  }

  ngroups_local = 0;
  for (i = 0; i < parallelData->nGroupsGlobal; i++) {
    if (GroupToProcMap[i] == parallelData->myNode) {
      ngroups_local++;      
    }
  }

  MPI_Allreduce(&nmol_local,&nmol_global,1,MPI_INT,MPI_SUM, 
		MPI_COMM_WORLD);

  MPI_Allreduce(&natoms_local,&natoms_global,1,MPI_INT,
		MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&ngroups_local,&ngroups_global,1,MPI_INT,
		MPI_SUM, MPI_COMM_WORLD);
  
  if( nmol_global != entryPlug->n_mol ){
    sprintf( painCave.errMsg,
             "The sum of all nmol_local, %d, did not equal the "
             "total number of molecules, %d.\n",
             nmol_global, entryPlug->n_mol );
    painCave.isFatal = 1;
    simError();
  }
  
  if( natoms_global != entryPlug->n_atoms ){
    sprintf( painCave.errMsg,
             "The sum of all natoms_local, %d, did not equal the "
             "total number of atoms, %d.\n",
             natoms_global, entryPlug->n_atoms );
    painCave.isFatal = 1;
    simError();
  }

  if( ngroups_global != entryPlug->ngroup ){
    sprintf( painCave.errMsg,
             "The sum of all ngroups_local, %d, did not equal the "
             "total number of cutoffGroups, %d.\n",
             ngroups_global, entryPlug->ngroup );
    painCave.isFatal = 1;
    simError();
  }

  sprintf( checkPointMsg,
	   "Successfully divided the molecules among the processors.\n" );
  MPIcheckPoint();

  parallelData->nMolLocal = nmol_local;
  parallelData->nAtomsLocal = natoms_local;
  parallelData->nGroupsLocal = ngroups_local;

  globalAtomIndex.resize(parallelData->nAtomsLocal);
  globalToLocalAtom.resize(parallelData->nAtomsGlobal);
  local_index = 0;
  for (i = 0; i < parallelData->nAtomsGlobal; i++) {
    if (AtomToProcMap[i] == parallelData->myNode) {
      globalAtomIndex[local_index] = i;

      globalToLocalAtom[i] = local_index;
      local_index++;
      
    }
    else
       globalToLocalAtom[i] = -1;
  }

  globalGroupIndex.resize(parallelData->nGroupsLocal);
  globalToLocalGroup.resize(parallelData->nGroupsGlobal);
  local_index = 0;
  for (i = 0; i < parallelData->nGroupsGlobal; i++) {
    if (GroupToProcMap[i] == parallelData->myNode) {
      globalGroupIndex[local_index] = i;

      globalToLocalGroup[i] = local_index;
      local_index++;
      
    }
    else
       globalToLocalGroup[i] = -1;
  }

  globalMolIndex.resize(parallelData->nMolLocal);
  globalToLocalMol.resize(parallelData->nMolGlobal);  
  local_index = 0;
  for (i = 0; i < parallelData->nMolGlobal; i++) {
    if (MolToProcMap[i] == parallelData->myNode) {
      globalMolIndex[local_index] = i;
      globalToLocalMol[i] = local_index;
      local_index++;
    }
    else
      globalToLocalMol[i] = -1;
  }
  
}


void mpiSimulation::mpiRefresh( void ){

  int isError, i;
  int *localToGlobalAtomIndex = new int[parallelData->nAtomsLocal];
  int *localToGlobalGroupIndex = new int[parallelData->nGroupsLocal];

  // Fortran indexing needs to be increased by 1 in order to get the 2
  // languages to not barf

  for(i = 0; i < parallelData->nAtomsLocal; i++) 
    localToGlobalAtomIndex[i] = globalAtomIndex[i] + 1;

  for(i = 0; i < parallelData->nGroupsLocal; i++) 
    localToGlobalGroupIndex[i] = globalGroupIndex[i] + 1;
  
  isError = 0;

  setFsimParallel( parallelData, 
                   &(parallelData->nAtomsLocal), localToGlobalAtomIndex, 
                   &(parallelData->nGroupsLocal),  localToGlobalGroupIndex, 
                   &isError );

  if( isError ){

    sprintf( painCave.errMsg,
	     "mpiRefresh errror: fortran didn't like something we gave it.\n" );
    painCave.isFatal = 1;
    simError();
  }

  delete[] localToGlobalGroupIndex;
  delete[] localToGlobalAtomIndex;


  sprintf( checkPointMsg,
	   " mpiRefresh successful.\n" );
  MPIcheckPoint();
}
 

#endif // is_mpi
