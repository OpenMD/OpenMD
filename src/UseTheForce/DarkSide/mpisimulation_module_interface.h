/*
 *  mpisimulation_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#define __C
#include "config.h"

extern "C" {
  
  typedef void (setFortranMPI)( mpiSimData* the_mpiPlug,
                                int* nLocal, 
                                int* globalAtomIndex,
                                int* nGroupsLocal,
                                int* globalGroupIndex,
                                int* isError ){
    F90_FUNC(setfortranmpi, SETFORTRANMPI)(the_mpiPlug,
                                           nLocal,
                                           globalAtomIndex,
                                           nGroupsLocal,
                                           globalGroupIndex);
    
    
  }
}