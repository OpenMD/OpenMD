/*
 *  mpisimulation_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef USETHEFORCE_DARKSIDE_MPISIMULATION_INTERFACE_H
#define USETHEFORCE_DARKSIDE_MPISIMULATION_INTERFACE_H

#define __C
#include "config.h"

extern "C" {
  void F90_FUNC(setfortranmpi, SETFORTRANMPI)( mpiSimData* the_mpiPlug,
                                               int* nLocal, 
                                               int* globalAtomIndex,
                                               int* nGroupsLocal,
                                               int* globalGroupIndex,
                                               int* isError );
  
  void setFortranMPI( mpiSimData* the_mpiPlug,
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
#endif
