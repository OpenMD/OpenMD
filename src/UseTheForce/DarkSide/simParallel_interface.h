/*
 *  simParallel_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef USETHEFORCE_DARKSIDE_SIMPARALLEL_INTERFACE_H
#define USETHEFORCE_DARKSIDE_SIMPARALLEL_INTERFACE_H

#define __C
#include "config.h"

#define setFsimParallel F90_FUNC(setfsimparallel, SETFSIMPARALLEL)

extern "C" {
  void setFsimParallel( mpiSimData* the_mpiPlug,
                        int* nLocal, 
                        int* globalAtomIndex,
                        int* nGroupsLocal,
                        int* globalGroupIndex,
                        int* isError );
}
#endif
