/*
 *  simulation_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */

#ifndef USETHEFORCE_DARKSIDE_SIMULATION_INTERFACE_H
#define USETHEFORCE_DARKSIDE_SIMULATION_INTERFACE_H

#define __C
#include "brains/fSimulation.h"
#include "config.h"

#define setFortranSim F90_FUNC(setfortransim, SETFORTRANSIM)
#define setFortranBox F90_FUNC(setfortranbox, SETFORTRANBOX)

extern "C"{
  void setFortranSim( simtype* the_Info,
                      int* nGlobal, 
                      int* nLocal, 
                      int* identArray,
                      int* nLocalExcludes,
                      int* excludesLocalArray,
                      int* nGlobalExcludes,
                      int* excludesGlobalArray,
                      int* molMembershipArray,
                      double* mfact,
                      int* ngroup,
                      int* globalGroupMembership,
                      int* isError );
  
  void setFortranBox( double *Hmat,
                      double *HmatI,
                      int* orthoRhombic );
}
#endif
