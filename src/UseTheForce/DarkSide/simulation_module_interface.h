/*
 *  simulation_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */


#define __C
#include "config.h"
extern "C"{
  typedef void (setFortranSim)( simtype* the_Info,
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
                                int* isError ){
    F90_FUNC(setfortransim, SETFORTRANSIM)( nGlobal,
                                            nLocal,
                                            identArray,
                                            nLocalExcludes,
                                            excludesLocalArray,
                                            nGlobalExcludes,
                                            excludesGlobalArray,
                                            molMembershipArray,
                                            mfact,
                                            ngroup,
                                            globalGroupMembership,
                                            isError);
  }
  typedef void (setFortranBox) ( double *Hmat,
                                 double *HmatI,
                                 int* orthoRhombic ){
    F90_FUNC(setfortranbox,SETFORTRANBOX)( Hmat,
                                           HmatI,
                                           orthoRhombic)
    
  }
  
}