/*
 *  atype_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef USETHEFORCE_DARKSIDE_ATYPE_INTERFACE_H
#define USETHEFORCE_DARKSIDE_ATYPE_INTERFACE_H

#define __C

#include "config.h"

#define makeAtype F90_FUNC(makeatype, MAKEATYPE)

extern "C" {
  void makeAtype(int* unique_ident, 
                 int* isLJ, 
                 int* isSticky,
                 int* isDipole, 
                 int* isGB, 
                 int* isEAM,
                 int* isCharge,
                 double* lj_epslon, 
                 double* lj_sigma, 
                 double* charge,
                 double* dipole_moment, 
                 int* status );
}  
#endif
