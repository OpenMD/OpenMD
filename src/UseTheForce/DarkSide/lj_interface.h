/*
 *  lj_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/21/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */




#ifndef USETHEFORCE_DARKSIDE_LJ_INTERFACE_H
#define USETHEFORCE_DARKSIDE_LJ_INTERFACE_H

#define __C

#include "config.h"

#define newLJtype F90_FUNC(newljtype, NEWLJTYPE)
#define useGeometricMixing F90_FUNC(usegeometricmixing, USEGEOMETRICMIXING)

extern "C"{
  void newLJtype( int* ident,
                  double* sigma,
                  double* epsilon
                  int* status);

  void useGeometricMixing( );
}  
#endif

