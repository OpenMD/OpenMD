/*
 *  dipole_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/21/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */




#ifndef USETHEFORCE_DARKSIDE_DIPOLE_INTERFACE_H
#define USETHEFORCE_DARKSIDE_DIPOLE_INTERFACE_H

#define __C

#include "config.h"

#define newDipoleType F90_FUNC(newdipoletype, NEWDIPOLETYPE)

extern "C"{
  void newDipoleType( int* ident,
                      double* dipole_moment,
                      int* status);
}  
#endif

