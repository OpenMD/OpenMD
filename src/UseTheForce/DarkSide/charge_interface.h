/*
 *  charge_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/21/04.
 *  Copyright 2004 University of Notre Dame. All rights reserved.
 *
 */




#ifndef USETHEFORCE_DARKSIDE_CHARGE_INTERFACE_H
#define USETHEFORCE_DARKSIDE_CHARGE_INTERFACE_H

#define __C

#include "config.h"

#define newChargeType F90_FUNC(newchargetype, NEWCHARGETYPE)

extern "C"{
  void newChargeType( int* ident,
                      double* charge,
                      int* status);
}  
#endif

