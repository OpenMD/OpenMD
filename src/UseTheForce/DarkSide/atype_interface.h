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
#include "types/AtomTypeProperties.h"

#define makeAtype F90_FUNC(makeatype, MAKEATYPE)

extern "C" {
  void makeAtype(AtomTypeProperties atp, int* status );
}  
#endif
