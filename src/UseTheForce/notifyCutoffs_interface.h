/*
 *  notifycutoffs_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef USETHEFORCE_NOTIFYCUTOFFS_INTERFACE_H
#define USETHEFORCE_NOTIFYCUTOFFS_INTERFACE_H

#define __C
#include "config.h"

#define notifyFortranCutoffs F90_FUNC(notifyfortrancutoffs, NOTIFYFORTRANCUTOFFS)

extern "C"{
  
 void notifyFortranCutoffs( double *rCut,
                           double *rSw,
                           double *rList );
  
}
#endif
