/*
 *  notifycutoffs_module_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#define __C
#include "config.h"
extern "C"{
  typedef void (notifyFortranCutOff) ( double *rCut,
                                       double *rSw,
                                       double *rList ){
    F90_FUNC(notifyfortrancutoff,NOTIFYFORTRANCUTOFF)( rCut,
                                                       rSw,
                                                       rList)
  }
  
}