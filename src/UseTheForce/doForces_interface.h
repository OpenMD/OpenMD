/*
 *  doForces_interface.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef USETHEFORCE_DOFORCES_INTERFACE_H
#define USETHEFORCE_DOFORCES_INTERFACE_H

#define __C
#include "config.h"

#define initFortranFF F90_FUNC(initfortranff, INITFORTRANFF)
#define doForceLoop F90_FUNC(doforceloop, DOFORCELOOP)

extern "C"{
  
  void initFortranFF( int* useReactionField,
                      int *isError );        

  
  void doForceLoop( double* positionArray,
                    double* rcArray,
                    double* RotationMatrixArray,
                    double* unitVectorArray_l,
                    double* forceArray,
                    double *torqueArray,
                    double* StressTensor, 
                    double* potentialEnergy, 
                    short int* doPotentialCalc, 
                    short int* doStressCalc,
                    int* isError );
}
#endif
