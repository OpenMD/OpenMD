/*
 *  do_Forces.h
 *  oopse
 *
 *  Created by Charles Vardeman II on 10/19/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef DO_FORCES_H
#define DO_FORCES_H

#define __C
#include "config.h"
extern "C"{
  
  void F90_FUNC(initFortranff,INITFORTRANFF)( int* LJ_mix_policy, 
                                              int* useReactionField,
                                              int *isError );        
  void (initFortranFF)( int* LJ_mix_policy, 
                        int* useReactionField,
                        int *isError ){           
    F90_FUNC(initFortranff,INITFORTRANFF)( LJ_mix_policy, 
                                           useReactionField,
                                           isError 
                                           );
  }
  
  void F90_FUNC(doforceloop,DOFORCELOOP)( double* positionArray,
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
  
  void (doForceLoop)( double* positionArray,
                      double* rcArray,
                      double* RotationMatrixArray,
                      double* unitVectorArray_l,
                      double* forceArray,
                      double *torqueArray,
                      double* StressTensor, 
                      double* potentialEnergy, 
                      short int* doPotentialCalc, 
                      short int* doStressCalc,
                      int* isError ){
    F90_FUNC(doforceloop,DOFORCELOOP)( positionArray,
                                       rcArray,
                                       RotationMatrixArray,
                                       unitVectorArray_l,
                                       forceArray,
                                       torqueArray,
                                       StressTensor, 
                                       potentialEnergy, 
                                       doPotentialCalc, 
                                       doStressCalc,
                                       isError );
  }
  
  
}

#endif