/**
 * @file shapes_module.h
 * @author Dan Gezelter
 * @date 10/19/2004
 * @version 1.0
 */

#ifndef USETHEFORCE_DARKSIDE_SHAPES_MODULE_H
#define USETHEFORCE_DARKSIDE_SHAPES_MODULE_H

#define __C

#include "config.h"

#define SH_COS 0
#define SH_SIN 1

extern "C" {

  void F90_FUNC(makeshape, MAKESHAPE)(int* nContactFuncs, 
                                      int* ContactFuncLValue, 
                                      int* ContactFuncMValue, 
                                      int* ContactFunctionType, 
                                      double* ContactFuncCoefficient, 
                                      int* nRangeFuncs, 
                                      int* RangeFuncLValue,
                                      int* RangeFuncMValue, 
                                      int* RangeFunctionType, 
                                      double* RangeFuncCoefficient, 
                                      int* nStrengthFuncs, 
                                      int* StrengthFuncLValue, 
                                      int* StrengthFuncMValue, 
                                      int* StrengthFunctionType, 
                                      double* StrengthFuncCoefficient,
                                      int* myAtid, 
                                      int* status);

  
  void makeShape(int* nContactFuncs, 
                 int* ContactFuncLValue, 
                 int* ContactFuncMValue, 
                 int* ContactFunctionType, 
                 double* ContactFuncCoefficient, 
                 int* nRangeFuncs, 
                 int* RangeFuncLValue,
                 int* RangeFuncMValue, 
                 int* RangeFunctionType, 
                 double* RangeFuncCoefficient, 
                 int* nStrengthFuncs, 
                 int* StrengthFuncLValue, 
                 int* StrengthFuncMValue, 
                 int* StrengthFunctionType, 
                 double* StrengthFuncCoefficient,
                 int* myAtid, 
                 int* status) {
  
    F90_FUNC(makeshape, MAKESHAPE)( nContactFuncs, 
                                    ContactFuncLValue, 
                                    ContactFuncMValue, 
                                    ContactFunctionType, 
                                    ContactFuncCoefficient, 
                                    nRangeFuncs, 
                                    RangeFuncLValue,
                                    RangeFuncMValue, 
                                    RangeFunctionType, 
                                    RangeFuncCoefficient, 
                                    nStrengthFuncs, 
                                    StrengthFuncLValue, 
                                    StrengthFuncMValue, 
                                    StrengthFunctionType, 
                                    StrengthFuncCoefficient,
                                    myAtid, 
                                    status);
  }
}
    
#endif
