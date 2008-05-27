/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
#ifndef USETHEFORCE_DOFORCES_INTERFACE_H
#define USETHEFORCE_DOFORCES_INTERFACE_H

#define __C
#include "config.h"

#define initFortranFF FC_FUNC(initfortranff, INITFORTRANFF)
#define doForceLoop FC_FUNC(doforceloop, DOFORCELOOP)
#define getAccumulatedBoxDipole FC_FUNC(getaccumulatedboxdipole, GETACCUMULATEDBOXDIPOLE)
#define setAccumulateBoxDipole FC_FUNC(setaccumulateboxdipole, SETACCUMULATEBOXDIPOLE)
#define setFortranElectrostaticMethod FC_FUNC(setfortranelectrostaticmethod, SETFORTRANELECTROSTATICMETHOD)
#define notifyFortranCutoffPolicy FC_FUNC(notifyfortrancutoffpolicy, NOTIFYFORTRANCUTOFFPOLICY)
#define notifyFortranSkinThickness FC_FUNC(notifyfortranskinthickness, NOTIFYFORTRANSKINTHICKNESS)
#define notifyFortranCutoffs FC_FUNC(notifyfortrancutoffs, NOTIFYFORTRANCUTOFFS)
#define notifyFortranYouAreOnYourOwn FC_FUNC(notifyfortranyouareonyourown, NOTIFYFORTRANYOUAREONYOUROWN)

extern "C"{
  
  void initFortranFF( int* isError );        
 
  void doForceLoop( RealType* positionArray,
                    RealType* rcArray,
                    RealType* RotationMatrixArray,
                    RealType* unitVectorArray_l,
                    RealType* forceArray,
                    RealType *torqueArray,
                    RealType* StressTensor, 
                    RealType* potentialEnergy,
		    RealType* particlePotArray,
                    short int* doPotentialCalc, 
                    short int* doStressCalc,
                    int* isError );

  void getAccumulatedBoxDipole( RealType* boxDipole );

  void setAccumulateBoxDipole();

  void setFortranElectrostaticMethod( int* electrostaticMethod );

  void notifyFortranCutoffPolicy( int* cutPolicy );

  void notifyFortranSkinThickness( RealType *skinThickness );

  void notifyFortranCutoffs( RealType *rCut,
			     RealType *rSw,
			     bool *ljsp,
			     bool *ljsf);

  void notifyFortranYouAreOnYourOwn( );

}
#endif
