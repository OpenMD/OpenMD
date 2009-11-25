/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
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
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#ifndef USETHEFORCE_DARKSIDE_ELECTROSTATIC_INTERFACE_H
#define USETHEFORCE_DARKSIDE_ELECTROSTATIC_INTERFACE_H

#define __OPENMD_C

#include "config.h"
#include "types/AtomTypeProperties.h"
 
#define setElectrostaticSummationMethod FC_FUNC(setelectrostaticsummationmethod, SETELECTROSTATICSUMMATIONMETHOD)
#define setScreeningMethod FC_FUNC(setscreeningmethod, SETSCREENINGMETHOD)
#define setElectrostaticCutoffRadius FC_FUNC(setelectrostaticcutoffradius, SETELECTROSTATICCUTOFFRADIUS)
#define setDampingAlpha FC_FUNC(setdampingalpha, SETDAMPINGALPHA)
#define setReactionFieldDielectric FC_FUNC(setreactionfielddielectric, SETREACTIONFIELDDIELECTRIC)

#define newElectrostaticType FC_FUNC(newelectrostatictype, NEWELECTROSTATICTYPE)
#define setCharge FC_FUNC(setcharge, SETCHARGE)
#define setDipoleMoment FC_FUNC(setdipolemoment, SETDIPOLEMOMENT)
#define setSplitDipoleDistance FC_FUNC(setsplitdipoledistance, SETSPLITDIPOLEDISTANCE)

#define setQuadrupoleMoments FC_FUNC(setquadrupolemoments, SETQUADRUPOLEMOMENTS)

#define destroyElectrostaticTypes FC_FUNC(destroyelectrostatictypes,DESTROYELECTROSTATICTYPES)
extern "C"{

  void setElectrostaticSummationMethod( int* theESM );
  void setScreeningMethod( int* theSM );
  void setElectrostaticCutoffRadius( RealType* theECR, RealType* theRSW );
  void setDampingAlpha( RealType* theDA );
  void setReactionFieldDielectric( RealType* theDielectric );

  void newElectrostaticType( AtomTypeProperties* atp,
                             int* status);
  
  void setCharge( int* c_ident, 
                  RealType* charge,
                  int* status);
  
  void setDipoleMoment( int* c_ident,
                        RealType* dipole_moment,
                        int* status);
  
  void setSplitDipoleDistance( int* c_ident,
                               RealType* split_dipole_distance,
                               int* status);
  
  void setQuadrupoleMoments( int* c_ident,
                             RealType* quadrupole_moments,
                             int* status);
	
  void destroyElectrostaticTypes(void);
}  
#endif

