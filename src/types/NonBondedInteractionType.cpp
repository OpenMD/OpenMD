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
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "types/NonBondedInteractionType.hpp"
#include "utils/simError.h"

namespace OpenMD {

  void NonBondedInteractionType::setAtomTypes(std::pair<AtomType*, AtomType*> ats) { 
    atomTypes_=ats; 
  }

  std::pair<AtomType*, AtomType*> NonBondedInteractionType::getAtomTypes() {
    return atomTypes_;
  }
   
  void NonBondedInteractionType::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
  }
  
  void NonBondedInteractionType::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);  
  }
  
  void NonBondedInteractionType::clearProperties() {
    properties_.clearProperties();
  }
  
  std::vector<std::string> NonBondedInteractionType::getPropertyNames() {
    return properties_.getPropertyNames();  
  }
  
  std::vector<GenericData*> NonBondedInteractionType::getProperties() { 
    return properties_.getProperties(); 
  }
  
  GenericData* NonBondedInteractionType::getPropertyByName(const std::string& propName) {
      return properties_.getPropertyByName(propName); 
  }  
  
  void NonBondedInteractionType::setLennardJones() {
    nbitp.is_LennardJones = 1;
  }
  
  bool NonBondedInteractionType::isLennardJones() {
    return nbitp.is_LennardJones;
  }

  void NonBondedInteractionType::setMorse() {
    nbitp.is_Morse = 1;
  }
  
  bool NonBondedInteractionType::isMorse() {
    return nbitp.is_Morse;
  }
  
  void NonBondedInteractionType::setEAM() {
    nbitp.is_EAM = 1;
  }
  
  bool NonBondedInteractionType::isEAM() {
    return nbitp.is_EAM;
  }
  
  bool NonBondedInteractionType::isSC() {
    return nbitp.is_SC;
  }
  
  void NonBondedInteractionType::setSC() {
    nbitp.is_SC = 1;
  }
  
  bool NonBondedInteractionType::isMetal() {
    return isSC() || isEAM();
  }  
}
