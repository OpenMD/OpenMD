/*
 * Copyright (c) 2007 The University of Notre Dame. All Rights Reserved.
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
 
#include "io/NonBondedInteractionsSectionParser.hpp"
#include "types/AtomType.hpp"
#include "types/ShiftedMorseInteractionType.hpp"
#include "types/MAWInteractionType.hpp"
#include "types/LennardJonesInteractionType.hpp"
#include "types/RepulsiveMorseInteractionType.hpp"
#include "types/RepulsivePowerInteractionType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/simError.h"
namespace OpenMD {

  NonBondedInteractionsSectionParser::NonBondedInteractionsSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("NonBondedInteractions");
    
    stringToEnumMap_["MAW"] =  MAW;                
    stringToEnumMap_["ShiftedMorse"] =  ShiftedMorse;
    stringToEnumMap_["LennardJones"] = LennardJones;
    stringToEnumMap_["RepulsiveMorse"] = RepulsiveMorse;
    stringToEnumMap_["RepulsivePower"] = RepulsivePower;
    
  }
  
  void NonBondedInteractionsSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    NonBondedInteractionType* nbiType = NULL;
    int nTokens = tokenizer.countTokens();
    
    if (nTokens < 3) {
      sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    }
    
    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string itype = tokenizer.nextToken();
    
    NonBondedInteractionTypeEnum nbit = getNonBondedInteractionTypeEnum(itype);
    nTokens -= 3;
    NonBondedInteractionType* interactionType;
    
    //switch is a nightmare to maintain
    switch(nbit) {
    case MAW :
      if (nTokens < 5) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r_e = tokenizer.nextTokenAsDouble();
        RealType D_e = tokenizer.nextTokenAsDouble();
        RealType beta = tokenizer.nextTokenAsDouble();
        RealType ca1 = tokenizer.nextTokenAsDouble();
        RealType cb1 = tokenizer.nextTokenAsDouble();
        interactionType = new MAWInteractionType(D_e, beta, r_e, ca1, cb1);
      }
      break;
      
    case ShiftedMorse :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0 = tokenizer.nextTokenAsDouble();
        RealType D0 = tokenizer.nextTokenAsDouble();
        RealType beta0 = tokenizer.nextTokenAsDouble();
        interactionType = new ShiftedMorseInteractionType(D0, beta0, r0);
      }
      break;
      
    case RepulsiveMorse :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0 = tokenizer.nextTokenAsDouble();
        RealType D0 = tokenizer.nextTokenAsDouble();
        RealType beta0 = tokenizer.nextTokenAsDouble();
        interactionType = new RepulsiveMorseInteractionType(D0, beta0, r0);
      }
      break;
      
    case LennardJones :
      if (nTokens < 2) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma = tokenizer.nextTokenAsDouble();
        RealType epsilon = tokenizer.nextTokenAsDouble();
        interactionType = new LennardJonesInteractionType(sigma, epsilon);
      }
      break;

    case RepulsivePower :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma = tokenizer.nextTokenAsDouble();
        RealType epsilon = tokenizer.nextTokenAsDouble();
        int nRep = tokenizer.nextTokenAsInt();
        interactionType = new RepulsivePowerInteractionType(sigma, epsilon, nRep);
      }
      break;
      
    case Unknown :
    default:
      sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Unknown Interaction Type at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
      
      break;
            
    }
    
    if (interactionType != NULL) {
      ff.addNonBondedInteractionType(at1, at2, interactionType);
    }
    
  }
  
  NonBondedInteractionsSectionParser::NonBondedInteractionTypeEnum NonBondedInteractionsSectionParser::getNonBondedInteractionTypeEnum(const std::string& str) {
    std::map<std::string, NonBondedInteractionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);
    
    return i == stringToEnumMap_.end() ? Unknown : i->second;
  }
  
} //end namespace OpenMD

