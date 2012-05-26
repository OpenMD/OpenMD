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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "io/InversionTypesSectionParser.hpp"
#include "types/InversionType.hpp"
#include "types/ImproperCosineInversionType.hpp"
//#include "types/ImproperHarmonicInversionType.hpp"
//#include "types/CentralAtomHeightInversionType.hpp"
//#include "types/DreidingInversionType.hpp"
#include "types/AmberImproperTorsionType.hpp"
#include "types/PolynomialInversionType.hpp"
//These two are added by me. Maybe it is wrong.
#include "brains/ForceField.hpp"
#include "utils/NumericConstant.hpp"

namespace OpenMD {
  
  InversionTypesSectionParser::InversionTypesSectionParser(ForceFieldOptions& options) : options_(options){
    
    setSectionName("InversionTypes");
    stringToEnumMap_["AmberImproper"] = itAmberImproper;
    stringToEnumMap_["ImproperCosine"] = itImproperCosine;
    stringToEnumMap_["ImproperHarmonic"] = itImproperHarmonic;
    stringToEnumMap_["CentralAtomHeight"] = itCentralAtomHeight;
    stringToEnumMap_["Dreiding"] = itDreiding;
  }
  
  void InversionTypesSectionParser::parseLine(ForceField& ff,
					      const std::string& line, 
					      int lineNo){
    StringTokenizer tokenizer(line);
    InversionType* inversionType = NULL;
    
    int nTokens = tokenizer.countTokens();
    
    if (nTokens < 5) {
      return;
    }
    
    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string at3 = tokenizer.nextToken();
    std::string at4 = tokenizer.nextToken();
    InversionTypeEnum it = getInversionTypeEnum(tokenizer.nextToken());

    nTokens -= 5;

    switch(it) {
      
    case InversionTypesSectionParser::itImproperCosine :

      if (nTokens < 3 || nTokens % 3 != 0) {

      } else {
        int nSets = nTokens / 3;
  
        std::vector<ImproperCosineInversionParameter> parameters;             
        for (int i = 0; i < nSets; ++i) {
          ImproperCosineInversionParameter currParam;
          currParam.kchi = tokenizer.nextTokenAsDouble();
          currParam.n = tokenizer.nextTokenAsInt();
	  currParam.delta = tokenizer.nextTokenAsDouble() / 180.0 * NumericConstant::PI; //convert to rad
          parameters.push_back(currParam);
        }

        inversionType = new ImproperCosineInversionType(parameters);

      }
    break;

    case InversionTypesSectionParser::itAmberImproper:
            
      if (nTokens < 1) {

      } else {
	RealType v2 = tokenizer.nextTokenAsDouble();
        
	inversionType = new AmberImproperTorsionType(v2);
      }

      break;
            



      /*      
    case InversionTypesSectionParser::itImproperHarmonic :
      if (nTokens < 2) {
	
      } else {
	
	RealType k = tokenizer.nextTokenAsDouble();
        
	inversionType = new ImproperHarmonicInversionType(k);
      }
      break;
       
    case InversionTypesSectionParser::itCentralAtomHeight :
      if (nTokens < 1) {
	
      } else {
	
	RealType k = tokenizer.nextTokenAsDouble();
        
	inversionType = new CentralAtomHeightInversionType(k);
      }
      break;                 
    case InversionTypesSectionParser::itDreiding :
      if (nTokens < 3) {
	
      } else {
	
	RealType k = tokenizer.nextTokenAsDouble();
        
	inversionType = new CentralAtomHeightInversionType(k);
      }
      break;                 
      */
    case InversionTypesSectionParser::itUnknown :
    default:
      
      break;
      
    }
    
    if (inversionType != NULL) {
      ff.addInversionType(at1, at2, at3, at4, inversionType);
    }

  }
  
  InversionTypesSectionParser::InversionTypeEnum InversionTypesSectionParser::getInversionTypeEnum(const std::string& str) {
    std::map<std::string, InversionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);
    
    return i == stringToEnumMap_.end() ? InversionTypesSectionParser::itUnknown : i->second;
  }
  
} //end namespace OpenMD



