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
 
#include "io/TorsionTypesSectionParser.hpp"
#include "types/TorsionType.hpp"
#include "types/CubicTorsionType.hpp"
#include "types/QuarticTorsionType.hpp"
#include "types/PolynomialTorsionType.hpp"
#include "types/CharmmTorsionType.hpp"
#include "types/OplsTorsionType.hpp"
#include "types/TrappeTorsionType.hpp"
#include "UseTheForce/ForceField.hpp"
#include "utils/NumericConstant.hpp"
//#include <fstream>

namespace OpenMD {

  TorsionTypesSectionParser::TorsionTypesSectionParser(ForceFieldOptions& options) : options_(options){

    setSectionName("TorsionTypes");
    stringToEnumMap_["GhostTorsion"] = ttGhostTorsion;
    stringToEnumMap_["Cubic"] = ttCubic;
    stringToEnumMap_["Quartic"] = ttQuartic;
    stringToEnumMap_["Polynomial"] = ttPolynomial;
    stringToEnumMap_["Charmm"] =  ttCharmm;
    stringToEnumMap_["Opls"] =  ttOpls;
    stringToEnumMap_["Trappe"] =  ttTrappe;
  }

  void TorsionTypesSectionParser::parseLine(ForceField& ff,
					    const std::string& line, 
					    int lineNo){
    StringTokenizer tokenizer(line);
    TorsionType* torsionType = NULL;
    
    int nTokens = tokenizer.countTokens();

    if (nTokens < 5) {
      return;
    }
    
    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string at3 = tokenizer.nextToken();
    std::string at4 = tokenizer.nextToken();
    TorsionTypeEnum tt = getTorsionTypeEnum(tokenizer.nextToken());

    //std::cout << at1 << "-" << at2 << "-" << at3 << "-" << at4 << std::endl;
    nTokens -= 5;

    switch(tt) {

    case TorsionTypesSectionParser::ttGhostTorsion:
      if (nTokens < 4) {

      } else {

	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
	torsionType = new CubicTorsionType(k3, k2, k1, k0);
      }
      break;
            
    case TorsionTypesSectionParser::ttCubic :
      if (nTokens < 4) {

      } else {

	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
	torsionType = new CubicTorsionType(k3, k2, k1, k0);
      }
      break;
            
    case TorsionTypesSectionParser::ttQuartic:
      if (nTokens < 5) {

      } else {

	RealType k4 = tokenizer.nextTokenAsDouble();
	RealType k3 = tokenizer.nextTokenAsDouble();
	RealType k2 = tokenizer.nextTokenAsDouble();
	RealType k1 = tokenizer.nextTokenAsDouble();
	RealType k0 = tokenizer.nextTokenAsDouble();
                
	torsionType = new QuarticTorsionType( k4, k3, k2, k1, k0);
      }
      break;

        
    case TorsionTypesSectionParser::ttPolynomial:
      if (nTokens < 2 || nTokens % 2 != 0) {

      } else {
	int nPairs = nTokens / 2;
	int power;
	RealType coefficient;
	torsionType = new PolynomialTorsionType();

        PolynomialTorsionType* ptt = dynamic_cast<PolynomialTorsionType*>(torsionType); 
                
	for (int i = 0; i < nPairs; ++i) {
	  power = tokenizer.nextTokenAsInt();
	  coefficient = tokenizer.nextTokenAsDouble();
	  ptt->setCoefficient(power, coefficient);
	}
      }
            
      break;
             
    case TorsionTypesSectionParser::ttCharmm:
            
      if (nTokens < 3 || nTokens % 3 != 0) {

      } else {
	int nSets = nTokens / 3;
  
        std::vector<CharmmTorsionParameter> parameters;             
	for (int i = 0; i < nSets; ++i) {
          CharmmTorsionParameter currParam;
	  currParam.kchi = tokenizer.nextTokenAsDouble();
	  currParam.n = tokenizer.nextTokenAsInt();
	  currParam.delta = tokenizer.nextTokenAsDouble() / 180.0 * NumericConstant::PI; //convert to rad
          parameters.push_back(currParam);
	}

	torsionType = new CharmmTorsionType(parameters);

      }

      break;

    case TorsionTypesSectionParser::ttOpls:
            
      if (nTokens < 3) {

      } else {
	RealType v1 = tokenizer.nextTokenAsDouble();
	RealType v2 = tokenizer.nextTokenAsDouble();
	RealType v3 = tokenizer.nextTokenAsDouble();
        
	torsionType = new OplsTorsionType(v1, v2, v3);
      }

      break;
            

    case TorsionTypesSectionParser::ttTrappe:
            
      if (nTokens < 4) {

      } else {

	RealType c0 = tokenizer.nextTokenAsDouble();
	RealType c1 = tokenizer.nextTokenAsDouble();
	RealType c2 = tokenizer.nextTokenAsDouble();
	RealType c3 = tokenizer.nextTokenAsDouble();
	torsionType = new TrappeTorsionType(c0, c1, c2, c3);
      }

      break;
            
    case TorsionTypesSectionParser::ttUnknown :
    default:

      break;
            
    }

    if (torsionType != NULL) {
      ff.addTorsionType(at1, at2, at3, at4, torsionType);
    }

  }

  TorsionTypesSectionParser::TorsionTypeEnum TorsionTypesSectionParser::getTorsionTypeEnum(const std::string& str) {
    std::map<std::string, TorsionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);
    
    return i == stringToEnumMap_.end() ? TorsionTypesSectionParser::ttUnknown : i->second;
  }
  
} //end namespace OpenMD



